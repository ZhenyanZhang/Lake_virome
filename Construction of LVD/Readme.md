# Construction of LVD

Prior to running VIBRANT, simplify the contig IDs by removing any characters after spaces in the FASTA headers.

```bash
# Simplify FASTA headers for all contig files in the input directory
for f in 2.contigs/*.fa; do
    sed -i 's/ .*//' "$f"
done
```

## 1. Run VIBRANT for Virus Identification

Execute the custom shell script ([run_vibrant.sh](./run_vibrant.sh)) to identify viral sequences:

```bash
bash run_vibrant.sh 2.contigs 4.vibrant_results
```

## 2. Quality Assessment

Run CheckV using the custom batch script ([run_checkV.sh](./run_checkV.sh)) to assess the quality of the identified viral contigs:

```bash
bash run_checkV.sh
```

## 3. Merge Sequences and Quality Reports

Consolidate all identified viral sequences and their corresponding quality summaries into single files:

```bash
# Merge viral FASTA sequences
find 4.vibrant_results -type f -path "*/CheckV_results/viruses.fna" -print0 | sort -z | xargs -0 cat > 4.vibrant_results/all.viruses.fna

# Merge CheckV quality summaries
find 4.vibrant_results -type f -path "*/CheckV_results/quality_summary.tsv" -print0 | sort -z | xargs -0 cat > 4.vibrant_results/all.quality_summary.tsv
```

## 4. Filter Contigs by Quality

Extract the contig IDs classified as Complete, High-quality, and Medium-quality, then filter the FASTA file accordingly using `seqtk`:

```bash
# Extract IDs for target quality tiers
awk -F'\t' '
    NR==1 {for(i=1;i<=NF;i++) if($i=="checkv_quality") col=i} 
    $col=="Complete" || $col=="High-quality" || $col=="Medium-quality" {print $1}
' 4.vibrant_results/all.quality_summary.tsv > 4.vibrant_results/filtered_id.list

# Filter the FASTA sequences based on the extracted ID list
seqtk subseq 4.vibrant_results/all.viruses.fna 4.vibrant_results/filtered_id.list > 4.vibrant_results/filtered_viruses.fna
```

## 5. Species-Level Clustering (vOTUs)

Consolidate the filtered viral contigs into species-level viral operational taxonomic units (vOTUs) using CD-HIT-EST:

```bash
cd-hit-est -i 4.vibrant_results/filtered_viruses.fna -o 4.vibrant_results/vOTUs.fna -c 0.95 -aS 0.85 -d 0 -T 100 -M 0
```

## 6. Read Mapping, Rarefaction, and Coverage Filtering

To accurately estimate viral abundances, reads are mapped back to the vOTUs index, followed by subsampling (rarefaction) for depth normalization and filtering based on breadth coverage thresholds.

### 6.1 Build Bowtie2 Index

Prior to read mapping, build a Bowtie2 index using the species-level vOTUs FASTA file generated from the previous clustering step:

```bash
bowtie2-build 4.vibrant_results/vOTUs.fna 5.abundance/vOTUs_index
```

### 6.2 Bowtie2 Alignment

Map the paired-end fastq reads to the established database index. Execute the custom parallel script ([run_bowtie2.sh](./run_bowtie2.sh)) specifying the inputs, database prefix, valid sample list, and allowed concurrent jobs:

```bash
bash run_bowtie2.sh 1.cleandata 5.abundance/bam_output 5.abundance/vOTUs_index sample_list.txt 3
```

### 6.3 Subsampling and CoverM Quantification

Filter out samples that do not meet the minimum sequencing depth (e.g., 50M reads). The remaining BAM files are subsampled to the target depth using `samtools`, quantified using `coverm contig`, and then filtered based on a specific breadth threshold (e.g., 0.75).

Run the pipeline using the custom filtering script ([run_subsample_filter.sh](./run_subsample_filter.sh)):

```bash
bash run_subsample_filter.sh 5.abundance/bam_output 0.75 5.abundance/final_results
```

### 6.4 Merge Abundance Results

Combine all individual sample filtering results into a single comprehensive abundance table. The following command merges all output files while preserving only the first header row:

```bash
awk 'NR==1 || FNR>1' 5.abundance/final_results/*.final_stats.txt > 5.abundance/5.abundance/all_abundance.txt
```

## 7. Rarefaction Analysis

To assess the coverage of lake viruses in the lake virome database, rarefaction analysis was performed using the `rtk` package in R (v4.3.1) with 100 random permutations.

The custom R script processes the combined abundance table into a wide matrix format, converts it to a presence/absence matrix, executes the permutation calculation, and plots the sample accumulation curve with a 95% confidence interval.

I executed this R script ([rarefaction_analysis.R](./rarefaction_analysis.R)) in R Studio in my Desktop (Windows 10).


## 8. Comparison of the lake virome database with other large viral databases

To determine the proportion of previously uncharacterized vOTUs in our lake virome database, we performed pairwise comparisons against major public viral databases, including the Global Oceans Viromes Database (v2.0), Gut Phage Database, Gut Virome Database, and IMG/VR v4. 

### 8.1 Build BLAST Databases

Before running the alignments, index the downloaded FASTA files into standard nucleotide BLAST databases. 

```bash
# Navigate to the db directory and index each database file
for db_file in db/*.fna; do
    makeblastdb -in "$db_file" -dbtype nucl -parse_seqids
done
```

### 8.2 Run Sequential BLASTn Alignments

Compare the newly identified vOTUs against each public database. The thresholds are strictly set to sequence similarity ≥ 90% and query coverage ≥ 75%. 

```bash
for db_file in db/*.fna; do
    db_name=$(basename -- "$db_file" .fna)
    blastn -query 4.vibrant_results/vOTUs.fna \
           -db "$db_file" \
           -out "4.vibrant_results/compare_other_db/blast_vs_${db_name}.tsv" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
           -perc_identity 90 \
           -qcov_hsp_perc 75 \
           -num_threads 64    
    echo "[Done] Finished ${db_name}"
done

echo "All pairwise comparisons are complete."
```

## 9. Gene Prediction and Taxonomic Clustering of vOTUs

To classify the vOTUs into genus-level and family-level taxonomic groups, we utilized a protein-sharing network approach based on the methodology described by Nayfach et al. (DOI: [10.1038/s41564-021-00928-6](https://doi.org/10.1038/s41564-021-00928-6)). The specific Python clustering scripts (`amino_acid_identity.py` and `filter_aai.py`) were adapted from the MGV repository ([https://github.com/snayfach/MGV/tree/master/aai_cluster](https://github.com/snayfach/MGV/tree/master/aai_cluster)).

### 9.1 ORF Prediction and Sequence Cleaning

Open reading frames (ORFs) were predicted using Prodigal. The resulting protein FASTA file was then cleaned to simplify sequence headers and remove stop codon asterisks (`*`) to ensure compatibility with downstream alignment tools.

```bash
# Predict proteins
prodigal -i 4.vibrant_results/vOTUs.fna -a 4.vibrant_results/gene_and_protein/vOTUs.faa -d 4.vibrant_results/gene_and_protein/vOTUs.genes.fna -p meta

# Clean FASTA headers and remove stop codons
awk '/^>/{print $1; next} {gsub(/\*/,""); print}' 4.vibrant_results/gene_and_protein/vOTUs.faa > 4.vibrant_results/gene_and_protein/vOTUs_clean.faa
```

### 9.2 All-vs-All Protein Alignment

An all-vs-all protein sequence alignment was performed using DIAMOND.

```bash
# Build DIAMOND database
diamond makedb --in 4.vibrant_results/gene_and_protein/vOTUs_clean.faa --db 4.vibrant_results/gene_and_protein/viral_proteins --threads 64

# Run BLASTp
diamond blastp --query 4.vibrant_results/gene_and_protein/vOTUs_clean.faa --db 4.vibrant_results/gene_and_protein/viral_proteins --out 4.vibrant_results/gene_and_protein/blastp.tsv --outfmt 6 --evalue 1e-5 --max-target-seqs 10000 --query-cover 50 --subject-cover 50 --threads 64
```

### 9.3 Protein-Sharing Network Clustering

Amino acid identity (AAI) and shared protein fractions between vOTUs were calculated based on the BLASTp results. The network edges were then filtered using specific sequence similarity and coverage thresholds to establish genus-level and family-level relationships. Finally, the Markov Cluster Algorithm (MCL) was utilized to generate the taxonomic clusters.

```bash
# Calculate AAI between genomes
python amino_acid_identity.py --in_faa 4.vibrant_results/gene_and_protein/vOTUs_clean.faa --in_blast 4.vibrant_results/gene_and_protein/blastp.tsv --out_tsv 4.vibrant_results/gene_and_protein/aai.tsv

# Filter edges for Genus-level clustering
python filter_aai.py --in_aai 4.vibrant_results/gene_and_protein/aai.tsv --min_percent_shared 20 --min_num_shared 16 --min_aai 40 --out_tsv 4.vibrant_results/gene_and_protein/genus_edges.tsv

# Filter edges for Family-level clustering
python filter_aai.py --in_aai 4.vibrant_results/gene_and_protein/aai.tsv --min_percent_shared 10 --min_num_shared 8 --min_aai 20 --out_tsv 4.vibrant_results/gene_and_protein/family_edges.tsv

# MCL Clustering
mcl 4.vibrant_results/gene_and_protein/genus_edges.tsv -te 8 -I 2.0 --abc -o 4.vibrant_results/gene_and_protein/genus_clusters.txt
mcl 4.vibrant_results/gene_and_protein/family_edges.tsv -te 8 -I 1.2 --abc -o 4.vibrant_results/gene_and_protein/family_clusters.txt
```
