# Virus Identification

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
bash run_bowtie2.sh 1.cleandata 5.abundance/bam_output vOTUs_index sample_list.txt 3
```

### 6.3 Subsampling and CoverM Quantification

Filter out samples that do not meet the minimum sequencing depth (e.g., 50M reads). The remaining BAM files are subsampled to the target depth using `samtools`, quantified using `coverm contig`, and then filtered based on a specific breadth threshold (e.g., 0.75).

Run the pipeline using the custom filtering script ([run_subsample_filter.sh](./run_subsample_filter.sh)):

```bash
bash run_subsample_filter.sh 5.abundance/bam_output 0.75 5.abundance/final_Results
```

### 6.4 Merge Abundance Results

Combine all individual sample filtering results into a single comprehensive abundance table. The following command merges all output files while preserving only the first header row:

```bash
awk 'NR==1 || FNR>1' 5.abundance/final_Results/*.final_stats.txt > 5.abundance/merged_abundance_table.tsv
```
