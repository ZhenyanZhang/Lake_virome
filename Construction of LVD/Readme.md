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
cd-hit-est -i all.phage.scaffold.MH.new.fna -o vOTUs_new.fna -c 0.95 -aS 0.85 -d 0 -T 100 -M 0
```
