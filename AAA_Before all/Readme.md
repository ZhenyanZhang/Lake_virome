# Raw Data Processing and Metagenome-Assembled Genomes (MAGs) Recovery

Before viral identification, raw sequencing data was processed through a standardized pipeline including quality control, de novo assembly, contig binning, and genome dereplication to obtain high-quality Metagenome-Assembled Genomes (MAGs).

## 1. Raw Data Quality Control

Raw paired-end sequencing reads were processed using `fastp` to remove adapters, filter low-quality bases, and trim reads. 

Execute the custom parallel script ([run_fastp.sh](./run_fastp.sh)) specifying the input directory, output directory, and the maximum number of parallel jobs[cite: 17]:

```bash
mkdir -p 1.cleandata
bash run_fastp.sh 0.rawdata 1.cleandata 5
```

## 2. Metagenome Assembly

Clean reads from each sample were individually assembled into contigs using `MEGAHIT` (v1.2.9)([run_megahit.sh](./run_megahit.sh)):

```bash
mkdir -p 2.assembly
bash run_megahit.sh 1.cleandata 2.temp_assembly 3

# Extract and pool all generated contig fasta files into a unified directory for downstream analyses
mkdir -p 2.contigs
find 2.assembly -name "*.contigs.fa" -exec cp {} 2.contigs/ \;
```

## 3. Contig Binning

Contigs were mapped back to the clean reads using `BWA-MEM` to calculate coverage depths. Subsequently, `MetaBAT2` (v2.12.1) was used to cluster the contigs into bins representing putative MAGs.

Execute the binning script ([run_metabat2.sh](./run_metabat2.sh)):

```bash
mkdir -p 3.MAG
bash run_metabat2.sh 2.contigs 1.cleandata 3.MAG/depth 3.MAG/bins 3

# Consolidate all generated bins from individual sample folders into a single directory for QC
mkdir -p 3.MAG/all_bins
find 3.MAG/bins -name "*.fa" -exec cp {} 3.MAG/all_bins/ \;
```

## 4. MAG Quality Control and Filtering

To assess the quality of the recovered genomes, `CheckM` (v1.1.3) was utilized to estimate the completeness and contamination of each bin based on single-copy marker genes. Following the filtering strategy described by Nayfach et al. (DOI: [10.1038/s41586-024-07891-2](https://doi.org/10.1038/s41586-024-07891-2)), a comprehensive quality score was calculated for each MAG (Quality score = Completeness - 5 × Contamination). 

Only MAGs meeting the medium or high-quality criteria (Completeness ≥ 50%, Contamination ≤ 10%, and Quality score ≥ 50) were extracted and moved into a new directory for downstream analysis.

```bash
# Run CheckM on all collected bins and output results as a tab-separated table
checkm lineage_wf -t 64 -x fa 3.MAG/all_bins 3.MAG/checkm_out --tab_table -f 3.MAG/checkm_out/checkm_results.tsv

# Create a directory for the filtered high-quality MAGs
mkdir -p 3.MAG/filtered_bins

# Use AWK to calculate the quality score and filter the MAGs 
# (Note: Col 12 is typically Completeness, Col 13 is Contamination in CheckM tab output)
awk -F'\t' 'NR>1 {
    bin_id = $1;
    comp = $12;
    cont = $13;
    qs = comp - (5 * cont);
    
    if (comp >= 50 && cont <= 10 && qs >= 50) {
        print bin_id
    }
}' 3.MAG/checkm_out/checkm_results.tsv > 3.MAG/checkm_out/passed_bins.list

# Copy the MAGs that passed the quality filter into the new folder
for bin in $(cat 3.MAG/checkm_out/passed_bins.list); do
    cp 3.MAG/all_bins/${bin}.fa 3.MAG/filtered_bins/
done
```

## 5. Dereplication

To obtain a non-redundant catalog of host genomes, the pre-filtered MAGs were dereplicated using `dRep` (v3.4.5) with a 95% average nucleotide identity (ANI) threshold. Since quality filtering was already strictly performed in the previous step, the internal genome quality assessment in `dRep` was bypassed.

```bash
# Dereplicate genomes using dRep (skipping internal quality checks)
dRep dereplicate 3.MAG/bin_drep \
    -g 3.MAG/filtered_bins/*.fa \
    -p 64 \
    -sa 0.95 \
    --ignoreGenomeQuality
```
