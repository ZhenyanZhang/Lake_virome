# Taxonomic assignment of vOTUs

The taxonomic classification pipeline utilizes a consensus-based approach integrating multiple databases. The custom classification scripts used in this section were modified from the TYMEFLIES_Viral repository by the Anantharaman Lab (DOI: [10.1038/s41564-024-01876-7](https://doi.org/10.1038/s41564-024-01876-7), GitHub: [https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification)).

## 1. Prepare Protein-to-vOTU Map

First, create a necessary working directory and generate the mapping file that links each predicted viral protein to its parent vOTU contig sequence.

```bash
# Create the output directory for taxonomic classification
mkdir -p 6.Taxonomic_classification

# Extract protein ID and corresponding vOTU ID from the cleaned FASTA file
awk 'BEGIN{OFS="\t"} /^>/ {id=substr($1,2); votu=id; sub(/_[^_]+$/,"",votu); print id, votu}' \
    4.vibrant_results/gene_and_protein/vOTUs_clean.faa > 6.Taxonomic_classification/protein2genome_map.txt
```

## 2. Virus Taxonomic Annotation

Three independent taxonomic classification methods are executed in parallel:
1. DIAMOND BLASTp against the NCBI RefSeq viral protein database.
2. HMMsearch against the VOG (Virus Orthologous Groups) marker database.
3. Taxonomic annotation using geNomad.

```bash
# Method 1: NCBI RefSeq viral protein searching
nohup perl Taxonomic_classification.step1.run_diamond_to_NCBI_RefSeq_viral.pl \
    4.vibrant_results/gene_and_protein/vOTUs_clean.faa \
    6.Taxonomic_classification/protein2genome_map.txt \
    > 6.Taxonomic_classification/diamond_tax.out 2>&1 &

# Method 2: VOG marker HMM searching
nohup perl Taxonomic_classification.step2.run_hmmsearch_to_VOG_marker.pl \
    4.vibrant_results/gene_and_protein/vOTUs_clean.faa \
    6.Taxonomic_classification/protein2genome_map.txt \
    > 6.Taxonomic_classification/diamond_vog.out 2>&1 &

# Method 3: geNomad classifying
nohup genomad annotate 4.vibrant_results/vOTUs.fna \
    6.Taxonomic_classification/genomad_output db/viruse/genomad_db \
    -t 64 > 6.Taxonomic_classification/genomad_tax2.out 2>&1 &
```

*(Note: Wait for all three background `nohup` jobs to complete before proceeding to the next integration step).*

## 3. Integrate Annotation Results

To finalize the taxonomic classification, the integration script consolidates the outputs from the three methods and infers taxonomy for unclassified vOTUs based on their genus-level protein-sharing network clusters.

```bash
# Copy the genus cluster file to the current root directory as required by the script
cp 4.vibrant_results/gene_and_protein/genus_clusters.txt ./

# Run the integration script
perl Taxonomic_classification.step4.integrate_all_results.pl
```
