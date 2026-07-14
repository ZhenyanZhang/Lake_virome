# Identifying AMGs carried by lake viruses

We annotated and identified AMGs in the vOTUs using the DRAM-v pipeline (v1.5.0), as recommended in previous study. Specifically, viral sequences were formatted using VirSorter2 (v2.2.4; `--prep-for-dramv`, `--min-length 1000`, `--include-groups dsDNAphage, ssDNA`) and functionally annotated with the default viral mode of DRAM-v (v1.5.0). According to the usage tutorial of DRAM-v, only viral genes containing the metabolic flag “M” were retained as putative AMGs. 

To ensure the reliability of the identified AMGs, these candidates were manually filtered as follows: 
(1) only genes flanked on both sides by viral or virus-like genes were retained, corresponding to DRAM-v auxiliary scores of 1 and 2; 
(2) genes carrying an “F” flag (located near the ends of scaffolds) were excluded; 
(3) genes likely involved in core viral processes or commonly misinterpreted as AMGs, including those associated with nucleotide metabolism, DNA modification, ribosomal functions, methyltransferases, glycosyltransferases, adenylyltransferases, and broadly acting carbohydrate-active enzymes, were excluded unless additional evidence supported a plausible auxiliary role. Functional categories of the retained AMGs were assigned by DRAM-v and manually verified against previously reported AMGs. The abundance of AMGs was calculated using the same method as for vOTUs.

## 1. Sequence Preparation with VirSorter2

First, format the viral sequences and prepare the specific inputs required by DRAM-v using VirSorter2. This step isolates the sequences and predicts the necessary viral features.

```bash
# Create the output directory
mkdir -p 8.AMG/vs2_out

# Run VirSorter2 to prepare inputs for DRAM-v
virsorter run \
    --use-conda-off \
    --db-dir ./db/vs2_db/db \
    --prep-for-dramv \
    -w 8.AMG/vs2_out \
    -i 4.vibrant_results/vOTUs.fna \
    --min-length 1000 \
    --include-groups dsDNAphage,ssDNA \
    -j 64
```

## 2. AMG Annotation with DRAM-v

Next, functionally annotate the processed viral contigs using the `DRAM-v.py annotate` module. This module assigns annotations to genes and determines their likelihood of being viral or cellular.

```bash
# Run DRAM-v annotation (skipping tRNAscan to speed up processing if tRNAs are not the focus)
DRAM-v.py annotate \
    -i 8.AMG/vs2_out/for-dramv/final-viral-combined-for-dramv.fa \
    -v 8.AMG/vs2_out/for-dramv/viral-affi-contigs-for-dramv.tab \
    -o 8.AMG/dramv_annotate \
    --threads 64 \
    --skip_trnascan
```

## 3. Distillation and Preliminary Filtering

Finally, summarize the annotations and perform the initial filtering using `DRAM-v.py distill`. Setting the `--max_auxiliary_score` to 2 automatically fulfills the first criteria of your manual filtering process (retaining only scores 1 and 2, which are flanked by viral/virus-like genes).

```bash
# Distill the annotations and retain only highly confident AMG candidates (scores 1 and 2)
DRAM-v.py distill \
    -i 8.AMG/dramv_annotate/annotations.tsv \
    -o 8.AMG/dramv_distill \
    --max_auxiliary_score 2
```

*(Note: The remaining strict filtering steps—excluding terminal genes with the "F" flag and removing genes associated with core viral processes or misidentified functions—were performed manually in Excel.)*

## 4. Visualizing the Genomic Context of AMGs

To further validate the identified AMGs and assess their genomic neighborhoods, we visualized the gene arrangement of AMG-carrying vOTUs. For each representative vOTU, the coordinates, strand, and functional annotation of all predicted genes were extracted from the DRAM-v output. Each gene was assigned to one of the following categories: AMGs, viral or virus-like genes, other functional genes, hypothetical proteins, and unannotated genes. 

Genome organization diagrams were generated in R (v4.3.1) using the `gggenes` package. Genes were rendered as strand-oriented arrows and colored according to their assigned functional category. To highlight the metabolic potential, only AMGs were annotated with their specific gene names (KO identifiers), whereas other genes were displayed without labels.


Executed the data preparation script ([01_prepare_context.R](./01_prepare_context.R)) to extract the features, align the lists, and assign the plotting categories.
Executed the visualization script ([02_plot_gggenes.R](./02_plot_gggenes.R)) to automatically construct and export the `gggenes` diagrams for each vOTU into the `context_plots` folder.

*(Note: Both R scripts were performed in my Desktop (Windows10))*
