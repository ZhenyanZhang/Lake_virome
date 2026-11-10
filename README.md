# Global biogeography and ecological significance of lake virome

This repository contains the bioinformatic workflows, custom scripts, parameter configurations, and step-by-step instructions to reproduce the principal analyses in our manuscript: **"Global biogeography and ecological significance of lake virome"**.

## 📌 Overview
Our analytical framework processes 706 global lake metagenomes to construct the lake virome database (42,141 species-level vOTUs), assign viral taxonomy, define community-level "virometypes" along climatic and anthropogenic gradients via machine learning and SHAP, infer host–virus interactions and virus-mediated horizontal gene transfer (HGT), and characterize putative auxiliary metabolic genes (AMGs).

## 🛠️ System Requirements & Environments
* **Operating Systems:** Ubuntu 22.04.3 LTS (Linux HPC environment for Bash, Perl, and Python workflows) and Windows 10 (for R-based statistical analyses and visualization).
* **Programming Languages:** Python (v3.9), R (v4.3.1), Perl (v5.34), and Bash.

### Software Versions
* **Quality Control & Assembly:** fastp (v0.23.4), MEGAHIT (v1.2.9)
* **Viral Identification, Dereplication & Abundance:** VIBRANT (v1.2.1), CheckV (v1.0.1), CD-HIT (v4.7), Bowtie2 (v2.5.4), SAMtools (v1.3.1), CoverM (v0.5.0)
* **Gene Calling, Clustering & Taxonomy:** Prodigal (v2.6.3), DIAMOND (v2.1.9.163), MCL (v22-282), HMMER (v3.4), geNomad (v1.9.1), BLASTn (v2.1.5 / v2.17.0+)
* **Host Recovery, Prediction & HGT:** MetaBAT2 (v2.12.1), dRep (v3.4.5), GTDB-Tk (v2.4.0), iPHoP (v1.4.1), MetaCHIP (v1.10.13), Ranger-DTL (v2.0), eggNOG-mapper (v2.1.13), Kraken2 (v2.1.2)
* **AMG Annotation, Phylogeny & Structure:** VirSorter2 (v2.2.4), DRAM-v (v1.5.0), MAFFT (v7.471), trimAl (v1.4.rev15), IQ-TREE (v2.0.3), AlphaFold3 (AlphaFold Server), PyMOL (v3.2.0a), TM-align
* **Statistical & Visualization Tools:** Python packages (`scikit-learn`, `scikit-learn-extra`, `SciPy`, `statsmodels`, `xgboost`, `shap`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `Basemap`), R packages (`rtk`, `ade4`, `phyloseq`, `geosphere`, `adespatial`, `vegan`, `psych`, `ggClusterNet`, `igraph`, `gggenes`, `ggplot2`), Cytoscape (v3.8.2), iTOL (v6.8.1), TBtools (v2.0.42), GraphPad Prism (v8.0.2), IBM SPSS Statistics (v20.0.0)

### Reference Database Versions
* **NCBI RefSeq Viral Protein Database:** v208
* **Viral Orthologous Groups (VOG) Database:** v231 (587 marker HMM profiles)
* **geNomad Reference Database:** v1.9
* **External Viral Databases for Comparison:** IMG/VR v4, GOV 2.0, GPD, GVD
* **GTDB-Tk Reference Database:** release226
* **iPHoP Database:** Jun_2025_pub (supplemented with 3,973 lake MAGs recovered in this study)
* **eggNOG Database:** v5.0.2
* **Kraken2 Standard Database:** Kraken2-PlusPF (release 20260226)

---

## 📂 Repository Structure & Analytical Workflow

Each directory below corresponds to a core analytical module in the manuscript and contains its own detailed `README`, scripts, and parameter settings[cite: 2]:

1. **[AAA_Before all](./AAA_Before%20all)**[cite: 2]  
   Prerequisites, environment preparation, and general notes before running the workflow.
2. **[Construction of LVD](./Construction%20of%20LVD)**[cite: 2]  
   Raw read quality control (`fastp`), metagenomic assembly (`MEGAHIT`), viral contig identification (`VIBRANT` and `CheckV`), species-level vOTU dereplication (`CD-HIT`), 50-million-read subsampling (`Bowtie2` and `SAMtools`), abundance quantification (`CoverM`), and rarefaction/database comparisons.
3. **[Taxonomic assignment of vOTUs](./Taxonomic%20assignment%20of%20vOTUs)**[cite: 2]  
   Gene calling (`Prodigal`), protein-sharing network clustering (`DIAMOND` and `MCL`), and the integrated hierarchical taxonomic workflow combining NCBI RefSeq homology, VOG HMM profiling, `geNomad`, and genus-level LCA assignment.
4. **[Definition of the virometype for lake samples](./Definition%20of%20the%20virometype%20for%20lake%20samples)**[cite: 2]  
   Jensen-Shannon distance (JSD) calculation, PAM clustering, Calinski-Harabasz index evaluation, between-class analysis, and the six-dimensional sensitivity analyses (alternative metrics, algorithms, filtering thresholds, normalizations, geographic/study-level subsampling, and repeated subsampling).
5. **[Machine learning and SHAP](./Machine%20learning%20and%20SHAP)**[cite: 2]  
   Multicollinearity filtering (VIF), hyperparameter tuning with 10-fold cross-validation (`Random Forest`, `XGBoost`, `Logistic Regression`), `SHAP` value computation (`TreeExplainer`), LOESS threshold fitting, leave-one-lake-out spatial cross-validation, and variation partitioning analysis (VPA).
6. **[Host prediction and HGT detection](./Host%20prediction%20and%20HGT%20detection)**[cite: 2]  
   Host MAG recovery and dereplication (`MetaBAT2`, `dRep`, `GTDB-Tk`), custom `iPHoP` host prediction, cross-genus HGT inference (`MetaCHIP` and `Ranger-DTL`), and identification of virus-carried host genes and virus-mediated HGT events.
7. **[AMG Identification](./AMG%20Identification)**[cite: 2]  
   Formatting and putative AMG annotation (`VirSorter2` and `DRAM-v`), manual curation rules, genomic context visualization (`gggenes`), phylogenetic tree reconstruction (`MAFFT`, `trimAl`, `IQ-TREE`), and protein structural alignment (`PyMOL` and `TM-align`).
8. **[Visualization](./Visualization)**[cite: 2]  
   Scripts for microbial community profiling, cross-domain and bacterial co-occurrence network construction (`ggClusterNet` and `igraph`), network robustness tests, and statistical plotting.

---

## 🚀 Usage Instructions
Clone this repository to your local machine or HPC cluster:
```bash
git clone [https://github.com/ZhenyanZhang/Lake_virome.git](https://github.com/ZhenyanZhang/Lake_virome.git)
cd Lake_virome
