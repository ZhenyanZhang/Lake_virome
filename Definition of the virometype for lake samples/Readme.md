# Definition of the virometype

We characterized the shared patterns of viruses in global lakes by comparing the detection of vOTUs in each pair of lakes. To further explore the biogeographic patterns of lake viromes, we grouped the lake samples into distinct “virometypes” by their viral community composition, following the “enterotype” concept in microbiome research. The pairwise Jensen-Shannon distance (JSD) was calculated based on the abundance of vOTUs to quantify the differences in overall viral community structure between samples. Partitioning around medoids (PAM) clustering was then applied to the resulting JSD matrix, and the optimal number of clusters was determined by the Calinski-Harabasz index. To statistically validate the clustering partitions, between-class analysis was performed, built upon the principal coordinate analysis of the JSD matrix. The significance of the between-class separation was further assessed by a Monte-Carlo permutation test (P = 0.001; 999 permutations). Analyses related to PAM clustering were performed in Python (v3.9) using the scikit-learn and SciPy packages, except between-class analysis, which was performed in R (v4.3.1) using the ade4 package.

## 1. JSD Calculation and PAM Clustering

Execute the Python script ([pam_clustering.py](./pam_clustering.py)) to calculate the Jensen-Shannon distance (JSD) matrix based on the filtered abundance of vOTUs, perform PAM clustering, and determine the optimal number of clusters.

## 2. Between-Class Analysis (BCA) and Visualization

Execute the R script ([bca_analysis.R](./bca_analysis.R)) to perform principal coordinate analysis (PCoA) and between-class analysis (BCA). The script also generates the BCA spider plot and conducts a Monte-Carlo permutation test.

*(Both scripts were performed in my Desktop (Windows 10).)
