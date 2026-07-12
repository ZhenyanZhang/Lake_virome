import pandas as pd
from scipy.spatial.distance import pdist, squareform, jensenshannon
from sklearn_extra.cluster import KMedoids
from sklearn.metrics import calinski_harabasz_score, silhouette_score
import warnings

# Ignore clustering warnings
warnings.filterwarnings("ignore")


def main():
    # Load abundance data
    df = pd.read_csv("D:/onedrive/zzy/isme_revised/1.overview/filtered_abundance_matrix.csv", index_col=0)

    # Filter low-abundance vOTUs
    occurrence_rate = (df > 0).mean(axis=1)
    rel_df = df.div(df.sum(axis=0), axis=1)
    max_rel_abund = rel_df.max(axis=1)
    df_filtered = df[(occurrence_rate >= 0.05) | (max_rel_abund >= 0.001)]
    df_filtered.to_csv("D:/onedrive/zzy/isme_revised/1.overview/PAM/vOTU_filtered_abundance.csv")

    # Transpose matrix for ML input
    X_raw = df_filtered.T.values.astype(float)
    sample_names = df_filtered.columns.tolist()

    # Calculate JSD distance matrix
    row_sums = X_raw.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    X_rel = X_raw / row_sums

    jsd_dist_array = pdist(X_rel, metric=lambda u, v: jensenshannon(u, v))
    jsd_dist_matrix = squareform(jsd_dist_array)

    # Perform PAM clustering and evaluation
    k_range = range(2, 11)
    evaluation_results = []
    cluster_assignments_dict = {}

    for k in k_range:
        pam = KMedoids(n_clusters=k, metric='precomputed', method='pam', init='heuristic', random_state=42)
        labels = pam.fit_predict(jsd_dist_matrix)

        cluster_assignments_dict[f"JSD_k{k}"] = labels
        ch_index = calinski_harabasz_score(X_rel, labels)
        sil_score = silhouette_score(jsd_dist_matrix, labels, metric="precomputed")

        evaluation_results.append({
            'K': k,
            'Calinski_Harabasz_Index': ch_index,
            'Silhouette_Score': sil_score
        })

    # Save results
    eval_df = pd.DataFrame(evaluation_results)
    eval_df.to_csv("D:/onedrive/zzy/isme_revised/1.overview/PAM/PAM_JSD_evaluation_metrics.csv", index=False)

    assignments_df = pd.DataFrame(cluster_assignments_dict, index=sample_names)
    assignments_df.index.name = "Sample_ID"

    best_row = eval_df.loc[eval_df['Calinski_Harabasz_Index'].idxmax()]
    best_k = int(best_row['K'])

    assignments_df['Best_Cluster_Label'] = assignments_df[f"JSD_k{best_k}"]
    assignments_df.to_csv("D:/onedrive/zzy/isme_revised/1.overview/PAM/PAM_JSD_cluster_assignments.csv")

    pd.DataFrame(jsd_dist_matrix, index=sample_names, columns=sample_names).to_csv(
        "D:/onedrive/zzy/isme_revised/redo/overview/PAM/Distance_Matrix_JSD.csv")


if __name__ == "__main__":
    main()