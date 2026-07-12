import os
import joblib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
import shap
import statsmodels.api as sm
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['font.family'] = 'Arial'

# Path configurations and data loading
data_path = 'D:/onedrive/zzy/isme_revised/1.overview/factor/data_for_ml.csv'
model_path = 'D:/onedrive/zzy/isme_revised/1.overview/factor/models/Random_Forest_best_model.pkl'

base_out_dir = 'D:/onedrive/zzy/isme_revised/1.overview/factor/results_shap'
dirs = {
    'global': os.path.join(base_out_dir, '01_Global_Importance'),
    'beeswarm': os.path.join(base_out_dir, '02_Beeswarm_Plots'),
    'dependence': os.path.join(base_out_dir, '03_Dependence_Plots')
}
for d in dirs.values():
    os.makedirs(d, exist_ok=True)

data = pd.read_csv(data_path)
X = data.iloc[:, 1:14]
y = data.iloc[:, 14]
feature_names = X.columns.tolist()

le = LabelEncoder()
y_encoded = le.fit_transform(y)
target_names = le.classes_

# Load model and compute SHAP
print("Loading model and calculating SHAP values...")
best_pipeline = joblib.load(model_path)
scaler = best_pipeline.named_steps['scaler']
rf_model = best_pipeline.named_steps['randomforest']

X_scaled_df = pd.DataFrame(scaler.transform(X), columns=feature_names)
explainer = shap.TreeExplainer(rf_model)
shap_values = explainer.shap_values(X_scaled_df)

def get_class_shap_values(shap_vals, class_idx):
    if isinstance(shap_vals, list):
        return shap_vals[class_idx]
    elif isinstance(shap_vals, np.ndarray) and len(shap_vals.shape) == 3:
        return shap_vals[:, :, class_idx]
    else:
        return shap_vals[:, :, class_idx].values

# Global Ecological Drivers
print("\nGenerating Global Importance Plot...")
total_imp = pd.Series(0.0, index=feature_names)
for i in range(len(target_names)):
    total_imp += np.abs(get_class_shap_values(shap_values, i)).mean(axis=0)

imp = total_imp.sort_values(ascending=False)
feats, vals = imp.index.tolist(), imp.values
pcts = 100 * vals / vals.sum()

colors = [plt.cm.RdBu_r(i / max(len(feats) - 1, 1)) for i in range(len(feats))]
fig, ax = plt.subplots(figsize=(10, 8))
ax.barh(np.arange(len(feats)), vals, color=colors, edgecolor='none')
ax.set_yticks(np.arange(len(feats)))
ax.set_yticklabels(feats, fontsize=11)
ax.invert_yaxis()
ax.set_xlabel('Mean |SHAP| value', fontsize=14)
for s in ['top', 'right']: ax.spines[s].set_visible(False)

# Donut chart inset
ax_pie = ax.inset_axes([0.42, 0.04, 0.56, 0.56])
wedges, _ = ax_pie.pie(vals, colors=colors, startangle=90, counterclock=False, wedgeprops=dict(width=0.42, edgecolor='white', linewidth=1.0))
ax_pie.text(0, 0, 'SHAP %', ha='center', va='center', fontsize=15)
for w, p in zip(wedges, pcts):
    ang = np.deg2rad((w.theta1 + w.theta2) / 2)
    ax_pie.text(1.18 * np.cos(ang), 1.18 * np.sin(ang), f'{p:.1f}%', ha='center', va='center', fontsize=7.5)
ax_pie.set(aspect='equal')

plt.tight_layout()
plt.savefig(os.path.join(dirs['global'], 'Global_Importance_Bar.pdf'), format='pdf', dpi=300, bbox_inches='tight')
plt.close()

# Beeswarm Plots
print("\nGenerating Beeswarm Plots...")
for i, vt_name in enumerate(target_names):
    plt.figure(figsize=(10, 8))
    shap.summary_plot(get_class_shap_values(shap_values, i), X, show=False)
    plt.title(f'Ecological Drivers for {vt_name} Assembly', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(dirs['beeswarm'], f'Beeswarm_{vt_name}.pdf'), format='pdf', dpi=300)
    plt.close()

# Dependence Plots (Top 5)
print("\nGenerating Dependence Plots...")
from matplotlib.transforms import blended_transform_factory

for i, vt_name in enumerate(target_names):
    shap_vals_class = get_class_shap_values(shap_values, i)
    top_5_idx = np.argsort(np.abs(shap_vals_class).mean(axis=0))[-5:][::-1]

    for rank, feat_idx in enumerate(top_5_idx):
        feat_name = feature_names[feat_idx]
        x_vals, y_vals = X[feat_name].values, shap_vals_class[:, feat_idx]

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(x_vals, y_vals, alpha=0.35, s=12, color='#4a90d9', edgecolors='none', zorder=1)
        ax.axhline(0, color='black', linestyle='--', linewidth=1.0, zorder=2)

        try:
            z = sm.nonparametric.lowess(y_vals, x_vals, frac=0.3, it=2)
            z_x, z_y = z[:, 0], z[:, 1]
            ax.plot(z_x, z_y, color='#d62728', linewidth=2.2, zorder=3)

            # Detect zero-crossings
            sign_arr = np.sign(z_y)
            sign_arr[sign_arr == 0] = 1
            for ci in np.where(np.diff(sign_arr))[0]:
                if abs(z_y[ci + 1] - z_y[ci]) >= 1e-12:
                    x_cross = z_x[ci] - z_y[ci] * (z_x[ci + 1] - z_x[ci]) / (z_y[ci + 1] - z_y[ci])
                    ax.axvline(x_cross, color='#888888', linestyle=':', linewidth=1.2, alpha=0.9, zorder=2)
                    ax.text(x_cross, 0.97, f'{x_cross:.2f}', transform=blended_transform_factory(ax.transData, ax.transAxes),
                            fontsize=8.5, ha='center', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
        except Exception as e:
            print(f"  [WARN] LOESS failed for {feat_name}: {e}")

        ax.set_title(f'{vt_name} — {feat_name}', fontsize=13, fontweight='bold')
        ax.set_xlabel(feat_name, fontsize=11)
        ax.set_ylabel('SHAP value', fontsize=11)
        sns.despine(ax=ax)
        plt.tight_layout()
        plt.savefig(os.path.join(dirs['dependence'], f'Dependence_{vt_name}_Rank{rank + 1}_{feat_name}.pdf'), format='pdf', dpi=300, bbox_inches='tight')
        plt.close()