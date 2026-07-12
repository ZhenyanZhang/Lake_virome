import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor

# Load data
data_df = pd.read_csv('D:/onedrive/zzy/isme_revised/1.overview/factor/factor_for_vif.csv')

def calculate_vif(data_frame):
    vif_data = pd.DataFrame()
    vif_data["Variable"] = data_frame.columns
    vif_data["VIF"] = [variance_inflation_factor(data_frame.values, i) for i in range(data_frame.shape[1])]
    return vif_data

# Iteratively remove variables with VIF >= 10
while True:
    vif_result = calculate_vif(data_df)
    max_vif_variable = vif_result.loc[vif_result['VIF'].idxmax(), 'Variable']
    max_vif = vif_result['VIF'].max()

    if max_vif < 10:
        break

    print(f"Removing variable with highest VIF: {max_vif_variable} (VIF={max_vif:.2f})")
    data_df = data_df.drop(max_vif_variable, axis=1)

print("Final variables (VIF < 10):")
print(data_df.columns.tolist())

# Save results
data_df.to_excel('D:/onedrive/zzy/isme_revised/1.overview/factor/final_data_with_vif.xlsx', index=False)