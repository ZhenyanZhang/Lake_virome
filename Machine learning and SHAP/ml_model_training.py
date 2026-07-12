import os
import time
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.metrics import accuracy_score, classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import Pipeline
from joblib import dump

# Global configurations
RANDOM_STATE = 123
BASE_DIR = 'D:/onedrive/zzy/isme_revised/1.overview/factor'
DATA_PATH = f'{BASE_DIR}/data_for_ml.csv'
MODEL_SAVE_PATH = f'{BASE_DIR}/models'
RESULTS_PATH = f'{BASE_DIR}/ml_results'

os.makedirs(MODEL_SAVE_PATH, exist_ok=True)
os.makedirs(RESULTS_PATH, exist_ok=True)

# Load and preprocess data
data = pd.read_csv(DATA_PATH)
X = data.iloc[:, 1:14]
y = data.iloc[:, 14]

le = LabelEncoder()
y_encoded = le.fit_transform(y)

print("Label Encoding Map:", dict(zip(le.classes_, le.transform(le.classes_))))

# Data splitting
X_train, X_test, y_train, y_test = train_test_split(
    X, y_encoded, test_size=0.2, random_state=RANDOM_STATE, stratify=y_encoded
)

# Define models
models = {
    'XGBoost': XGBClassifier(random_state=RANDOM_STATE),
    'Random Forest': RandomForestClassifier(random_state=RANDOM_STATE),
    'Logistic Regression': LogisticRegression(solver='saga', random_state=RANDOM_STATE)
}

# Define hyperparameter grids
param_grids = {
    'XGBoost': {
        'xgboost__n_estimators': [100, 200, 300, 400, 500],
        'xgboost__max_depth': [3, 6, 9],
        'xgboost__learning_rate': [0.01, 0.1, 0.2],
        'xgboost__subsample': [0.8, 1.0]
    },
    'Random Forest': {
        'randomforest__n_estimators': [100, 200, 300, 400, 500],
        'randomforest__max_depth': [None, 10, 20],
        'randomforest__min_samples_split': [2, 5, 10]
    },
    'Logistic Regression': {
        'logisticregression__penalty': ['l1', 'l2'],
        'logisticregression__C': [0.001, 0.01, 0.1, 1, 10, 100]
    }
}

cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=RANDOM_STATE)

best_results = []
all_cv_results_list = []

# Training and evaluation
for model_name, model in models.items():
    print(f"Training {model_name}...")
    start = time.time()

    pipeline_step_name = model_name.lower().replace(' ', '')
    pipeline = Pipeline([('scaler', StandardScaler()), (pipeline_step_name, model)])

    grid_search = GridSearchCV(
        pipeline, param_grids[model_name],
        cv=cv, scoring='accuracy', refit=True, n_jobs=-1
    )
    grid_search.fit(X_train, y_train)

    cv_res = pd.DataFrame(grid_search.cv_results_)
    cv_res['Model'] = model_name
    all_cv_results_list.append(cv_res[['Model', 'params', 'mean_test_score', 'std_test_score']])

    y_pred = grid_search.predict(X_test)
    acc = accuracy_score(y_test, y_pred)

    print(f"Best Parameters: {grid_search.best_params_}")

    # Save best model
    model_filename = f"{MODEL_SAVE_PATH}/{model_name.replace(' ', '_')}_best_model.pkl"
    dump(grid_search.best_estimator_, model_filename)

    best_results.append({
        'Model': model_name,
        'Best Parameters': str(grid_search.best_params_),
        'Test Accuracy': acc
    })
    print(f"{model_name} Best Test Accuracy: {acc:.4f} (Time: {time.time() - start:.1f}s)\n")

# Export results
best_results_df = pd.DataFrame(best_results)
best_results_df.to_excel(f'{RESULTS_PATH}/classification_best_results.xlsx', index=False)

all_cv_results_df = pd.concat(all_cv_results_list, ignore_index=True)
params_df = all_cv_results_df['params'].apply(pd.Series)
params_df.columns = [col.split('__')[-1] for col in params_df.columns]

all_cv_results_df = pd.concat([all_cv_results_df.drop('params', axis=1), params_df], axis=1)
all_cv_results_df.rename(columns={'mean_test_score': 'CV_Mean_Accuracy', 'std_test_score': 'CV_Std_Accuracy'},
                         inplace=True)
all_cv_results_df.to_excel(f'{RESULTS_PATH}/all_hyperparameter_combinations_accuracy.xlsx', index=False)