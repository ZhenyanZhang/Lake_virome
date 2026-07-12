import os
import joblib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, label_binarize
from sklearn.metrics import confusion_matrix, roc_curve, auc

# Data preparation
data = pd.read_csv('D:/onedrive/zzy/isme_revised/1.overview/factor/data_for_ml.csv')
X = data.iloc[:, 1:14]
y = data.iloc[:, 14]

le = LabelEncoder()
y_encoded = le.fit_transform(y)
classes = le.classes_

X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.2, random_state=42, stratify=y_encoded)

y_test_bin = label_binarize(y_test, classes=np.unique(y_encoded))
n_classes = y_test_bin.shape[1]

# Load model and evaluate
model_path = 'D:/onedrive/zzy/isme_revised/1.overview/factor/models'
model_files = [f for f in os.listdir(model_path) if f.endswith('Random_Forest_best_model.pkl')]

for file in model_files:
    print(f"\n--- Evaluating: {file} ---")
    model = joblib.load(os.path.join(model_path, file))

    y_pred = model.predict(X_test)
    y_score = model.predict_proba(X_test)

    # Confusion Matrix
    cm = confusion_matrix(y_test, y_pred, normalize='true')
    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt='.3f', cmap='Blues', xticklabels=classes, yticklabels=classes)
    plt.title(f'Confusion Matrix: {file}')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.show()

    # ROC Curves (One-vs-Rest)
    plt.figure(figsize=(8, 6))
    for i in range(n_classes):
        fpr, tpr, _ = roc_curve(y_test_bin[:, i], y_score[:, i])
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f'Class {classes[i]} (AUC = {roc_auc:.2f})')

    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curves: {file}')
    plt.legend(loc="lower right")
    plt.show()