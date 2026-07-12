# Machine learning and SHAP

## 1. Feature Selection

To minimize the impact of multicollinearity between collected variables in subsequent analyses, variance inflation factors (VIF) were calculated. Variables were iteratively removed until only those with a VIF < 10 were retained.

Execute the Python script ([vif_calculation.py](./vif_calculation.py)) to process the data and output the reserved variables.


## 2. Machine Learning Construction and Evaluation

Machine learning was subsequently applied to construct an accurate classification model and evaluate the simultaneous effects of the 13 variables on the biogeographic pattern of the lake virome (specifically, the virometype of each sample). 

Hyperparameter tuning was performed for multiple algorithms (random forest, logistic regression, and XGBoost) using 10-fold cross-validation to select the best combination of algorithms and hyperparameters. 

Execute the model training script ([ml_model_training.py](./ml_model_training.py)) to run the cross-validation grid search and save the best classification models.

After training, execute the evaluation script ([ml_model_evaluation.py](./ml_model_evaluation.py)) to confirm the performance of the final model by generating the receiver operating characteristic (ROC) curves and the confusion matrix.


## 3. Interpreting Model Outputs via SHAP

To quantify the contribution of anthropogenic activities to the lake virome, we applied SHapley Additive exPlanations (SHAP). SHAP values were computed with the `TreeExplainer` algorithm to obtain class-specific values indicating each predictor's contribution to the predicted probability of each virometype.

The overall importance was determined by summing the class-specific mean absolute SHAP values. Beeswarm summary plots were used to visualize the direction and class specificity of each predictor’s effect, and SHAP dependence plots fitted with locally weighted regression (LOESS) were generated for the most important predictors to characterize response thresholds.

Execute the SHAP analysis script ([shap_analysis.py](./shap_analysis.py)) to interpret the Random Forest model and automatically generate all related visualizations.

