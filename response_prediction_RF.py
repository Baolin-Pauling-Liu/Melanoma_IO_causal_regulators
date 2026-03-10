import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.preprocessing import LabelEncoder

df = pd.read_csv('/home/pauling/projects/01.melanoma/09.data.new/22.all.causal.feature.csv')
df = df[df["Timepoint"] == "pre"].copy()    # Or use all samples

X = df.drop(columns=["sample", "Timepoint","response"]) 
y = df['response']   
label_encoder = LabelEncoder()
y = label_encoder.fit_transform(y)

accuracies = []
auc_scores = []
rf_imps = []

for i in range(100):
    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=0.3, stratify=y, random_state=i
    )

    rf = RandomForestClassifier(n_estimators=200, random_state=i)
    rf.fit(X_train, y_train)

    y_pred = rf.predict(X_val)
    accuracies.append(accuracy_score(y_val, y_pred))

    y_pred_proba = rf.predict_proba(X_val)[:, 1]
    auc_scores.append(roc_auc_score(y_val, y_pred_proba))

    rf_imps.append(rf.feature_importances_)  # 每次都存下来

accuracies = np.array(accuracies)
auc_scores = np.array(auc_scores)

rf_imps = np.vstack(rf_imps)  # shape: (100, n_features)
feature_names = X.columns if hasattr(X, "columns") else [f"f{i}" for i in range(X.shape[1])]

imp_df = pd.DataFrame(rf_imps, columns=feature_names)
imp_long = imp_df.reset_index(names="iter").melt(
    id_vars="iter", var_name="feature", value_name="importance"
)

imp_long.to_csv("/home/liubaoli/projects/01.melanoma/02.data/19.classification/imp_long_all_samples.csv", index=False)

auc_df = pd.DataFrame({
    "iter": range(len(auc_scores)),
    "auc": auc_scores
})
auc_df.to_csv("/home/liubaoli/projects/01.melanoma/02.data/19.classification/auc_scores_all_samples.csv", index=False)
