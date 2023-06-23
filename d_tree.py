import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix
from sklearn.tree import plot_tree
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
from sklearn.impute import SimpleImputer
import joblib

data = pd.read_csv('dataset.csv')

ptm_yes_count = (data['ptm'] == 'yes').sum()
print(f"The number of 'yes' values in the ptm column is: {ptm_yes_count}")
#data = data.dropna()

data['ptm'] = data['ptm'].map({'yes': 1, 'no': 0})
data = pd.get_dummies(data, columns=['pH'])

data = pd.get_dummies(data, columns=['Motif'])

data['Is in Beta Sheet'] = data['Is in Beta Sheet'].map({'Yes': 1, 'No': 0})

feature_cols = ['SASA', 'C_gamma', 'Is in Beta Sheet'] + [col for col in data.columns if col.startswith('Motif') or col.startswith('pH')] + ['Phi', 'Psi', 'Chi1', 'Chi2', 'Hydrogen bonds to side chain', 'Temperature']
X = data[feature_cols]
y = data['ptm']
print(y.sum())
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5)
print(y_test.sum())
#
# clf = DecisionTreeClassifier(class_weight='balanced', ccp_alpha=0.01)
# clf.fit(X_train, y_train)
#
# y_pred = clf.predict(X_test)
imputer = SimpleImputer(strategy='mean')

# Fit the imputer to the training data
imputer.fit(X_train)

# Transform the training and test data
X_train_imputed = imputer.transform(X_train)
X_test_imputed = imputer.transform(X_test)
rf_clf = RandomForestClassifier(class_weight='balanced', ccp_alpha=0.01)
rf_clf.fit(X_train_imputed, y_train)
joblib.dump(rf_clf, 'rf_clf.joblib')
y_pred = rf_clf.predict(X_test_imputed)


accuracy = accuracy_score(y_test, y_pred)
print(f'Accuracy: {accuracy}')
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)
auc_roc = roc_auc_score(y_test, y_pred)

tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()

print(f'Precision: {precision}')
print(f'Recall: {recall}')
print(f'F1-score: {f1}')
print(f'AUC-ROC: {auc_roc}')
print(f'TP: {tp}, TN: {tn}, FP: {fp}, FN: {fn}')

# plt.figure(figsize=(20,10))
# plot_tree(clf, feature_names=feature_cols, class_names=['0', '1'], filled=True)
# plt.show()
plt.figure(figsize=(20,10))
plot_tree(rf_clf.estimators_[0], feature_names=feature_cols, class_names=['0', '1'], filled=True)
plt.show()