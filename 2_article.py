import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix

df = pd.read_csv("ab_chains_with_indices_with_features .csv")

y_true = df['ptm?'].copy()

df['ptm?'] = 'No'
for index, row in df.iterrows():
    if row['Is in Beta Sheet'] == True:
        if (row['pH'] != 'Both' and float(row['pH']) < 7.0 or row['pH'] == 'Both') and row['Temperature'] <= 100:
            df.at[index, 'ptm?'] = 'No'
        elif (row['pH'] != 'Both' and float(row['pH']) >= 7.0 or row['pH'] == 'Both') and row['Temperature'] > 100:
            df.at[index, 'ptm?'] = 'No'
    elif row['Is in Beta Sheet'] == False:
        if row['C_gamma'] > 3.4:
            if (row['pH'] != 'Both' and float(row['pH']) < 7.0 or row['pH'] == 'Both') and row['Temperature'] <= 100:
                df.at[index, 'ptm?'] = 'No'
            elif (row['pH'] != 'Both' and float(row['pH']) >= 7.0 or row['pH'] == 'Both') and row['Temperature'] > 100:
                if row['Motif'] in ['NG', 'NH', 'NS']:
                    df.at[index, 'ptm?'] = 'Yes'
                else:
                    df.at[index, 'ptm?'] = 'No'
        elif row['C_gamma'] <= 3.4:
            if row['SASA'] <= 10:
                if (row['pH'] != 'Both' and float(row['pH']) < 7.0 or row['pH'] == 'Both') and row['Temperature'] <= 100:
                    df.at[index, 'ptm?'] = 'No'
                elif (row['pH'] != 'Both' and float(row['pH']) >= 7.0 or row['pH'] == 'Both') and row['Temperature'] > 100:
                    df.at[index, 'ptm?'] = 'Yes'
            elif row['SASA'] > 10:
                if (row['pH'] != 'Both' and float(row['pH']) < 7.0 or row['pH'] == 'Both') and row['Temperature'] <= 100:
                    df.at[index, 'ptm?'] = 'Yes'
                elif (row['pH'] != 'Both' and float(row['pH']) >= 7.0 or row['pH'] == 'Both') and row['Temperature'] > 100:
                    df.at[index, 'ptm?'] = 'Yes'

y_pred = df['ptm?']

y_true = y_true.map({True: 1, False: 0})
y_pred = y_pred.map({'Yes': 1, 'No': 0})

print(y_true)
print(y_pred)
accuracy = accuracy_score(y_true, y_pred)
precision = precision_score(y_true, y_pred )
recall = recall_score(y_true, y_pred)
f1 = f1_score(y_true, y_pred)
auc_roc = roc_auc_score(y_true, y_pred)

tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

print(f'Accuracy: {accuracy}')
print(f'Precision: {precision}')
print(f'Recall: {recall}')
print(f'F1-score: {f1}')
print(f'AUC-ROC: {auc_roc}')
print(f'TP: {tp}, TN: {tn}, FP: {fp}, FN: {fn}')


