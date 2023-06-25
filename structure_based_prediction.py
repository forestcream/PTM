import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score,  confusion_matrix


def predict_based_on_struvture(input_path="ab_chains_with_indices_with_features .csv"):
    '''
    Предсказывает ПТМ, основываясь на статье https://doi.org/10.1080/19420862.2018.1478646
    :param input_path: путь к датасету
    :return: None
    '''
    df = pd.read_csv(input_path)
    y_true = df['ptm?'].copy()

    df['ptm?'] = 'No'
    for index, row in df.iterrows():
        if row['Is in Beta Sheet']:
            if (row['pH'] != 'Both' and float(row['pH']) < 7.0 or row['pH'] == 'Both') and row['Temperature'] <= 100:
                df.at[index, 'ptm?'] = 'No'
            elif (row['pH'] != 'Both' and float(row['pH']) >= 7.0 or row['pH'] == 'Both') and row['Temperature'] > 100:
                df.at[index, 'ptm?'] = 'No'
        elif not row['Is in Beta Sheet']:
            if row['C_gamma'] > 3.4:
                if (row['pH'] != 'Both' and float(row['pH']) < 7.0 or row['pH'] == 'Both') and row[
                    'Temperature'] <= 100:
                    df.at[index, 'ptm?'] = 'No'
                elif (row['pH'] != 'Both' and float(row['pH']) >= 7.0 or row['pH'] == 'Both') and row[
                    'Temperature'] > 100:
                    if row['Motif'] in ['NG', 'NH', 'NS']:
                        df.at[index, 'ptm?'] = 'Yes'
                    else:
                        df.at[index, 'ptm?'] = 'No'
            elif row['C_gamma'] <= 3.4:
                if row['SASA'] <= 10:
                    if (row['pH'] != 'Both' and float(row['pH']) < 7.0 or row['pH'] == 'Both') and row[
                        'Temperature'] <= 100:
                        df.at[index, 'ptm?'] = 'No'
                    elif (row['pH'] != 'Both' and float(row['pH']) >= 7.0 or row['pH'] == 'Both') and row[
                        'Temperature'] > 100:
                        df.at[index, 'ptm?'] = 'Yes'
                elif row['SASA'] > 10:
                    if (row['pH'] != 'Both' and float(row['pH']) < 7.0 or row['pH'] == 'Both') and row[
                        'Temperature'] <= 100:
                        df.at[index, 'ptm?'] = 'Yes'
                    elif (row['pH'] != 'Both' and float(row['pH']) >= 7.0 or row['pH'] == 'Both') and row[
                        'Temperature'] > 100:
                        df.at[index, 'ptm?'] = 'Yes'

    y_pred = df['ptm?']

    y_true = y_true.map({True: 1, False: 0})
    y_pred = y_pred.map({'Yes': 1, 'No': 0})

    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    print(f'Accuracy: {accuracy}')
    print(f'Precision: {precision}')
    print(f'Recall: {recall}')
    print(f'TP: {tp}, TN: {tn}, FP: {fp}, FN: {fn}')


predict_based_on_struvture()
