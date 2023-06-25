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


def print_how_many_ptm_in_dataset(data, ptm_column = 'ptm?'):
    '''
    Считает количесвто ПТМ в таблице
    :param data: Таблица
    :param ptm_column: Название колонки с данными о ПТМ
    :return: None
    '''
    ptm_yes_count = (data[ptm_column] == True).sum()
    print(f"The number of 'yes' values in the ptm column is: {ptm_yes_count}")


def print_metrics(y_test, y_pred):
    '''
    Выводит различные метрки
    :param y_test:
    :param y_pred:
    :return:
    '''
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


def visualise_tree(model, feature_cols):
    '''
    строит визуализацию дерева решений (случайного леса)
    :param model: модель
    :param feature_cols: название фичей
    :return: Nont
    '''
    plt.figure(figsize=(20, 10))
    plot_tree(model.estimators_[0], feature_names=feature_cols, class_names=['0', '1'], filled=True)
    plt.show()

def run_on_article_data(input_file='mmc2(3) (2).xlsx'):
    '''
    Строит случайный лес по таблице из статьи
    :param input_file:
    :return: None
    '''
    data = pd.read_excel(input_file)
    print_how_many_ptm_in_dataset(data, ptm_column='is_deam')

    data['ptm'] = data['is_deam'].map({True: 1, False: 0})

    data['Is in Beta Sheet'] = data['SHEET']
    data = pd.get_dummies(data, columns=['N+1'])  # to numerical from categorical
    feature_cols = ['SASA', 'C_NAttackDistance', 'Is in Beta Sheet', 'phi', 'psi', 'chi1', 'chi2', 'NumH_bonds']

    for index, row in data.iterrows():
        if index < 10:
            print(row)
    X = data[feature_cols]
    y = data['ptm']
    print(y.sum())
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5)
    joblib.dump((X_train, X_test, y_train, y_test), 'test_data_from_second_article.joblib')

    imputer = SimpleImputer(strategy='mean')

    imputer.fit(X_train)

    X_train_imputed = imputer.transform(X_train)
    X_test_imputed = imputer.transform(X_test)

    rf_clf = RandomForestClassifier(class_weight='balanced', ccp_alpha=0.01)
    rf_clf.fit(X_train_imputed, y_train)
    joblib.dump(rf_clf, 'rf_clf_on_article_data_with_cpp.joblib')
    y_pred = rf_clf.predict(X_test_imputed)

    print_metrics(y_test, y_pred)
    visualise_tree(rf_clf, feature_cols)

def build_random_forest(input_file='ab_chains_with_indices_with_features .csv', size_of_test=0.5,
                        use_loaded_model=False, path_to_loaded_model='rf_clf_on_article_data_2_balanced.joblib',
                        use_loaded_test=False):
    '''
    Строит случайный лес на посчитанных данных
    :param input_file: путь до датасета
    :param size_of_test: велечина тестовой выборки
    :param use_loaded_model: использовать или нет скаченную модель (по умолчанию False)
    :param path_to_loaded_model: путь до скаченной модели
    :param use_loaded_test: использовать ли уже скаченные тыстовые данные
    :return:
    '''
    data = pd.read_csv(input_file)
    print_how_many_ptm_in_dataset(data)

    data['ptm?'] = data['ptm?'].map({True: 1, False: 0})
    data = pd.get_dummies(data, columns=['pH'])
    data = pd.get_dummies(data, columns=['Motif'])

    data['Is in Beta Sheet'] = data['Is in Beta Sheet'].map({True: 1, False: 0})

    feature_cols = ['SASA', 'C_gamma', 'Is in Beta Sheet'] + ['Phi', 'Psi', 'Chi1', 'Chi2',
                                                              'Hydrogen bonds to side chain']

    # feature_cols = ['SASA', 'C_gamma', 'Is in Beta Sheet'] + [col for col in data.columns if col.startswith('Motif') or col.startswith('pH')] + ['Phi', 'Psi', 'Chi1', 'Chi2', 'Hydrogen bonds to side chain', 'Temperature']
    X = data[feature_cols]
    y = data['ptm?']
    if use_loaded_test == True:
        X_train, X_test, y_train, y_test = joblib.load('test_data_from_my_data.joblib')
    else:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=size_of_test)
        joblib.dump((X_train, X_test, y_train, y_test), 'test_data_from_my_data_50.joblib')

    imputer = SimpleImputer(strategy='mean')

    imputer.fit(X_train)

    X_train_imputed = imputer.transform(X_train)
    X_test_imputed = imputer.transform(X_test)

    if use_loaded_model == True:
        loaded_model = joblib.load(path_to_loaded_model)
        y_pred = loaded_model.predict(X_test_imputed)
    else:
        rf_clf = RandomForestClassifier(class_weight='balanced', ccp_alpha=0.01)
        rf_clf.fit(X_train_imputed, y_train)
        joblib.dump(rf_clf, 'rf_clf_train_my_data_test_my_data_with_cpp.joblib')
        y_pred = rf_clf.predict(X_test_imputed)

    print_metrics(y_test, y_pred)
    visualise_tree(rf_clf, feature_cols)


build_random_forest(input_file='ab_chains_with_indices_with_features .csv', size_of_test=0.5,
                    use_loaded_model=False, path_to_loaded_model='rf_clf_on_article_data_2_balanced.joblib')
