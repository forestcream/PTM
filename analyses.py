

def get_pH_and_temperature(folder_path="ab_pdb"):
    """
    Добавляем колонку в файл с различными метриками pdb при какой температуре и pH они были сделаны
    :param folder_path: папка с pdb
    :return: None
    """
    import os
    import pandas as pd
    pdb_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]

    df = pd.read_csv('result.csv')
    df['pH'] = None
    df['Temperature'] = None
    for pdb_file in pdb_files:
        print(pdb_file)
        ph = None
        temp = None
        with open(f'ab_pdb\\{pdb_file}', 'r') as f:
            for line in f:
                if line.startswith('REMARK'):
                    if ' PH ' in line.upper():
                        if ph is None:
                            ph_line = line.split()[-1]
                            if ph_line == 'NULL':
                                ph = 'Both'
                            else:
                                try:
                                    ph = float(ph_line)

                                except ValueError:
                                    ph = None

                        else:
                            break
                    elif 'TEMPERATURE' in line.upper():
                        if temp is None:
                            temp_line = line.split()[-1]
                            if temp_line == 'NULL':
                                temp = 'None'
                            else:
                                temp = float(temp_line)
                        else:
                            continue
        print(ph, temp)
        base_name = os.path.basename(pdb_file)
        pdb_id = os.path.splitext(base_name)[0]
        pdb_id = pdb_id.upper()
        mask = df['PDB ID'] == pdb_id
        df.loc[mask, 'pH'] = ph
        df.loc[mask, 'Temperature'] = temp

    df.to_csv('result.csv', index=False)


def do_file_without_duplecates(path_input, path_output):
    """
    Создает файл без повторов одинаковых сиквенсов
    :param path_input: путь файла с сиквенсами
    :param path_output: путь файла для записи
    :return: None
    """
    import csv

    with open(path_input, 'r') as input_file, open(path_output, 'w') as output_file:
        reader = csv.DictReader(input_file)
        writer = csv.DictWriter(output_file, fieldnames=reader.fieldnames, lineterminator='\n')
        writer.writeheader()

        seen = set()
        for row in reader:
            key = (row['PDB ID'], row['Chain'], row['Duplicate Number'])
            if key not in seen:
                writer.writerow(row)
                seen.add(key)




def write_asn_residue_numbers(folder_path="ab_pdb", output_file="all_asn.csv"):
    """
    Создает файл со всеми аспарагинами со всех pdb
    :param folder_path: путь к папке с pdb
    :param output_file: путь к файлу для записи
    :return: None
    """
    from Bio import PDB
    import os
    parser = PDB.PDBParser()
    with open(output_file, 'w') as f:
        for filename in os.listdir(folder_path):
            if filename.endswith('.pdb'):
                pdb_id = filename.split('.')[0]
                structure = parser.get_structure(pdb_id, os.path.join(folder_path, filename))
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.get_resname() == 'ASN':
                                res_num = residue.get_id()[1]
                                f.write(f'{pdb_id}\t{res_num}\n')




def make_dataset():
    """
    Создаает датасет
    :return: None
    """
    import pandas as pd

    result = pd.read_csv('updated_file.csv')
    tmp = pd.read_csv('tmp.csv')

    def pH_equal(pH_result, pH_tmp):
        if pH_tmp == "High pH" and (pH_result == "Both" or float(pH_result) >= 7):
            return True
        elif pH_tmp == "Low pH" and (pH_result == "Both" or float(pH_result) < 7):
            return True
        elif pH_tmp == "Both" and pH_result == "Both":
            return True
        else:
            return False

    dataset = pd.DataFrame()

    for i, row in result.iterrows():
        ptm = 'no'
        for j, row2 in tmp.iterrows():
            # print("i am here")
            # if row['Drug Name'] == row2['name'] and abs(row['updated indexes'] - row2['index']) <= 5 and row[
            #     'Motif'] == row2[
            #     ' motif'] and pH_equal(row['pH'], row2['pH']):
            # if row['Drug Name'] == row2['name'] and pH_equal(row['pH'], row2['pH']) and (row['updated indexes'] == row2['index']):
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if row['Drug Name'] == row2['name'] and (
                    int(row['updated indexes']) == row2['index']):
                ptm = 'yes'
                break
        new_row = row.to_frame().T.assign(ptm=ptm)
        dataset = pd.concat([dataset, new_row], ignore_index=True)



def count_the_same_names():
    """
    Считает количество уникальных названий лекарств
    :return: None
    """
    import pandas as pd
    df1 = pd.read_csv('tmp.csv')
    df2 = pd.read_csv('result.csv')
    print('Column names in df1:', df1.columns)
    print('Column names in df2:', df2.columns)
    merged_df = pd.merge(df1, df2, left_on='name', right_on='Drug Name')
    unique_count = merged_df[['name', 'Drug Name']].nunique()
    unique_count = df1['name'].nunique()
    print(unique_count)

import pandas as pd
import difflib

def longest_common_subsequence(s1, s2):
    matcher = difflib.SequenceMatcher(None, s1, s2)

    matching_blocks = matcher.get_matching_blocks()

    lcs = ''.join(s1[block.a:(block.a + block.size)] for block in matching_blocks)

    return lcs

def drop_from_dataset():
    """
    Выкидывает из датасета сиквенсы не относящиеся к вариабельному дамену или выходящие за его пределы
    :return: None
    """
    import pandas as pd

    dataset = pd.read_csv('dataset.csv')
    updated_indices = pd.read_csv('updated_indices_190623 (2).csv')

    for index, row in dataset.iterrows():

        sequence = row['Sequence']
        ab_seq = max(updated_indices['ab_seq'], key=lambda x: len(longest_common_subsequence(x, sequence)))

        start_index = sequence.find(ab_seq)
        end_index = start_index + len(ab_seq) - 1

        res_num = row['Res_Num']

        if res_num < start_index or res_num > end_index:
            # Drop the row from the dataset
            dataset.drop(index, inplace=True)
    dataset = dataset[dataset['Sequence'].isin(updated_indices['seq'])]
    dataset.to_csv('updated_dataset.csv', index=False)


#do_file_without_duplecates('dataset.csv', 'pdb_without_duplicates.csv')
# split_datatset('result.csv', 'train_with_dup.csv', 'test_with_dup.csv')
# split_datatset('pdb_without_duplicates.csv', 'train_without_dup.csv', 'test_without_dup.csv')
# second_article('train_with_dup.csv', 'train_2_ar_with_dup.csv')
# second_article('train_without_dup.csv', 'train_2_ar_without_dup.csv')
# find_tp_tn_fp_fn('tmp.csv', 'no_ptm_train_2_arr.csv')
# print("Without duplicates-train")
# find_tp_tn_fp_fn('tmp.csv', 'train_2_ar_without_dup.csv')
# print("With duplicates-train")
# find_tp_tn_fp_fn('tmp.csv', 'train_2_ar_with_dup.csv')
# write_asn_residue_numbers()
make_dataset()
#count_the_same_names()
# get_pH_and_temperature()
# compare_folder_and_file()
# file = pd.read_csv('pdb_without_duplicates.csv')
drop_from_dataset()