import pandas as pd


# TODO

# вторая статья
# 1) и/или - посмотреть внимательнее
# 2.1) отдельный датасет только из тех антител, для которых есть пдб  и есть птм хоть один(посмотреть статью, откуда мы брали данные о птм, понять, правда ли, что те анитела, которых нет в таблице, не выдали никаких птмов)
# 2.2) отдельный датасет только из тех антител, для которых есть пдб
# 3.1) попробовать учесть pH антител для датасета из 2.1 (low ph - acidic, high - basic)
# 3.2) попробовать учесть pH антител для датасета из 2.2
#
# третья статья
# 1) сделать decision tree/random forest без pphl
# 1.1) а отдельном только из тех антител, для которых есть пдб  и есть птм хоть один (взять из предыдущего сообщения)
# 1.2) отдельный датасет только из тех антител, для которых есть пдб (взять из предыдущего сообщения)
# 2) сделать decision tree/random forest c pphl
# 2.1) а отдельном только из тех антител, для которых есть пдб  и есть птм хоть один (взять из предыдущего сообщения)
# 2.2) отдельный датасет только из тех антител, для которых есть пдб (взять из предыдущего сообщения)

def get_pH_and_temperature(folder_path="ab_pdb"):
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
                    # print(line)
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


def second_article(path_input, path_output):
    import pandas as pd
    df = pd.read_csv(path_input)
    df['ptm'] = 'No'
    for index, row in df.iterrows():
        if row['Is in Beta Sheet'] == True:
            if (row['pH'] < 7.0 or row['pH'] == 'Both') and row['Temperature'] < 100:
                df.at[index, 'ptm'] = 'No'
            elif (row['pH'] >= 7.0 or row['pH'] == 'Both') and row['Temperature'] >= 100:
                df.at[index, 'ptm'] = 'No'
        elif row['Is in Beta Sheet'] == False:
            if row['C_gamma'] > 3.4:
                if (row['pH'] < 7.0 or row['pH'] == 'Both') and row['Temperature'] <= 100:
                    df.at[index, 'ptm'] = 'No'
                elif (row['pH'] >= 7.0 or row['pH'] == 'Both') and row['Temperature'] > 100:
                    if row['Motif'] in ['NG', 'NH', 'NS']:
                        df.at[index, 'ptm'] = 'Yes'
                    else:
                        df.at[index, 'ptm'] = 'No'
            elif row['C_gamma'] <= 3.4:
                if row['SASA'] <= 10:
                    if (row['pH'] < 7.0 or row['pH'] == 'Both') and row['Temperature'] <= 100:
                        df.at[index, 'ptm'] = 'No'
                    elif (row['pH'] >= 7.0 or row['pH'] == 'Both') and row['Temperature'] > 100:
                        df.at[index, 'ptm'] = 'Yes'
                elif row['SASA'] > 10:
                    if (row['pH'] < 7.0 or row['pH'] == 'Both') and row["Temperature"] < 100:
                        df.at[index, 'ptm'] = 'Yes'
                    elif (row["pH"] >= 7.0 or row['pH'] == 'Both') and row["Temperature"] >= 100:
                        df.at[index, 'ptm'] = 'Yes'
    columns_to_keep = ['PDB ID', 'Drug Name', 'Chain', 'Duplicate Number', 'Sequence', 'SASA', 'C_gamma', 'Res_Num',
                       'Is in Beta Sheet', 'Motif', 'ptm']
    new_df = df[columns_to_keep]
    new_df.to_csv(path_output, index=False)
    no_df = df[df['ptm'] == 'Yes']
    no_df.to_csv('no_ptm_train_2_arr.csv', index=False)


def do_file_without_duplecates(path_input, path_output):
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


def split_datatset(path_input, path_train, path_test):
    import csv
    import random

    with open(path_input, 'r') as input_file:
        reader = csv.DictReader(input_file)
        data = list(reader)

    random.shuffle(data)

    split_index = int(0.8 * len(data))
    train_data = data[:split_index]
    test_data = data[split_index:]

    with open(path_train, 'w') as train_file:
        writer = csv.DictWriter(train_file, fieldnames=reader.fieldnames, lineterminator='\n')
        writer.writeheader()
        writer.writerows(train_data)

    with open(path_test, 'w') as test_file:
        writer = csv.DictWriter(test_file, fieldnames=reader.fieldnames, lineterminator='\n')
        writer.writeheader()
        writer.writerows(test_data)


from Bio.PDB import PDBParser


def write_asn_residue_numbers(folder_path="ab_pdb", output_file="all_asn.csv"):
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


def get_asn_residue_numbers(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)
    return {residue.get_id()[1] for residue in structure.get_residues() if residue.get_resname() == 'ASN'}


def compare_folder_and_file(pdb_folder="ab_pdb", drug_file="pdb_to_drug.txt"):
    import os
    pdb_files = set(os.listdir(pdb_folder))

    with open(drug_file, 'r') as f:
        lines = f.readlines()

    with open(drug_file, 'w') as f:
        for line in lines:
            drug, pdb_id = line.strip().split()
            if pdb_id + '.pdb' in pdb_files:
                f.write(f'{drug} {pdb_id} yes\n')
            else:
                f.write(f'{drug} {pdb_id} no\n')


def make_dataset():
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




def find_tp_tn_fp_fn(path1, path2):
    import csv
    import os
    import warnings
    from Bio.PDB.PDBExceptions import PDBConstructionWarning

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)

        pdb_folder = 'D:\\Study\\нир\\ab_pdb'
        pdb_files = [os.path.join(pdb_folder, f) for f in os.listdir(pdb_folder) if f.endswith('.pdb')]
        tn = 0
        tp = 0
        fp = 0
        fn = 0

        with open(path1, 'r') as f:
            reader = csv.DictReader(f)
            file1_data = [row for row in reader]

        with open(path2, 'r') as f:
            reader = csv.DictReader(f)
            file2_data = [row for row in reader]

        for pdb_file in pdb_files:
            asn_numbers = get_asn_residue_numbers(pdb_file)
            for asn_number in asn_numbers:
                file1_row = next((row for row in file1_data if int(row['index']) - 1 == asn_number), None)
                file2_row = next((row for row in file2_data if int(row['Res_Num']) - 1 == asn_number), None)

                if (not file1_row or float(file1_row['% modified']) < 5) and not file2_row:
                    tn += 1
                elif file1_row and float(file1_row['% modified']) > 5 and not file2_row:
                    fn += 1
                elif (not file1_row or float(file1_row['% modified']) < 5) and file2_row:
                    fp += 1
                elif file1_row and float(file1_row['% modified']) > 5 and file2_row:
                    tp += 1

    recall = tp / (tp + fn) if (tp + fn) != 0 else 0
    precision = tp / (tp + fp) if (tp + fp) != 0 else 0
    accuracy = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) != 0 else 0

    print(f'TP: {tp}')
    print(f'FP: {fp}')
    print(f'FN: {fn}')
    print(f'TN: {tn}')
    print(f'Recall: {recall}')
    print(f'Precision: {precision}')
    print(f'Accuracy: {accuracy}')


def count_the_same_names():
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