def make_not_good_seq():
    import csv

    heavy_chains = {}
    with open('heavy_chains.txt', 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            heavy_chains[row[1]] = row[0]

    light_chains = {}
    with open('light_chains.txt', 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            light_chains[row[1]] = row[0]

    tmp_chains = {}
    with open('tmp.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            drug_name = row['name']
            if drug_name not in tmp_chains:
                tmp_chains[drug_name] = []
            tmp_chains[drug_name].append(row)

    with open('result.csv', 'r') as results_file, open('bad_seq.csv', 'w') as bad_seq_file:
        results_reader = csv.DictReader(results_file)
        fieldnames = ['PDB ID', 'Drug Name', 'Chain', 'Duplicate Number', 'Sequence', 'Chain from tmp', 'Indexes',
                      'Res_num', 'CDR from tmp']
        bad_seq_writer = csv.DictWriter(bad_seq_file, fieldnames=fieldnames)
        bad_seq_writer.writeheader()

        for row in results_reader:
            drug_name = row['Drug Name']
            sequence = row['Sequence']

            if drug_name in heavy_chains:
                if not sequence.startswith(heavy_chains[drug_name]):
                    new_row = {
                        'PDB ID': row['PDB ID'],
                        'Drug Name': drug_name,
                        'Chain': row['Chain'],
                        'Duplicate Number': row['Duplicate Number'],
                        'Sequence': sequence,
                        'Chain from tmp': heavy_chains[drug_name],
                        'Res_num': row['Res_Num']
                    }
                    if drug_name in tmp_chains:
                        indexes = [tmp_row['index'] for tmp_row in tmp_chains[drug_name]]
                        new_row['Indexes'] = ','.join(indexes)
                        cdrs = [tmp_row['CDR'] for tmp_row in tmp_chains[drug_name]]
                        new_row['CDR from tmp'] = ','.join(cdrs)
                    else:
                        new_row['Indexes'] = 'Not in tmp'
                        new_row['CDR from tmp'] = 'Not in tmp'
                    bad_seq_writer.writerow(new_row)
            elif drug_name in light_chains:
                if not sequence.startswith(light_chains[drug_name]):
                    new_row = {
                        'PDB ID': row['PDB ID'],
                        'Drug Name': drug_name,
                        'Chain': row['Chain'],
                        'Duplicate Number': row['Duplicate Number'],
                        'Sequence': sequence,
                        'Chain from tmp': light_chains[drug_name],
                        'Res_num': row['Res_Num']
                    }
                    if drug_name in tmp_chains:
                        indexes = [tmp_row['index'] for tmp_row in tmp_chains[drug_name]]
                        new_row['Indexes'] = ','.join(indexes)
                        cdrs = [tmp_row['CDR'] for tmp_row in tmp_chains[drug_name]]
                        new_row['CDR from tmp'] = ','.join(cdrs)
                    else:
                        new_row['Indexes'] = 'Not in tmp'
                        new_row['CDR from tmp'] = 'Not in tmp'
                    bad_seq_writer.writerow(new_row)


# make_not_good_seq()
import pandas as pd


# df = pd.read_csv('bad_seq.csv')
# df['light_chain'] = ''
# df['heavy_chain'] = ''
#
# for i, row in df.iterrows():
#     if row['Indexes'] != 'not in tmp':
#         indexes = row['Indexes'].split(',')
#         if row['CDR from tmp'] != 'not in tmp':
#             chains = row['CDR from tmp'].split(',')
#             light_chain_indexes = []
#             heavy_chain_indexes = []
#             for chain, index in zip(chains, indexes):
#                 if chain.startswith('L'):
#                     light_chain_indexes.append(index)
#                 elif chain.startswith('H'):
#                     heavy_chain_indexes.append(index)
#             df.at[i, 'light_chain'] = ', '.join(light_chain_indexes)
#             df.at[i, 'heavy_chain'] = ', '.join(heavy_chain_indexes)
#
# print(df)
# df.to_csv('bad_seq_changed.csv', index=False)

def calculate_updated_indexes(row):
    indices_rows = updated_indices[
        (updated_indices['seq'] == row['Sequence'])]
    print(indices_rows['start_ab_in_bad_seq'])
    if len(indices_rows) == 0:
        return None

    indices_row = indices_rows.iloc[0]

    return int(row['Res_Num']) + int(indices_row['how_much_cut']) - int(indices_rows['start_ab_in_bad_seq'].iloc[0])


    # if indices_row['new_light'] == 'NOTCHAIN' or indices_row['new_heavy'] == 'NOTCHAIN':
    #     return 'NOTCHAIN'
    # elif indices_row['new_light'] == 'same' or indices_row['new_heavy'] == 'same':
    #     return row['Res_Num']
    # elif pd.isna(indices_row['new_light']) or pd.isna(indices_row['new_heavy']):
    #     if pd.isna(indices_row['how_much_cut']) or indices_row['how_much_cut'] == '-':
    #         return "-"
    #     else:
    #         return row['Res_Num'] + int(indices_row['how_much_cut'])
    # else:
    #     return row['Res_Num'] + int(indices_row['how_much_cut'])


import pandas as pd

result = pd.read_csv('result.csv')
updated_indices = pd.read_csv('updated_indices_190623 (2).csv')

result['updated indexes'] = result.apply(calculate_updated_indexes, axis=1)
###!!!!!!!!!!!!!!!!!!!!!
result.loc[result['updated indexes'].isna(), 'updated indexes'] = result['Res_Num']

result.to_csv('updated_file.csv', index=False)
#
# import pandas as pd
#
# your_file = pd.read_csv('result.csv')
#
# other_file = pd.read_csv('updated_indices.csv')
#
# merged = pd.merge(your_file, other_file[['name', 'seq', 'new_light', 'new_heavy']], left_on=['Drug Name', 'Sequence'], right_on=['name', 'seq'], how='left')
#
# merged['is in pdb same light chain'] = merged.apply(lambda row: 'yes' if row['new_light'] == 'same' else 'no', axis=1)
# merged['is in pdb same heavy chain'] = merged.apply(lambda row: 'yes' if row['new_heavy'] == 'same' else 'no', axis=1)
#
#
# final = merged[['PDB ID', 'Drug Name', 'Chain', 'Duplicate Number', 'Sequence', 'is in pdb same light chain', 'is in pdb same heavy chain']].copy()
# final = final.drop_duplicates(subset=['Drug Name', 'Chain'])
# final.to_csv('chain_of_ab.csv', index=False)
