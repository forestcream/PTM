from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB import DSSP
from Bio.PDB.Polypeptide import is_aa, Polypeptide
from Bio.PDB.vectors import calc_dihedral
import numpy as np
import csv
import math


def get_distance(atom1, atom2):
    coord1 = atom1.get_vector()
    coord2 = atom2.get_vector()
    return np.linalg.norm(coord1 - coord2)


data_dict = {}
with open('pdb_names_additional.txt', 'r') as f:
    for line in f:
        if not line.strip():
            continue
        try:
            drug, pdb_id = line.split()
            data_dict[pdb_id] = drug
        except ValueError:
            print(f"Error: Could not split line into drug and PDB ID: {line}")

ids = list(data_dict.keys())

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

with open('result_added.csv', 'a', newline='') as cvsfile:
    writer = csv.writer(cvsfile, delimiter=',')
    #writer.writerow(['PDB ID','Drug Name', 'Chain', 'Duplicate Number', 'Sequence', 'SASA', 'C_gamma', 'Res_Num', 'Is in Beta Sheet',
                     #'Motif', 'Phi', 'Psi', 'Chi1', 'Chi2', 'Hydrogen bonds to side chain'])
    for id in ids:

        dn = 1
        p = PDBParser(QUIET=1)
        struct = p.get_structure("id", f'ab_pdb/{id}.pdb')

        model = struct[0]
        dssp = DSSP(model, f'ab_pdb/{id}.pdb',
                    dssp="C:/Users/Owner/AppData/Local/Programs/Python/Python311/Scripts/dssp.exe")
        seq = []
        sequence_data = {}
        for chain in model:
            print("Chain ID")
            print(chain.id)
            sequence = ""
            for residue in chain:
                resname = residue.get_resname()
                if resname in d:
                    sequence += d[resname]

            seq.append(sequence)
            # duplicates of chains
            if sequence not in sequence_data:
                sequence_data[sequence] = dn
                dn += 1

            residues = list(chain)
            poly = Polypeptide(chain)
            phi_psi = poly.get_phi_psi_list()

            for i, residue in enumerate(residues):
                isInBetaSheet = False
                if residue.get_resname() == "ASN":
                    # residue number
                    res_num = residue.get_id()[1]
                    print("residue number")
                    print(res_num)

                    # is in beta sheet
                    res_data = dssp[(residue.get_parent().id, residue.id)]
                    sec_struct = res_data[2]
                    if sec_struct == 'E':
                        isInBetaSheet = True

                    # sasa
                    sr = ShrakeRupley()
                    sr.compute(struct, level="R")
                    sasa = (round(residue.sasa, 2))

                    # C_gamma
                    for atom in residue:
                        if atom.get_name() == "CG":
                            cg_atom = atom
                            next_residue = residues[i + 1]
                            n_atoms = [atom for atom in next_residue if atom.get_name() == "N"]
                            if n_atoms:
                                dist = get_distance(cg_atom, *n_atoms)
                                c_gamma = dist
                            else:
                                continue

                    # phi and psi angles
                    phi, psi = phi_psi[i]

                    # chi1 angle
                    if "N" in residue and "CA" in residue and "CB" in residue and "CG" in residue:
                        N = residue["N"].get_vector()
                        CA = residue["CA"].get_vector()
                        CB = residue["CB"].get_vector()
                        CG = residue["CG"].get_vector()
                        chi1 = calc_dihedral(N, CA, CB, CG)

                    # chi2 angle
                    if "CA" in residue and "CB" in residue and "CG" in residue and "OD1" in residue:
                        CA = residue["CA"].get_vector()
                        CB = residue["CB"].get_vector()
                        CG = residue["CG"].get_vector()
                        OD1 = residue["OD1"].get_vector()
                        chi2 = calc_dihedral(CA, CB, CG, OD1)

                    #Hydrogen bonding

                    hbond_count = 0
                    for atom in residue:
                        if atom.element == 'H' or atom.element == 'N':
                            for other_atom in residue:
                                if other_atom.element == 'O':
                                    distance = atom - other_atom
                                    if distance < 3:
                                        hbond_count += 1
                                        if hbond_count == 4:
                                            break

                    if next_residue.get_resname() in d and phi and psi and chi1 and chi2:
                        writer.writerow([id, data_dict[id],  chain.id, sequence_data[sequence], sequence, sasa, c_gamma, res_num,
                                     "Yes" if isInBetaSheet else "No", "N" + d[next_residue.get_resname()],
                                     phi * 180 / math.pi, psi * 180 / math.pi, chi1 * 180 / math.pi, chi2 * 180/ math.pi, hbond_count])
