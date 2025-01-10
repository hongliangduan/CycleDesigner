from Bio.PDB import PDBParser
import numpy as np
import os

# parser
def get_structure(input_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)
    return structure

# calculate distance
def calc_distances(chain1, chain2, threshold=5.0):
    result = []
    contact_residues = set()
    
    for residue1 in chain1:
        if not residue1.has_id('CA'):
            continue
        for residue2 in chain2:
            if not residue2.has_id('CA'):
                continue

            # distance
            atom1 = residue1['CA']
            atom2 = residue2['CA']
            distance = atom1 - atom2

            if distance < threshold:
                result.append((residue1.get_id()[1], residue2.get_id()[1], distance))
                res_id = residue1.get_id()[1]
                contact_residues.add(f"{residue1}{res_id}")
                break  # 跳出 atom_short 循环
    
    return result, contact_residues

def main(input_pdb, chain_id1, chain_id2, threshold=5.0):
    structure = get_structure(input_pdb)
    chain1 = structure[0][chain_id1]
    chain2 = structure[0][chain_id2]

    distances, contact_residues = calc_distances(chain1, chain2, threshold)

    # 打印出距离小于阈值的残基序号和对应的距离
    for res1, res2, dist in distances:  # Unpacking the tuples correctly
        print(f"Chain {chain_id1} Residue {res1} - Chain {chain_id2} Residue {res2}: Distance = {dist:.2f} Å")

    print("Contact residues:", contact_residues)  # Optionally print contact residues


