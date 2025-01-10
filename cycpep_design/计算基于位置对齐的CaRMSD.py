from Bio import PDB
import numpy as np

def get_shorter_chain(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    chains = list(structure.get_chains())
    chain_lengths = [(chain, len(list(chain.get_residues()))) for chain in chains]
    shorter_chain = min(chain_lengths, key=lambda x: x[1])[0]
    return shorter_chain

def get_backbone_coords(chain):
    backbone_atoms = ['N', 'CA', 'C']
    coords = []
    residue_ids = []
    for residue in chain:
        if PDB.is_aa(residue):
            residue_coords = []
            for atom_name in backbone_atoms:
                try:
                    atom = residue[atom_name]
                    residue_coords.append(atom.get_coord())
                except KeyError:
                    residue_coords = []
                    break
            if residue_coords:
                coords.extend(residue_coords)
                residue_ids.append(residue.get_id())
    return np.array(coords), residue_ids

def get_residue_id_to_indices(residue_ids):
    residue_id_to_indices = {}
    for i, residue_id in enumerate(residue_ids):
        start_idx = i * 3
        end_idx = start_idx + 3
        residue_id_to_indices[residue_id] = (start_idx, end_idx)
    return residue_id_to_indices

def kabsch(P, Q):
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = np.linalg.det(np.dot(W.T, V.T))
    D = np.identity(3)
    D[2,2] = d
    U = np.dot(np.dot(W.T, D), V.T)
    return U

def superimpose_and_calculate_rmsd(coords1, coords2):
    P_centroid = np.mean(coords1, axis=0)
    Q_centroid = np.mean(coords2, axis=0)
    P_centered = coords1 - P_centroid
    Q_centered = coords2 - Q_centroid
    U = kabsch(P_centered, Q_centered)
    Q_rotated = np.dot(Q_centered, U)
    Q_superposed = Q_rotated + P_centroid
    diff = coords1 - Q_superposed
    rmsd = np.sqrt(np.sum(np.square(diff)) / len(coords1))
    return rmsd

def main(pdb_file1, pdb_file2):
    chain1 = get_shorter_chain(pdb_file1)
    chain2 = get_shorter_chain(pdb_file2)
    coords1, residue_ids1 = get_backbone_coords(chain1)
    coords2, residue_ids2 = get_backbone_coords(chain2)
    residue_id_to_indices1 = get_residue_id_to_indices(residue_ids1)
    residue_id_to_indices2 = get_residue_id_to_indices(residue_ids2)
    common_residue_ids = set(residue_ids1).intersection(set(residue_ids2))
    if not common_residue_ids:
        # print("No common residues found between the two chains.")
        return '/'
    coords1_common = []
    coords2_common = []
    for residue_id in common_residue_ids:
        idx1 = residue_id_to_indices1[residue_id]
        idx2 = residue_id_to_indices2[residue_id]
        coords1_common.append(coords1[idx1[0]:idx1[1]])
        coords2_common.append(coords2[idx2[0]:idx2[1]])
    coords1_common = np.concatenate(coords1_common, axis=0)
    coords2_common = np.concatenate(coords2_common, axis=0)
    rmsd = superimpose_and_calculate_rmsd(coords1_common, coords2_common)
    # print(f"Backbone RMSD between the shorter chains: {rmsd:.3f} Å")
    return rmsd

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print("Usage: python script.py pdb_file1 pdb_file2")
        sys.exit(1)
    pdb_file1 = sys.argv[1]
    pdb_file2 = sys.argv[2]
    main(pdb_file1, pdb_file2)
