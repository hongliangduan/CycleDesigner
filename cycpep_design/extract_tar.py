import os
from Bio import PDB

def extract_tar_chain(input_pdb):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("input_structure", input_pdb)

    # 假设只有两条链，获取所有链
    chains = list(structure[0].get_chains())

    if len(chains) != 2:
        raise ValueError("The PDB file must contain exactly two chains.")

    # 比较长度
    chain_1 = chains[0]
    chain_2 = chains[1]

    if len(chain_1) >= len(chain_2):
        longer_chain = chain_1
        shorter_chain = chain_2
    else:
        longer_chain = chain_2
        shorter_chain = chain_1

    # 靶点和肽长度，用于contigmap.contigs=...
    tar_len = len(longer_chain)
    pep_len = len(shorter_chain)
    print(f"Length of tar is {tar_len}, length of pep is {pep_len}")

    # 新建pdb，添加长链
    new_structure = PDB.Structure.Structure("target_structure")
    new_model = PDB.Model.Model(0)
    new_model.add(longer_chain)
    new_structure.add(new_model)

    # save
    file_prefix = os.path.basename(input_pdb)[:4]
    output_pdb = f"/home/yons/inputs/{file_prefix}_tar.pdb"
    io = PDB.PDBIO()
    io.set_structure(new_structure)
    io.save(output_pdb)
    print(f"Target extracted and saved to {output_pdb}")

# input_pdb = "/home/yons/lab/zch/cycpep_design/tar_pdb/input/3p8f_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_000.pdb"

# extract_tar_chain(input_pdb)
