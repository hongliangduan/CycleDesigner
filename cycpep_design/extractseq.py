import os
from Bio.PDB import PDBParser, Polypeptide
from Bio.SeqUtils import seq1

def extract_sequences_from_pdb(file_path):
    """从PDB文件提取所有链的序列"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(file_path).split('.')[0], file_path)
    sequences = []

    for model in structure:
        for chain in model:
            sequence = []
            for residue in chain:
                if Polypeptide.is_aa(residue, standard=True):  # 检查是否为标准氨基酸
                    aa = seq1(residue.resname)  # 转换三字母为一字母代码
                    sequence.append(aa)
            if sequence:
                sequences.append("".join(sequence))

    # 在第一条链的序列末尾添加符号 ":"
    if sequences:
        sequences[0] += ":"

    return sequences

def write_fasta(pdb_id, sequences, output_folder):
    """将提取的序列写入FASTA文件"""
    fasta_content = f">{pdb_id}\n" + "\n".join(sequences)
    output_path = os.path.join(output_folder, f"{pdb_id}.fasta")
    with open(output_path, "w") as f:
        f.write(fasta_content)

def main(input_folder, output_folder):
    """遍历PDB文件，提取序列并写入FASTA文件"""
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith(".pdb"):
            pdb_path = os.path.join(input_folder, filename)
            pdb_id = os.path.splitext(filename)[0]
            sequences = extract_sequences_from_pdb(pdb_path)
            write_fasta(pdb_id, sequences, output_folder)

if __name__ == "__main__":
    input_folder = "/home/yons/inputs/pdbfixer/rfd_native/fixed"
    output_folder = "/home/yons/inputs/pdbfixer/rfd_native/fixed/fasta"
    main(input_folder, output_folder)
