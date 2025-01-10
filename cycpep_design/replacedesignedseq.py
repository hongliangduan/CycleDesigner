import os
import re

def extract_pdbid_from_filename(fa_file):
    """从FA文件名中提取pdbid"""
    filename = os.path.basename(fa_file)
    pdbid = filename.split('_')[0]  # 提取'1sfi'部分
    return pdbid

def read_sequences_from_fa(fa_file):
    """读取FA文件中的所有有效序列"""
    sequences = []
    with open(fa_file, 'r') as f:
        content = f.read()
        # 匹配纯字母序列行
        all_sequences = re.findall(r'^[A-Z]+$', content, re.MULTILINE)
        
        # 过滤掉全为G的序列
        for seq in all_sequences:
            if set(seq) != {'G'}:  # 如果序列不全为G，则为有效序列
                sequences.append(seq)

    if len(sequences) < 5:
        raise ValueError(f"{fa_file} 中有效序列不足五个")
    return sequences[:5]

def read_fasta_sequences(fasta_path):
    """读取FASTA文件中的链序列"""
    sequences = []
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
        # 提取不以'>'开头的序列行
        for line in lines:
            if not line.startswith(">"):
                sequences.append(line.strip())
    return sequences

def replace_chain_and_save(pdb_id, original_sequences, new_sequence, output_folder, index, fa_index):
    """替换第二条链的序列并保存为新的FASTA文件"""
    modified_sequences = original_sequences.copy()
    if len(modified_sequences) > 1:
        modified_sequences[1] = new_sequence  # 替换第二条链

    # 生成新的FASTA文件内容
    fasta_content = f">{pdb_id}_{fa_index}_{index}\n" + "\n".join(modified_sequences)
    output_path = os.path.join(output_folder, f"{pdb_id}_{fa_index}_{index}.fasta")
    with open(output_path, "w") as f:
        f.write(fasta_content)
    # print(f"生成新文件: {output_path}")

def process_fa_file(fa_file, fasta_folder, output_folder, fa_index):
    """处理单个FA文件并替换匹配FASTA文件中的第二条链"""
    pdb_id = extract_pdbid_from_filename(fa_file)
    new_sequences = read_sequences_from_fa(fa_file)

    matched_file = f"{pdb_id}.fasta"
    fasta_path = os.path.join(fasta_folder, matched_file)

    if os.path.exists(fasta_path):
        # 读取原始FASTA文件中的链序列
        original_sequences = read_fasta_sequences(fasta_path)

        # 用五个新序列逐一替换第二条链并保存
        for i, new_seq in enumerate(new_sequences, start=1):
            replace_chain_and_save(pdb_id, original_sequences, new_seq, output_folder, i, fa_index)
    else:
        print(f"未找到匹配的FASTA文件: {matched_file}")

def main(fa_folder, fasta_folder):
    """主函数，遍历FA文件并处理匹配的FASTA文件"""
    output_folder = os.path.join(fa_folder, "new")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for filename in os.listdir(fa_folder):
        if re.match(r".+_\d\.fa$", filename):  # 匹配以数字结尾的.fa文件
            fa_file = os.path.join(fa_folder, filename)
            fa_index = filename.split('_')[-1].split('.')[0]  # 提取最后的数字部分
            process_fa_file(fa_file, fasta_folder, output_folder, fa_index)

if __name__ == "__main__":
    fa_folder = "/home/yons/inputs/pdbfixer/rfd_native/fixed/fasta/241023_50T"
    fasta_folder = "/home/yons/inputs/pdbfixer/rfd_native/fixed/fasta"
    main(fa_folder, fasta_folder)
