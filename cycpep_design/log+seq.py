import os
import pandas as pd

# 加载 CSV 文件
csv_file_path = '/home/yons/lab/zch/rfd2HF/outputdir/1104_H_25T/log_ipae_rmsd.csv'
data = pd.read_csv(csv_file_path)

# 生成标识列：protein_id_backbone_seqid
data['identifier'] = data['protein_id'].astype(str) + '_' + \
                     data['backbone'].astype(str) + '_' + \
                     data['seqid'].astype(str)

# 指定包含 .fasta 文件的文件夹路径
fasta_folder = '/home/yons/lab/zch/rfd2HF/inputdir/1104_H_25T'  # 替换为实际路径

# 用于存储序列的列表
sequences = []

# 遍历标识列，寻找对应的 .fasta 文件
for identifier in data['identifier']:
    fasta_file_path = os.path.join(fasta_folder, f"{identifier}.fasta")
    if os.path.exists(fasta_file_path):
        # 打开并读取 .fasta 文件的最后一行
        with open(fasta_file_path, 'r') as file:
            lines = file.readlines()
            # 获取最后一行，去掉换行符
            sequence = lines[-1].strip() if lines else ''
    else:
        # 如果文件不存在，记录为空值
        sequence = ''
    sequences.append(sequence)

# 添加 sequence 列到 DataFrame
data['sequence'] = sequences

# 保存更新后的 CSV 文件
output_csv_path = '/home/yons/lab/zch/rfd2HF/outputdir/1104_H_25T/log_ipae_rmsd-1120.csv'
data.to_csv(output_csv_path, index=False)

print(f"处理完成，文件已保存到：{output_csv_path}")
