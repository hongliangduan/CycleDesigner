import json
import os
import numpy as np

def fasta_read(fasta_file):
    with open(fasta_file, 'r') as f:
        list_info = f.readlines()
        peptide_length = len(list_info[1].strip())  # 计算肽链的长度
        return peptide_length

def interface_plddt_read(interface_file):
    with open(interface_file, 'r') as f:
        list_info = f.readlines()
        chainA = [int(item) - 1 for item in list_info[0].strip().split()]  # 将残基编号减1
        chainB = [int(item) - 1 for item in list_info[1].strip().split()]  # 将残基编号减1
        return chainA, chainB

# 设置路径
path_fasta = '/home/yons/lab/zch/rfd2HF/inputdir/1104_H_25T/1sfi_0_1.fasta'
path_json = '/home/yons/lab/zch/rfd2HF/outputdir/1104_H_25T/1sfi_0_1/1sfi_0_1_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json'
path_interface = '/path/to/interface_directory'

# 结果保存列表
result_data = []

for fasta_file in os.listdir(path_fasta):
    # 读取肽链长度和接口信息
    peptide_length = fasta_read(os.path.join(path_fasta, fasta_file))
    info = fasta_file.split('.')[0]
    chainA, chainB = interface_plddt_read(os.path.join(path_interface, f"{info.split('_')[0]}.txt"))
    
    # 读取JSON文件中的PAE信息
    json_file = os.path.join(path_json, info, f'rank_001.json')  # 假设使用rank_001的模型
    with open(json_file) as f:
        confidences = json.load(f)
        pae_matrix = np.array(confidences['pae'])

        # 计算interface_PAE
        peptide_receptor_pae = sum(pae_matrix[y][x] for y in chainA for x in chainB)
        receptor_peptide_pae = sum(pae_matrix[y][x] for y in chainB for x in chainA)
        num = len(chainA) * len(chainB)
        total_interface_pae = round((peptide_receptor_pae / num + receptor_peptide_pae / num) / 2, 2)

        result_data.append([info, total_interface_pae])

# 输出结果
for data in result_data:
    print(f"Peptide Name: {data[0]}, Interface PAE: {data[1]}")




import json
import numpy as np

# 读取 JSON 文件
json_file = "your_file.json"
with open(json_file) as f:
    confidences = json.load(f)
    pae_matrix = np.array(confidences['pae'])

# 设置序列长度
pep_len = ...  # 设置pep_len的值
tar_len = ...  # 设置tar_len的值

# 计算ipae
ipae = (np.mean(pae_matrix[:tar_len, tar_len:]) + np.mean(pae_matrix[tar_len:, :tar_len])) / 2

print("ipae:", ipae)
