import os
import re
import json
import numpy as np
import pandas as pd
from 计算基于位置对齐的CaRMSD import main

def extract_plddt_from_file(log_file, fasta_dir, native_dir):
    """
    从指定的日志文件中提取 PLDDT 和 IPAE 相关信息。
    """
    path = os.path.dirname(log_file)
    metrics = []
    
    with open(log_file, 'r') as fr:
        lines = fr.readlines()
        for line in lines:
            items = line.strip().split(' ')
            # 过滤满足 rank_00 开头的行
            if len(items) >= 3 and items[2].startswith('rank_00'):
                # 初始化记录
                m = [path, items[2], -1, -1, -1]
                
                # 提取 plddt, ptm 和 iptm
                if len(items) >= 4:
                    m[2] = items[3].split('=')[1]
                if len(items) >= 5:
                    m[3] = items[4].split('=')[1]
                if len(items) >= 6:
                    m[4] = items[5].split('=')[1]
                
                # 解析路径信息生成 protein_id 和其他字段
                protein_id = path.split('/')[-1].split('_')[0]
                backbone = path.split('_')[1]
                seqid = path.split('_')[2]
                rank = items[2].split('_')[1][2]
                model = items[2].split('_')[-3]
                
                # 添加 IPAE 信息提取
                json_file_path = f"{path}/{protein_id}_{backbone}_{seqid}_scores_rank_00{rank}_alphafold2_multimer_v3_model_{model}_seed_000.json"
                fasta_file_path = f"{fasta_dir}/{protein_id}_{backbone}_{seqid}.fasta"
                if os.path.exists(json_file_path):
                    with open(json_file_path) as f:
                        confidences = json.load(f)
                        pae_matrix = np.array(confidences['pae'])
                        with open(fasta_file_path) as fasta_file:
                            for fasta_line in fasta_file:
                                fasta_line = fasta_line.strip()
                                if fasta_line.endswith(":"):
                                    tar_len = len(fasta_line) - 1
                                    pae_interaction1 = np.mean(pae_matrix[:tar_len, tar_len:])
                                    pae_interaction2 = np.mean(pae_matrix[tar_len:, :tar_len])
                                    ipae = round((pae_interaction1 + pae_interaction2) / 2, 2)
                                    ipae = ipae/31
                                    # ipae = (np.mean(pae_matrix[:tar_len, tar_len:]) + np.mean(pae_matrix[tar_len:, :tar_len])) / 2
                                    m.append(ipae)

                # 添加rmsd计算
                pdb_file_path = f"{path}/{protein_id}_{backbone}_{seqid}_unrelaxed_rank_00{rank}_alphafold2_multimer_v3_model_{model}_seed_000.pdb"
                native_file_path = f"{native_dir}/{protein_id}_{backbone}.pdb"
                if os.path.exists(pdb_file_path) and os.path.exists(native_file_path):
                    rmsd = main(pdb_file_path, native_file_path)
                    m.append(rmsd)


                # # return native_file
                # match = re.search(r"_rank_00(\d+)", pdb_file_name)
                # for root, _, files in os.walk(pdb_dir):
                #     for pdb_file in files:
                #         if pdb_file.endswith(".pdb"):
                #             pdb_file_path = os.path.join(root, pdb_file)
                #             rank_number = extract_rank_number(pdb_file)

                #             # 找到匹配的 native 文件
                #             native_file = find_native_file(native_dir, pdb_file)
                #             if native_file:
                #                 try:
                #                     # 计算 RMSD
                #                     rmsd_value = main(pdb_file_path, native_file)
                #                     m.append([rmsd_value])
                #                 except Exception as e:
                #                     print(f"计算 {pdb_file} 和 {native_file.name} 时出错: {e}")

                
                # 合并信息并添加到结果列表
                m = [protein_id, backbone, seqid, rank] + m
                metrics.append(m)
    
    return metrics

def extract_plddt_from_folder(folder, fasta_dir, native_dir):
    """
    遍历指定文件夹中的所有日志文件，提取所有 PLDDT 和 IPAE 相关信息。
    """
    metrics = []
    
    for root, _, files in os.walk(folder):
        for f in files:
            if f == 'log.txt':
                fullname = os.path.join(root, f)
                metrics += extract_plddt_from_file(fullname, fasta_dir, native_dir)
    
    return metrics

def write_to_file(metrics, output_file):
    """
    将提取的 PLDDT 和 IPAE 信息写入 CSV 文件。
    """
    # Adjusted columns to match the structure of metrics
    df = pd.DataFrame(metrics, columns=[
        'protein_id', 'backbone', 'seqid', 'rank', 'folder', 'file', 'plddt', 'ptm', 'iptm', 'ipae', 'rmsd'
    ])
    df.to_csv(output_file, index=False)
    print(f"Data written to {output_file}")
    return "done"

