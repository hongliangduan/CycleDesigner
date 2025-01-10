import os
import csv
import re
from pathlib import Path
from 计算基于位置对齐的CaRMSD import main

def find_native_file(native_dir, pdb_file_name):
    """
    根据 PDB ID 找对应的 native 文件。
    """
    # prefix = pdb_file_name[:4]
    prefix = pdb_file_name[:6]
    native_file = next(Path(native_dir).glob(f"{prefix}*.pdb"), None)
    return native_file

def extract_rank_number(pdb_file_name):
    """
    提取rank
    """
    match = re.search(r"_rank_00(\d+)", pdb_file_name)
    return match.group(1) if match else None

def process_files_and_write_csv(pdb_dir, native_dir, output_file):
    # pdb_dir = "/home/yons/lab/zch/rfd2HF/outputdir/1023/1023_30T"
    # native_dir = "/home/yons/inputs/pdbfixer/rfd_native/fixed"
    # output_file = "rmsd_results.csv"

    results = []

    # 遍历 PDB 文件路径中的所有子文件夹和 PDB 文件
    for root, _, files in os.walk(pdb_dir):
        for pdb_file in files:
            if pdb_file.endswith(".pdb"):
                pdb_file_path = os.path.join(root, pdb_file)
                rank_number = extract_rank_number(pdb_file)

                # 找到匹配的 native 文件
                native_file = find_native_file(native_dir, pdb_file)
                if native_file:
                    try:
                        # 计算 RMSD
                        rmsd_value = main(pdb_file_path, native_file)
                        results.append([pdb_file, rank_number,native_file.name, rmsd_value])
                    except Exception as e:
                        print(f"计算 {pdb_file} 和 {native_file.name} 时出错: {e}")

    # 将结果写入 CSV 文件
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PDB File", "Rank", "Native File", "RMSD"])
        writer.writerows(results)

    print(f"RMSD 计算完成，结果已写入 {output_file}")

if __name__ == "__main__":
    process_files_and_write_csv()
