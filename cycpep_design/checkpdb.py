import os

def check_subfolders_for_pdb_files(parent_folder):
    # 遍历主文件夹中的每个子文件夹
    for subfolder in os.listdir(parent_folder):
        subfolder_path = os.path.join(parent_folder, subfolder)
        
        # 确保它是一个文件夹
        if os.path.isdir(subfolder_path):
            contains_unrelaxed_rank_file = False
            
            # 检查该子文件夹中的所有文件
            for filename in os.listdir(subfolder_path):
                if "unrelaxed_rank" in filename and filename.endswith(".pdb"):
                    contains_unrelaxed_rank_file = True
                    break
            
            # 如果没有找到符合条件的文件，打印子文件夹名称
            if not contains_unrelaxed_rank_file:
                print(f"子文件夹 '{subfolder}' 中不包含 'unrelaxed_rank' 的 PDB 文件")

# 示例用法
parent_folder_path = "/home/yons/lab/zch/rfd2HF/outputdir/1104_H_25T"
check_subfolders_for_pdb_files(parent_folder_path)
