import pandas as pd
import os
import shutil

def filter(csv_file, to_rosetta):
    file_path = csv_file
    output_path = to_rosetta
    os.makedirs(output_path, exist_ok=True)

    # read csv
    df = pd.read_csv(file_path)

    # 筛选符合条件的数据
    # filtered_df = df[(df['plddt'] > 80) & (df['iptm'] > 0.5) & (df['ipae'] < 0.35) & (df['rmsd'] < 3.5)]

    # 保存筛选后的数据到新 CSV 文件
    # filtered_df.to_csv('filtered_data.csv', index=False)

    # 筛选符合条件的数据
    filtered_df = df[(df['plddt'] > 80) & (df['iptm'] > 0.5) & (df['ipae'] < 0.35) & (df['rmsd'] < 3.5)]

    # 遍历筛选出的行并移动文件
    for _, row in filtered_df.iterrows():
        # 构建文件的完整路径
        src_path = os.path.join(
            row['folder'], f"{row['protein_id']}_{row['backbone']}_{row['seqid']}_unrelaxed_{row['file']}.pdb"
        )
        
        if os.path.exists(src_path):
            # 将文件移动到目标文件夹
            dest_path = os.path.join(output_path, os.path.basename(src_path))
            shutil.copy(src_path, dest_path)
            # print(f"已移动文件: {src_path} 到 {dest_path}")
        else:
            print(f"文件不存在，无法复制: {src_path}")

    return "done"
