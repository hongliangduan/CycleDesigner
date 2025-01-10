import math
import sys

# 从命令行获取输入和输出文件
cyx_sg_atoms_path = sys.argv[1]
disulfide_pairs_path = sys.argv[2]

cyx_sg_atoms = []

# 解析 cyx_sg_atoms.txt 文件，提取残基编号和坐标
with open(cyx_sg_atoms_path, 'r') as file:
    for line in file:
        if line.startswith("ATOM"):
            residue_number = line[22:26].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            cyx_sg_atoms.append((residue_number, x, y, z))

# 计算 SG 原子之间的距离
threshold = 2.5
disulfide_pairs = []

for i in range(len(cyx_sg_atoms)):
    for j in range(i + 1, len(cyx_sg_atoms)):
        res1, x1, y1, z1 = cyx_sg_atoms[i]
        res2, x2, y2, z2 = cyx_sg_atoms[j]
        distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
        if distance < threshold:
            disulfide_pairs.append((res1, res2))

# 将计算结果写入 disulfide_pairs.txt 文件
with open(disulfide_pairs_path, "w") as outfile:
    for res1, res2 in disulfide_pairs:
        outfile.write(f"{res1} {res2}\n")
