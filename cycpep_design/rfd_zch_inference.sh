#!/bin/bash

#  run：./run_inference.sh input_pdb_file

if [ "$#" -ne 1 ]; then
    echo "用法：$0 input_pdb_file"
    exit 1
fi

INPUT_PDB=$1

# 提取链 ID
CHAIN_IDS=($(grep '^ATOM' $INPUT_PDB | awk '{print $5}' | sort | uniq))
NUM_CHAINS=${#CHAIN_IDS[@]}

if [ "$NUM_CHAINS" -ne 2 ]; then
    echo "错误：PDB 文件必须包含两条链。"
    exit 1
fi

CHAIN1=${CHAIN_IDS[0]}
CHAIN2=${CHAIN_IDS[1]}

# 获取链的长度和序列范围
get_chain_info() {
    local pdb_file=$1
    local chain_id=$2
    local res_nums=$(grep '^ATOM' $pdb_file | awk -v chain=$chain_id '$5 == chain {print $6}' | sort -n | uniq)
    local start_res=$(echo $res_nums | awk '{print $1}')
    local end_res=$(echo $res_nums | awk '{print $NF}')
    echo "$start_res-$end_res"
}

CHAIN1_RANGE=$(get_chain_info $INPUT_PDB $CHAIN1)
CHAIN2_RANGE=$(get_chain_info $INPUT_PDB $CHAIN2)

# 获取链长度
get_chain_length() {
    local pdb_file=$1
    local chain_id=$2
    local res_nums=$(grep '^ATOM' $pdb_file | awk -v chain=$chain_id '$5 == chain {print $6}' | sort -n | uniq)
    echo $res_nums | wc -w
}

LEN1=$(get_chain_length $INPUT_PDB $CHAIN1)
LEN2=$(get_chain_length $INPUT_PDB $CHAIN2)

# 确定较长和较短的链
LEN1=$(echo $CHAIN1_RANGE | awk -F'-' '{print $2-$1+1}')
LEN2=$(echo $CHAIN2_RANGE | awk -F'-' '{print $2-$1+1}')

if [ "$LEN1" -ge "$LEN2" ]; then
    LONG_CHAIN=$CHAIN1
    SHORT_CHAIN=$CHAIN2
    LONG_CHAIN_RANGE=$CHAIN1_RANGE
    SHORT_CHAIN_RANGE=$CHAIN2_RANGE
    TAR_LEN=$LEN1
    PEP_LEN=$LEN2
else
    LONG_CHAIN=$CHAIN2
    SHORT_CHAIN=$CHAIN1
    LONG_CHAIN_RANGE=$CHAIN2_RANGE
    SHORT_CHAIN_RANGE=$CHAIN1_RANGE
    TAR_LEN=$LEN2
    PEP_LEN=$LEN1
fi

echo "较长的链：$LONG_CHAIN (范围 $LONG_CHAIN_RANGE)"
echo "较短的链：$SHORT_CHAIN (范围 $SHORT_CHAIN_RANGE)"

# 提取靶点并保存为 xxxx_tar.pdb
PDB_BASENAME=$(basename $INPUT_PDB .pdb)
OUTPUT_PDB="/home/yons/inputs/hf_pre/single_tar/${PDB_BASENAME}_tar.pdb"

grep '^ATOM' $INPUT_PDB | awk -v chain=$LONG_CHAIN '$5 == chain' > $OUTPUT_PDB

sleep 2

# # 一对一hotspot
# cat > compute_contacts.py <<EOF
# # compute_contacts.py
# import sys
# from Bio.PDB import PDBParser
# import numpy as np

# def get_structure(pdb_file):
#     parser = PDBParser(QUIET=True)
#     return parser.get_structure('protein', pdb_file)

# def calc_distances(chain1, chain2, threshold=5.0):
#     contact_residues = set()
#     distances = []
#     for residue1 in chain1:
#         if not residue1.has_id('CA'):
#             continue
#         for residue2 in chain2:
#             if not residue2.has_id('CA'):
#                 continue
#             atom1 = residue1['CA'].coord
#             atom2 = residue2['CA'].coord
#             distance = np.linalg.norm(atom1 - atom2)
#             distances.append((distance, residue1.get_id()[1]))  # 记录距离和残基ID
#             if distance < threshold:
#                 contact_residues.add(residue1.get_id()[1])  # 仅添加长链的氨基酸序号

#         # 如果 contact_residues 中元素少于3个，按距离最小的排序，取前3个残基
#     if len(contact_residues) < 3:
#         sorted_distances = sorted(distances, key=lambda x: x[0])  # 按距离排序
#         top3_residues = [res_id for _, res_id in sorted_distances[:3]]  # 取前3个res_id
#         contact_residues.update(top3_residues)  # 更新 contact_residues
#     return contact_residues

# if __name__ == "__main__":
#     pdb_file = sys.argv[1]
#     chain_id1 = sys.argv[2]
#     chain_id2 = sys.argv[3]
#     threshold = float(sys.argv[4])

#     structure = get_structure(pdb_file)
#     chain1 = structure[0][chain_id1]
#     chain2 = structure[0][chain_id2]

#     contact_residues = calc_distances(chain1, chain2, threshold)
    
#     for res in contact_residues:
#         print(f"{chain_id1}{res}")  # 输出格式为 A1
# EOF

# 一对多hotspot
cat > compute_contacts.py <<EOF
import sys
from Bio.PDBp import *

pdb_file = sys.argv[1]
chain_long = sys.argv[2]
chain_short = sys.argv[3]
distance_threshold = 5.0

parser = PDBParser(QUIET=True)
structure = parser.get_structure('structure', pdb_file)

atoms_long = []
residues_long = []
atoms_short = []
residues_short = []

for model in structure:
    for chain in model:
        if chain.id == chain_long:
            residues_long.extend([res for res in chain if is_aa(res)])
            atoms_long.extend(res.get_atoms() for res in chain if is_aa(res))
        elif chain.id == chain_short:
            residues_short.extend([res for res in chain if is_aa(res)])
            atoms_short.extend(res.get_atoms() for res in chain if is_aa(res))

from itertools import product
contact_residues = set()
for res_long in residues_long:
    for res_short in residues_short:
        for atom_long, atom_short in product(res_long, res_short):
            distance = atom_long - atom_short
            if distance < distance_threshold:
                res_id = res_long.get_id()[1]
                contact_residues.add(res_id)
                break

for res_id in sorted(contact_residues):
    print(f"{chain_long}{res_id}")
EOF

CONTACT_RESIDUES=$(python compute_contacts.py $INPUT_PDB $LONG_CHAIN $SHORT_CHAIN 5.0)

# 格式化 ppi.hotspot_res=[]
HOTSPOT_RES="["
FIRST=1
for res in $CONTACT_RESIDUES; do
    if [ $FIRST -eq 1 ]; then
        HOTSPOT_RES+="$res"
        FIRST=0
    else
        HOTSPOT_RES+=",$res"
    fi
done
HOTSPOT_RES+="]"

echo "ppi.hotspot_res=$HOTSPOT_RES"

# 写入 contigmap.contigs，使用实际的残基编号范围
CONTIGMAP_CONTIGS="[${LONG_CHAIN}${LONG_CHAIN_RANGE}/0 ${PEP_LEN}-${PEP_LEN}]"

echo "contigmap.contigs=$CONTIGMAP_CONTIGS"

# 运行脚本
INFERENCE_SCRIPT=./scripts/run_inference.py
INFERENCE_MODEL_DIR=/home/yons/models
INFERENCE_INPUT_PDB=$OUTPUT_PDB
INFERENCE_OUTPUT_PREFIX=/home/yons/outputs/rfd_out/$(date +%Y%m%d)/$PDB_BASENAME/$PDB_BASENAME
INFERENCE_LOG=$PDB_BASENAME:_$(date +%Y%m%d).log

CMD="$INFERENCE_SCRIPT 'contigmap.contigs=$CONTIGMAP_CONTIGS' 'ppi.hotspot_res=$HOTSPOT_RES' diffuser.T=25 inference.deterministic=True inference.num_designs=5 inference.model_directory_path=$INFERENCE_MODEL_DIR inference.input_pdb=$INFERENCE_INPUT_PDB inference.output_prefix=$INFERENCE_OUTPUT_PREFIX > $INFERENCE_LOG 2>&1"

echo "运行命令："
echo $CMD
eval $CMD
