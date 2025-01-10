#!/bin/bash

# 默认值为空
INPUT_PDB=""

# 解析命令行选项
while getopts ":p:" opt; do
  case $opt in
    p)
      INPUT_PDB=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# 检查是否提供了 PDB 文件
if [ -z "$INPUT_PDB" ]; then
  echo "Usage: $0 -p <input_pdb_file>"
  exit 1
fi

BASENAME=$(basename "$INPUT_PDB" .pdb)
CLEAN_PDB="${BASENAME}_clean.pdb"
PROTONATED_PDB="${BASENAME}_protonated.pdb"
DISULFIDE_PAIRS="${BASENAME}_disulfide_pairs.txt"
TLEAP_INPUT="${BASENAME}_tleap.in"
PRMTOP_FILE="com_sol.prmtop"
INPCRD_FILE="com_sol.inpcrd"

# Step 1: 使用 pdb4amber 进行初步处理
echo "Running pdb4amber to clean the PDB file..."
pdb4amber -i $INPUT_PDB -o $CLEAN_PDB -y

# Step 2: 使用 pdb2pqr 进行质子化处理
echo "Running pdb2pqr to protonate the PDB file..."
pdb2pqr --ff=AMBER --keep-chain --ffout=AMBER $CLEAN_PDB $PROTONATED_PDB

# Step 3: 提取 CYX 残基中的 SG 原子
echo "Extracting CYX SG atoms from the protonated PDB file..."
grep "CYX" $PROTONATED_PDB | grep " SG " > cyx_sg_atoms.txt

# Step 4: 调用 Python 脚本计算二硫键
echo "Calculating SG-SG distances using Python..."
python3 calculate_disulfide_bonds.py cyx_sg_atoms.txt $DISULFIDE_PAIRS

# Step 5: 自动检测第二条链的残基编号范围
echo "Detecting the residue range of the second chain..."
SECOND_CHAIN=$(grep "ATOM" $PROTONATED_PDB | awk '{print $5}' | uniq | sed -n '2p') # 提取第二条链的残基编号
FIRST_RES=$(grep "ATOM" $PROTONATED_PDB | grep " $SECOND_CHAIN " | awk '{print $6}' | head -n 1)
LAST_RES=$(grep "ATOM" $PROTONATED_PDB | grep " $SECOND_CHAIN " | awk '{print $6}' | tail -n 1)
RESIDUE_RANGE=":${FIRST_RES}-${LAST_RES}"

echo "Second chain residue range detected: $RESIDUE_RANGE"

# Step 6: 生成 tleap 输入文件，定义二硫键
echo "Generating tleap input file..."
cat > $TLEAP_INPUT << EOF
# 加载力场和溶剂模型，默认ff19SB+opc
source leaprc.protein.ff19SB
source leaprc.water.opc

# 加载 PDB 文件
mol = loadpdb ${PROTONATED_PDB}

# 定义二硫键
EOF

# 读取并写入二硫键对
while read -r res1 res2; do
  echo "bond mol.${res1}.SG mol.${res2}.SG" >> $TLEAP_INPUT
done < $DISULFIDE_PAIRS

# 完成 tleap 输入文件
cat >> $TLEAP_INPUT << EOF
solvateoct mol OPCBOX 10.0
addions mol Na+ 0
addions mol Cl- 0
set default PBRadii mbondi2
saveamberparm mol ${PRMTOP_FILE} ${INPCRD_FILE}
quit
EOF

# Step 7: 使用 tleap 生成拓扑和坐标文件
echo "Running tleap to generate topology and coordinate files..."
tleap -f $TLEAP_INPUT

# Step 8: 使用 ante-MMPBSA 剥离溶剂和生成复合物、受体、配体的无水拓扑文件
echo "Running ante-MMPBSA.py to generate no-water topology files for complex, receptor, and ligand..."
ante-MMPBSA.py -p ${PRMTOP_FILE} \
               -c com_nowat.parm7 \
               -r rep_nowat.parm7 \
               -l lig_nowat.parm7 \
               -s ":WAT,Na+,Cl-" \
               -n "$RESIDUE_RANGE"

echo "Topology files generated:"
echo " - Complex without water: com_nowat.parm7"
echo " - Receptor without water: rep_nowat.parm7"
echo " - Ligand without water: lig_nowat.parm7"
