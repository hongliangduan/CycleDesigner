import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def fix_missing_residues_and_gaps(input_pdb, output_pdb):
    # 加载PDB文件
    fixer = PDBFixer(filename=input_pdb)

    # 查找缺失的残基并补齐
    # print("查找并修复缺失的氨基酸残基以及不连续的序列...")
    fixer.findMissingResidues()
    
    if fixer.missingResidues:
        print(f"missing residues: {fixer.missingResidues}")
    else:
        print("No missing residues found.")

    # 补齐所有缺失残基
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # 处理不连续的残基（gaps），如果存在，需要填补这些间隙
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    # 添加缺失的氢原子
    # print("添加缺失的氢原子...")
    fixer.addMissingHydrogens()

    # 保存修复后的PDB文件
    print(f"fixed pdb saved to{output_pdb}")
    with open(output_pdb, 'w') as out_file:
        PDBFile.writeFile(fixer.topology, fixer.positions, out_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        # print("用法: python fix_pdb.py input_pdb output_pdb")
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]

    fix_missing_residues_and_gaps(input_pdb, output_pdb)
