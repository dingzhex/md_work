# 初始结构获取 & 2D -> 3D 转换 (RDKit)
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdMolDescriptors
import os

workpath = 'outdata/comp_3D'
os.makedirs(f"{workpath}", exist_ok=True)

# 方式1: 从SMILES字符串创建 (自己绘制)
smiles = "O=C(NCC(N1CCC[C@H]1C#N)=O)C2=C3C(C=CC(OCCCN)=C3)=NC=C2"
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)  # 添加氢原子
'''
rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
print(f"总可旋转键数量: {rotatable_bonds}")

# 获取每个可旋转键的详细信息
for bond in mol.GetBonds():
    if bond.GetBondType() == Chem.BondType.SINGLE:
        if not bond.IsInRing():  # 排除环内键
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()

            # 排除酰胺键
            if not ((begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 7 and
                     begin_atom.GetDegree() == 3 and any(n.GetAtomicNum() == 8 for n in begin_atom.GetNeighbors())) or
                    (end_atom.GetAtomicNum() == 6 and begin_atom.GetAtomicNum() == 7 and
                     end_atom.GetDegree() == 3 and any(n.GetAtomicNum() == 8 for n in end_atom.GetNeighbors()))):

                # 排除终端键（如甲基）
                if begin_atom.GetDegree() > 1 and end_atom.GetDegree() > 1:
                    print(f"可旋转键: {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()}",
                          f"类型: {begin_atom.GetSymbol()}-{end_atom.GetSymbol()}")
'''
# 方式2: 从文件读取 (如 SDF, MOL)
# mol = Chem.MolFromMolFile('your_ligand.mol')
# mol = Chem.MolFromMolBlock(mol_block)

# 生成3D构象 (使用ETKDGv3方法)
params = AllChem.ETKDGv3()
params.randomSeed = 0xf00d  # 可选的随机种子，确保可重复性
params.numThreads = 0  # 0表示使用所有可用核心
params.useSmallRingTorsions = True  # 改进小环构象采样
params.useMacrocycleTorsions = True  # 大环采样优化（若适用）
params.enforceChirality = True  # 确保手性正确
conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=200, params=params)  # 100-200个构象

# 使用MMFF94优化所有构象并选择最优  UFF力场不如MMFF94精确
min_energy = float('inf')
energies = []
for conf_id in conf_ids:
    # 使用MMFF94替代UFF
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    ff.Minimize(maxIts=500, energyTol=1e-5)
    energy = ff.CalcEnergy()
    # ff.CalcEnergy()
    energies.append(energy)

# 3. 构象聚类分析（在优化后！）
# 计算所有构象对之间的RMSD（忽略氢原子）
rmsd_matrix = []
heavy_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
for i in range(len(conf_ids)):
    for j in range(i + 1, len(conf_ids)):
        rmsd = AllChem.GetConformerRMS(
            mol, conf_ids[i], conf_ids[j],
            prealigned=True  # 自动对齐
        )
        rmsd_matrix.append(rmsd)

from rdkit.ML.Cluster import Butina
import numpy as np

# Butina聚类

dist_threshold = 1.5  # 紫杉醇推荐阈值(Å)

for dist_threshold in np.arange(1, 10 + 0.2 / 2, 0.2):

    clusters = Butina.ClusterData(
        rmsd_matrix,
        len(conf_ids),
        dist_threshold,
        isDistData=True
    )
    if len(clusters) <= 10:
        break
# 4. 选择每个簇的代表构象
representative_confs = []
representative_confs_energies = []
for cluster in clusters:
    # 取簇内能量最低的构象
    cluster_energies = [energies[i] for i in cluster]
    min_energy_idx = cluster[np.argmin(cluster_energies)]
    representative_confs.append(conf_ids[min_energy_idx])
    representative_confs_energies.append(min(cluster_energies))

print(f"从{len(conf_ids)}个构象中选出{len(representative_confs)}个代表构象")

for min_conf_id in representative_confs:
    # 获取最低能量构象
    min_conf = mol.GetConformer(min_conf_id)
    print(energies[min_conf_id])

    # 创建分子的深拷贝（包含原子和键信息）
    new_mol = Chem.Mol(mol)
    # 移除新分子中的所有构象
    new_mol.RemoveAllConformers()
    # 将最低能量构象添加到新分子中
    new_mol.AddConformer(min_conf)

    import os

    os.makedirs(f'outdata/comp_3D/ligand_3d_initial_{min_conf_id}', exist_ok=True)
    # 保存初步3D结构
    Chem.MolToMolFile(new_mol, f'outdata/comp_3D/ligand_3d_initial_{min_conf_id}/ligand_3d_initial_{min_conf_id}_.mol')

import pandas as pd

# pd.DataFrame({'liga_energies_sort_id':representative_confs}).to_csv('outdata/comp_3D/energies_id.csv', index=False)

pd.DataFrame({'liga_energies_sort_id': representative_confs, 'energy': representative_confs_energies}).sort_values(
    'energy').to_csv('outdata/comp_3D/energies_id.csv', index=False)

import subprocess

# 电荷分配与质子化状态优化
for file in os.listdir(f'outdata/comp_3D'):

    if 'ligand' not in file:
        continue

    tmp_file = os.listdir(os.path.join(f'outdata/comp_3D', file))
    input_path = tmp_file[0]
    output_path = 'obabel_ligand.mol2'
    # 构建命令
    command = [
        "obabel",
        input_path,  # 使用文件名而不是完整路径
        "-O", output_path,
        "-p", str(7.4),
        # '--gen3d',
        '--fast',
        '--minimize',
        '--steps 50',  # 50
        ' --crit 1e-6',  # 1e-6
        '--ff', "UFF",  # 用快速力场MMFF94s替代默认MMFF94  #UFF速度更快
    ]
    # 运行命令
    result = subprocess.run(
        command,
        cwd=str(os.path.join(f'outdata/comp_3D', file)),  # 指定工作目录
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding='utf-8'
    )

    # 使用'-c', 'bcc'重新计算更精确的电荷  （AM1-BCC）    或者不变  '-c', 'bcc'
    command = [
        'antechamber',
        '-i', 'obabel_ligand.mol2',
        '-fi', 'mol2',
        '-o', 'ligand.prepi',
        '-fo', 'prepi', '-c', 'bcc'  # , '-nc', '0'
    ]
    result = subprocess.run(
        command,
        cwd=str(os.path.join(f'outdata/comp_3D', file)),  # 指定工作目录
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding='utf-8'
    )
