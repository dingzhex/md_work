# 加氢,能量最小化


from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm.unit import picoseconds, kilojoule, mole, nanometer

# 1. 读取初始结构（假设为无水化或氢原子缺失的PDB）
input_pdb = 'outdata/aa_site_randsearch_/protein.pdb'
# input_pdb = 'outdata/protein/new_AF-P31751-F1-model_v4.pdb'
# input_pdb = 'outdata/protein/new_amb_AF-P31751-F1-model_v4.pdb'
pdb = PDBFile(input_pdb)  # 替换为你的输入文件

# 2. 使用Modeller添加氢原子
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')  # 选择力场（可替换）



modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield,pH=7.4)  # 根据力场添加氢原子
'''
# 3. 创建系统（隐式溶剂模型，适合快速最小化）
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=NoCutoff,  # 无截断（最小化常用）
    constraints=HBonds        # 约束所有键（加速最小化）
)

'''

# 创建带显式溶剂的系统
modeller.addSolvent(
    forcefield,
    padding=1.0*nanometers,   # 水层厚度
    ionicStrength=0.15*molar   # 离子浓度
)
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,       # 使用长程静电
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds
)



# 4. 设置能量最小化
integrator = VerletIntegrator(0.001*picoseconds)  # 极小时间步（仅最小化不需要积分）
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# 5. 执行最小化
print("Minimizing energy...")
simulation.minimizeEnergy(
    maxIterations=1000,       # 最大迭代次数
    tolerance=10.0*kilojoule_per_mole/nanometer  # 修正为梯度容差
)

# 6. 获取并保存最小化后的结构
minimized_positions = simulation.context.getState(
    getPositions=True
).getPositions()


os.makedirs(os.path.dirname('outdata/aa_site_randsearch_/minimized_all.pdb'),exist_ok=True)
PDBFile.writeFile(
    simulation.topology,
    minimized_positions,
    open('outdata/aa_site_randsearch_/minimized_all.pdb', 'w')
)

# 更精细 可以去除水分子和离子(或部分)， 并在下面 去除-d 参数

print("Done! Saved minimized structure to minimized.pdb")

# 在模拟后移除溶剂和离子
from simtk.openmm.app import PDBFile, Modeller
# 1. 加载溶剂化系统
pdb = PDBFile('outdata/aa_site_randsearch_/minimized_all.pdb')

# 2. 识别并移除水和离子
protein_atoms = []
for residue in pdb.topology.residues():
    # print(residue.name)
    if residue.name in ['HOH', 'SOL', 'WAT', 'NA', 'CL']:  # 保留蛋白质
        # print(residue.atoms())
        # print([a.index for a in residue.atoms()])
        protein_atoms.extend(residue.atoms())
        # protein_atoms.extend(a.index for a in residue.atoms())

# 3. 创建仅蛋白质的模型
modeller_protein = Modeller(pdb.topology, pdb.positions)
modeller_protein.delete(protein_atoms)  # 删除非蛋白部分

# 4. 保存纯净蛋白质
with open('outdata/aa_site_randsearch_/minimized.pdb', 'w') as f:
    PDBFile.writeFile(modeller_protein.topology, modeller_protein.positions, f)


###
from openmm.app import *
from openmm import *
from openmm.unit import *
import subprocess

# 2. 使用OpenBabel转换PDB到PDBQT
def convert_pdb_to_pdbqt(input_pdb, output_pdbqt):
    command = [
        'obabel',
        '-ipdb', input_pdb,
        '-opdbqt',
        '-xr',  # 移除非极性氢（对接标准）
        # '-d',
        '-O', output_pdbqt
    ]
    subprocess.run(command, check=True)

convert_pdb_to_pdbqt('outdata/aa_site_randsearch_/minimized.pdb', 'outdata/aa_site_randsearch_/output.pdbqt')

# obabel -ipdb outdata/mm_minimized/minimized.pdb -opdbqt -xr -O output_pdbqt outdata/mm_minimized/output.pdbqt

# 还差添加离子 生成拓扑文件 定义模拟盒子 添加水分子
