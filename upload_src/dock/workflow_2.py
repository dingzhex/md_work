# 全局盲对接：不指定结合位点，全蛋白扫描

import os
import subprocess
import glob

# from p_7 import aa_path

aa_path = 'outdata/aa_site_randsearch_'
os.makedirs(os.path.join(aa_path), exist_ok=True)
ligand_path = 'outdata/comp_3D'
ligand_list = glob.glob(os.path.join(ligand_path, 'ligand*'))


for ligand in ligand_list:
     os.makedirs(os.path.join(aa_path, 'ligand'), exist_ok=True)
     command = [
          "cp", '-r',
          ligand,  # 使用文件名而不是完整路径
          os.path.join(aa_path, 'ligand')
     ]
     result = subprocess.run(
          command,
          check=True,
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE,
          text=True,
          encoding='utf-8'
     )



ligand_path = os.path.join(aa_path, 'ligand')
ligand_list = glob.glob(os.path.join(ligand_path, 'ligand*'))
# 生成配体pdbqt文件
for path in ligand_list:
    command = [
        "obabel",
        'obabel_ligand.mol2',  # 使用文件名而不是完整路径
        "-O", 'ligand.pdbqt',
    ]

    result = subprocess.run(
        command,
        cwd=path,  # 指定工作目录
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding='utf-8'
    )


ligand_list = os.listdir(ligand_path)
os.makedirs('outdata/aa_site_randsearch_/res', exist_ok=True)
for ligand in ligand_list:
     1
     command_str = (
          f"vina --receptor 'output.pdbqt' "
          f"--ligand 'ligand/{ligand}/ligand.pdbqt' "
          f"--out 'res/{ligand}_output_file.pdbqt' "
          f"--size_x 30 --size_y 30 --size_z 30 --exhaustiveness=64 --num_modes=20"
          f"2>&1 | tee -a res/{ligand}_log_file.txt"
     )


     # 使用单字符串命令 + shell=True
     result = subprocess.run(
          command_str,
          cwd='outdata/aa_site_randsearch_',
          check=True,
          shell=True,
          text=True,
          encoding='utf-8'
     )


# vina --receptor protein.pdbqt --ligand ligand.pdbqt \
#      --size_x 60 --size_y 60 --size_z 60 \
#      --exhaustiveness=64 --num_modes=20


# 3. 计算对接盒子中心（自动计算）
# fpocket -f protein.pdb
#
#
#
#
# # 计算结合能
# best_energy=$(grep "Estimated Free Energy of Binding" dock.dlg | awk '{print $8}' | head -1)
# echo "最佳结合能: $best_energy kcal/mol"

# load ligand_3d_initial_91_output_file, ligand

