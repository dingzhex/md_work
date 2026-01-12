#!/bin/bash

#source split.sh

source convert_ligand.sh

export PGI_FASTMATH_CPU=sandybridge

# 记录当前目录
current_dir=$(pwd)

# 移动 ligand.mol2 和 ligand.chg 文件到 home 目录
mv ligand.mol2 "/home/ding/sobtop_1.0(dev4)"

# 进入 home 目录
cd "/home/ding/sobtop_1.0(dev4)"

# 运行 ligand.sh 脚本
./ligand_nochg.sh

# 移动文件回原目录
mv ligand.itp ligand.top ligand.gro ligand.mol2 "$current_dir"

cd "$current_dir"

# 提示操作完成
echo "Files moved back to the original directory: $current_dir"

gmx pdb2gmx -f protein.pdb -o protein.gro -p topol.top  -ignh -ff amber99sb -water tip3p 

#复制gro文件到protein组成复合物，以及原子数目（手动）
source ./complex.sh

#选择组system 1000kj/mol/nm2
echo system | gmx genrestr -f ligand.gro -o posre_ligand.itp

#下列语句插入到ligand.itp文件末尾,配体的itp文件引入整体的拓扑结构top,分子数目
source ./topol.sh

gmx editconf -f complex.gro -o complex~box.gro -d 1.0 -bt cubic 

gmx solvate -cp complex~box.gro -o complex~SOL.gro -p topol.top

gmx grompp -f em.mdp -c complex~SOL.gro -p topol.top -o em.tpr -maxwarn 99
echo SOL | gmx genion -s em.tpr -p topol.top -o system.gro -neutral -conc 0.15 -pname K -nname CL

gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr -maxwarn 99
gmx mdrun -v -deffnm em -pin on -ntmpi 1 -ntomp 4

