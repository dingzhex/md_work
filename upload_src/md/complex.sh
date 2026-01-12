#!/bin/bash

# 读取 ligand.gro 文件的内容并提取第三行到倒数第二行的数据
ligand_data=$(sed -n '3,$p' ligand.gro | sed '$d')

# 读取 protein.gro 文件的内容
protein_header=$(sed -n '1p' protein.gro)
protein_data=$(sed -n '3,$p' protein.gro | sed '$d')
protein_footer=$(tail -n 1 protein.gro)

# 获取 ligand.gro 和 protein.gro 的原子数目
ligand_atom_count=$(sed -n '2p' ligand.gro)
protein_atom_count=$(sed -n '2p' protein.gro)
total_atom_count=$((ligand_atom_count + protein_atom_count))

# 创建 complex.gro 文件并写入内容
{
  # 写入第一行（标题）
  echo "$protein_header"
  
  # 写入总原子数
  echo "$total_atom_count"

  # 写入 protein.gro 的数据
  echo "$protein_data"

  # 写入 ligand.gro 的数据
  echo "$ligand_data"

  # 写入 protein.gro 的最后一行
  echo "$protein_footer"
} > complex.gro

echo "complex.gro 文件已生成。"

