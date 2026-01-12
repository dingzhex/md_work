#!/bin/bash

# 定义输入文件和输出文件
input_file="protein_clean.pdb"
protein_output="protein.pdb"
ligand_output="ligand.pdb"

# 使用 grep 筛选并写入对应的文件
grep '^ATOM' $input_file > $protein_output
grep '^HETATM' $input_file > $ligand_output

# 提示用户操作完成
echo "ATOM lines written to $protein_output"
echo "HETATM lines written to $ligand_output"

