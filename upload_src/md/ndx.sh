#!/bin/bash

# 运行 gmx make_ndx 命令并传递输入
gmx make_ndx -f pr.gro <<EOF
1 | 13
!20
name 20 protein_lig
name 21 envir
q
EOF

echo "gmx make_ndx 命令已完成。"

