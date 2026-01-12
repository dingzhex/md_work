#!/bin/bash

# 执行 sobtop 命令并传递参数
./sobtop <<EOF
./ligand.mol2
2

7
10
./ligand.chg
0
1
2
4


0
EOF

