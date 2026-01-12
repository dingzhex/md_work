#!/bin/bash

# 1. 将三行插入 ligand.itp 文件的末尾
cat <<EOL >> ligand.itp
#ifdef POSRES
#include "posre_ligand.itp"
#endif
EOL

# 2. 将 #include "ligand.itp" 插入 topol.top 文件中 ; Include forcefield parameters 的下面两行处
awk '
/; Include forcefield parameters/ {
    print
    getline
    print
    getline
    print
    print "#include \"ligand.itp\""
    next
}
{print}
' topol.top > temp.top && mv temp.top topol.top

# 3. 在 topol.top 文件的最后一行插入 ligand              1
echo "ligand              1" >> topol.top

echo "文件修改完成。"

