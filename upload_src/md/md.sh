#ixnazhi
gmx grompp -f pr.mdp -c em.gro -p topol.top -r em.gro -o pr.tpr -maxwarn 99
gmx mdrun -v -deffnm pr -pin on -ntmpi 1 -ntomp 4 

#changgui  
source ./ndx.sh

gmx grompp -f md.mdp -c pr.gro -p topol.top -o md.tpr -n index.ndx -maxwarn 99
gmx mdrun -v -deffnm md -pin on -ntmpi 1 -ntomp 4 

# gpu
# gmx mdrun -v -deffnm md -pin on -ntmpi 1 -ntomp 4 \
#          -nb gpu -pme gpu -bonded gpu -update gpu \
#          -nstlist 200