#!/bin/bash
#SBATCH --job-name=gromacs-prod1-6
#SBATCH --partition=sonmi
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -w compute-0-0
#SBATCH -n 16
#SBATCH --time=10-00:00:00

module load gmx2025_gpu

cd prod_rep1
gmx grompp -f prod_30ns.mdp -c step4.2_equilibration_npt.gro -t step4.2_equilibration_npt.cpt -p topol.top -n index.ndx -o prod_30ns_rep1.tpr -maxwarn 2
gmx mdrun -v -deffnm prod_30ns_rep1 -nb gpu -pme gpu -bonded gpu -ntomp 24 -ntmpi 4 -npme 1 -pin on -dlb yes -cpt 5 -cpo prod_30ns_rep1.cpt


cd ../prod_rep2
gmx grompp -f prod_30ns.mdp -c step4.2_equilibration_npt.gro -t step4.2_equilibration_npt.cpt -p topol.top -n index.ndx -o prod_30ns_rep2.tpr -maxwarn 2
gmx mdrun -v -deffnm prod_30ns_rep2 -nb gpu -pme gpu -bonded gpu -ntomp 24 -ntmpi 4 -npme 1 -pin on -dlb yes -cpt 5 -cpo prod_30ns_rep2.cpt


cd ../prod_rep3
gmx grompp -f prod_30ns.mdp -c step4.2_equilibration_npt.gro -t step4.2_equilibration_npt.cpt -p topol.top -n index.ndx -o prod_30ns_rep3.tpr -maxwarn 2
gmx mdrun -v -deffnm prod_30ns_rep3 -nb gpu -pme gpu -bonded gpu -ntomp 24 -ntmpi 4 -npme 1 -pin on -dlb yes -cpt 5 -cpo prod_30ns_rep3.cpt

cd ../prod_rep4
gmx grompp -f prod_30ns.mdp -c step4.2_equilibration_npt.gro -t step4.2_equilibration_npt.cpt -p topol.top -n index.ndx -o prod_30ns_rep4.tpr -maxwarn 2
gmx mdrun -v -deffnm prod_30ns_rep4 -nb gpu -pme gpu -bonded gpu -ntomp 24 -ntmpi 4 -npme 1 -pin on -dlb yes -cpt 5 -cpo prod_30ns_rep4.cpt


cd ../prod_rep5
gmx grompp -f prod_30ns.mdp -c step4.2_equilibration_npt.gro -t step4.2_equilibration_npt.cpt -p topol.top -n index.ndx -o prod_30ns_rep5.tpr -maxwarn 2
gmx mdrun -v -deffnm prod_30ns_rep5 -nb gpu -pme gpu -bonded gpu -ntomp 24 -ntmpi 4 -npme 1 -pin on -dlb yes -cpt 5 -cpo prod_30ns_rep5.cpt


cd ../prod_rep6
gmx grompp -f prod_30ns.mdp -c step4.2_equilibration_npt.gro -t step4.2_equilibration_npt.cpt -p topol.top -n index.ndx -o prod_30ns_rep6.tpr -maxwarn 2
gmx mdrun -v -deffnm prod_30ns_rep6 -nb gpu -pme gpu -bonded gpu -ntomp 24 -ntmpi 4 -npme 1 -pin on -dlb yes -cpt 5 -cpo prod_30ns_rep6.cpt

exit
