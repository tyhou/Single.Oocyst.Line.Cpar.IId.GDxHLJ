#!/bin/bash
#SBATCH --job-name=gromacs
#SBATCH --partition=sonmi
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -w compute-0-0
#SBATCH -n 4
##SBATCH --time=24:00:00

module load gmx2025_gpu

gmx grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx
gmx mdrun -v -deffnm step4.0_minimization -gpu_id 0

gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx
gmx mdrun -v -deffnm step4.1_equilibration -nb gpu -pme gpu -bonded gpu -ntomp 4 -pin on -dlb yes

gmx grompp -f step4.2_equilibration_npt.mdp -o step4.2_equilibration_npt.tpr -c step4.1_equilibration.gro -p topol.top -n index.ndx 
gmx mdrun -v -deffnm step4.2_equilibration_npt -nb gpu -pme gpu -bonded gpu -ntomp 4 -pin on -dlb yes

exit
