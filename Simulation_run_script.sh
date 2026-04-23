#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=48:00:00
#SBATCH --job-name=wat_tryp_gpu
#SBATCH --error=tryp_gpu_error
#SBATCH --output=tryp_0_gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1  
 
module load cmake/3.14.3
module load compiler/gcc/7.3.0
module load compiler/intel/2018.2.199
module load compiler/cuda/10.1
module load python/conda-python/3.9

source /home/parsapeer20/softwares/gromacs2018_gpu/bin/GMXRC
 
export I_MPI_FABRICS=shm:dapl
export OMP_NUM_THREADS=40
 
gmx_mpi grompp -f step6.0_minimization.mdp -c step5_input.gro -r step5_input.gro -p topol.top -o em_1.tpr -maxwarn 2
gmx_mpi mdrun -deffnm em_1 -v -nb gpu -pme gpu
wait
gmx_mpi grompp -f step6.0_minimization_1.mdp -c em_1.gro -r em_1.gro -p topol.top -o em_2.tpr -maxwarn 2
gmx_mpi mdrun -deffnm em_2 -v -nb gpu -pme gpu
wait
gmx_mpi grompp -f annealing.mdp -c em_2.gro -r em_2.gro -p topol.top -o anneal.tpr -n index.ndx
gmx_mpi mdrun -deffnm anneal -v -nb gpu -pme gpu

#### EQUILIBRIATION #######
########################### NVT-LIPIDS ###################################################################################
# NVT Equilibriation at of position restraints [10000]
gmx_mpi grompp -f nvt_lipids.mdp -c em_2.gro -r em_2.gro -p topol.top -o nvt_lipids.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_lipids -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 8000]
gmx_mpi grompp -f nvt_lipids_1.mdp -c nvt_lipids.gro -r nvt_lipids.gro -p topol.top -o nvt_lipids_1.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_lipids_1 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 6000]
gmx_mpi grompp -f nvt_lipids_2.mdp -c nvt_lipids_1.gro -r nvt_lipids_1.gro -p topol.top -o nvt_lipids_2.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_lipids_2 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 4000]
gmx_mpi grompp -f nvt_lipids_3.mdp -c nvt_lipids_2.gro -r nvt_lipids_2.gro -p topol.top -o nvt_lipids_3.tpr -n index.ndx 
gmx_mpi mdrun -deffnm nvt_lipids_3 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 2000]
gmx_mpi grompp -f nvt_lipids_4.mdp -c nvt_lipids_3.gro -r nvt_lipids_3.gro -p topol.top -o nvt_lipids_4.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_lipids_4 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 1000]
gmx_mpi grompp -f nvt_lipids_5.mdp -c nvt_lipids_4.gro -r nvt_lipids_4.gro -p topol.top -o nvt_lipids_5.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_lipids_5 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 500]
gmx_mpi grompp -f nvt_lipids_6.mdp -c nvt_lipids_5.gro -r nvt_lipids_5.gro -p topol.top -o nvt_lipids_6.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_lipids_6 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 0]
gmx_mpi grompp -f nvt_lipids_7.mdp -c nvt_lipids_6.gro -r nvt_lipids_6.gro -p topol.top -o nvt_lipids_7.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_lipids_7 -v -nb gpu -pme gpu
wait 

# ########################### NVT-SIDECHAIN ###################################################################################

# # NVT Equilibriation at of position restraints [10000]
gmx_mpi grompp -f nvt_sc.mdp -c nvt_lipids_7.gro -r nvt_lipids_7.gro -p topol.top -o nvt_sc.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_sc -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 8000]
gmx_mpi grompp -f nvt_sc_1.mdp -c nvt_sc.gro -r nvt_sc.gro -p topol.top -o nvt_sc_1.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_sc_1 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 6000]
gmx_mpi grompp -f nvt_sc_2.mdp -c nvt_sc_1.gro -r nvt_sc_1.gro -p topol.top -o nvt_sc_2.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_sc_2 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 4000]
gmx_mpi grompp -f nvt_sc_3.mdp -c nvt_sc_2.gro -r nvt_sc_2.gro -p topol.top -o nvt_sc_3.tpr -n index.ndx 
gmx_mpi mdrun -deffnm nvt_sc_3 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 2000]
gmx_mpi grompp -f nvt_sc_4.mdp -c nvt_sc_3.gro -r nvt_sc_3.gro -p topol.top -o nvt_sc_4.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_sc_4 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 1000]
gmx_mpi grompp -f nvt_sc_5.mdp -c nvt_sc_4.gro -r nvt_sc_4.gro -p topol.top -o nvt_sc_5.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_sc_5 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 500]
gmx_mpi grompp -f nvt_sc_6.mdp -c nvt_sc_5.gro -r nvt_sc_5.gro -p topol.top -o nvt_sc_6.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_sc_6 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 0]
gmx_mpi grompp -f nvt_sc_7.mdp -c nvt_sc_6.gro -r nvt_sc_6.gro -p topol.top -o nvt_sc_7.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_sc_7 -v -nb gpu -pme gpu
wait

############################################ NVT-BACKBONE #############################################################
# NVT Equilibriation at of position restraints [10000]
gmx_mpi grompp -f nvt_bb.mdp -c nvt_sc_7.gro -r nvt_sc_7.gro -p topol.top -o nvt_bb.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_bb -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 8000]
gmx_mpi grompp -f nvt_bb_1.mdp -c nvt_bb.gro -r nvt_bb.gro -p topol.top -o nvt_bb_1.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_bb_1 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 6000]
gmx_mpi grompp -f nvt_bb_2.mdp -c nvt_bb_1.gro -r nvt_bb_1.gro -p topol.top -o nvt_bb_2.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_bb_2 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 4000]
gmx_mpi grompp -f nvt_bb_3.mdp -c nvt_bb_2.gro -r nvt_bb_2.gro -p topol.top -o nvt_bb_3.tpr -n index.ndx 
gmx_mpi mdrun -deffnm nvt_bb_3 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 2000]
gmx_mpi grompp -f nvt_bb_4.mdp -c nvt_bb_3.gro -r nvt_bb_3.gro -p topol.top -o nvt_bb_4.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_bb_4 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 1000]
gmx_mpi grompp -f nvt_bb_5.mdp -c nvt_bb_4.gro -r nvt_bb_4.gro -p topol.top -o nvt_bb_5.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_bb_5 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 500]
gmx_mpi grompp -f nvt_bb_6.mdp -c nvt_bb_5.gro -r nvt_bb_5.gro -p topol.top -o nvt_bb_6.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_bb_6 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 0]
gmx_mpi grompp -f nvt_bb_7.mdp -c nvt_bb_6.gro -r nvt_bb_6.gro -p topol.top -o nvt_bb_7.tpr -n index.ndx
gmx_mpi mdrun -deffnm nvt_bb_7 -v -nb gpu -pme gpu
wait

############################################ NPT-LIPIDS #############################################################

# NVT Equilibriation at of position restraints [10000]
gmx_mpi grompp -f npt_lipid.mdp -c nvt_bb_7.gro -r nvt_bb_7.gro -p topol.top -o npt_lipid.tpr -n index.ndx -maxwarn 2
gmx_mpi mdrun -deffnm npt_lipid -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 8000]
gmx_mpi grompp -f npt_lipid_1.mdp -c npt_lipid.gro -r npt_lipid.gro -p topol.top -o npt_lipid_1.tpr -n index.ndx -maxwarn 2
gmx_mpi mdrun -deffnm npt_lipid_1 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 6000]
gmx_mpi grompp -f npt_lipid_2.mdp -c npt_lipid_1.gro -r npt_lipid_1.gro -p topol.top -o npt_lipid_2.tpr -n index.ndx -maxwarn 2
gmx_mpi mdrun -deffnm npt_lipid_2 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 4000]
gmx_mpi grompp -f npt_lipid_3.mdp -c npt_lipid_2.gro -r npt_lipid_2.gro -p topol.top -o npt_lipid_3.tpr -n index.ndx -maxwarn 2
gmx_mpi mdrun -deffnm npt_lipid_3 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 2000]
gmx_mpi grompp -f npt_lipid_4.mdp -c npt_lipid_3.gro -r npt_lipid_3.gro -p topol.top -o npt_lipid_4.tpr -n index.ndx -maxwarn 2
gmx_mpi mdrun -deffnm npt_lipid_4 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 1000]
gmx_mpi grompp -f npt_lipid_5.mdp -c npt_lipid_4.gro -r npt_lipid_4.gro -p topol.top -o npt_lipid_5.tpr -n index.ndx -maxwarn 2
gmx_mpi mdrun -deffnm npt_lipid_5 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 500]
gmx_mpi grompp -f npt_lipid_6.mdp -c npt_lipid_5.gro -r npt_lipid_5.gro -p topol.top -o npt_lipid_6.tpr -n index.ndx -maxwarn 2
gmx_mpi mdrun -deffnm npt_lipid_6 -v -nb gpu -pme gpu
wait 
# NVT equilibriation [ 0]
gmx_mpi grompp -f npt_lipid_7.mdp -c npt_lipid_6.gro -r npt_lipid_6.gro -p topol.top -o npt_lipid_7.tpr -n index.ndx -maxwarn 2
gmx_mpi mdrun -deffnm npt_lipid_7 -v -nb gpu -pme gpu
wait

# ################################################ NPT-SIDECHAINS ###########################################################

# NPT equilibriation at of position restraints [10000] [side chain]
gmx_mpi grompp -f npt_sc.mdp -c npt_lipid_7.gro -r npt_lipid_7.gro -p topol.top -o npt_sc.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_sc -v -nb gpu -pme gpu
wait 
# NPT equilibriation [ 8000]
gmx_mpi grompp -f npt_sc_1.mdp -c npt_sc.gro -r npt_sc.gro -p topol.top -o npt_sc_1.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_sc_1 -v -nb gpu -pme gpu
wait 
# NPT equilibriation [ 6000]
gmx_mpi grompp -f npt_sc_2.mdp -c npt_sc_1.gro -r npt_sc_1.gro -p topol.top -o npt_sc_2.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_sc_2 -v -nb gpu -pme gpu
wait 
# NPT equilibriation [ 4000]
gmx_mpi grompp -f npt_sc_3.mdp -c npt_sc_2.gro -r npt_sc_2.gro -p topol.top -o npt_sc_3.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_sc_3 -v -nb gpu -pme gpu
wait 
# NPT equilibriation [ 2000]
gmx_mpi grompp -f npt_sc_4.mdp -c npt_sc_3.gro -r npt_sc_3.gro -p topol.top -o npt_sc_4.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_sc_4 -v -nb gpu -pme gpu
wait 
# NPT equilibriation [ 1000]
gmx_mpi grompp -f npt_sc_5.mdp -c npt_sc_4.gro -r npt_sc_4.gro -p topol.top -o npt_sc_5.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_sc_5 -v -nb gpu -pme gpu
wait 
# NPT equilibriation [ 500]
gmx_mpi grompp -f npt_sc_6.mdp -c npt_sc_5.gro -r npt_sc_5.gro -p topol.top -o npt_sc_6.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_sc_6 -v -nb gpu -pme gpu
wait 
# NPT equilibriation [ 0]
gmx_mpi grompp -f npt_sc_7.mdp -c npt_sc_6.gro -r npt_sc_6.gro -p topol.top -o npt_sc_7.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_sc_7 -v -nb gpu -pme gpu
wait 

################################################ NPT-BACKBONE ##############################################

#NPT equilibriation at of position restraints [10000] [backbone]
gmx_mpi grompp -f npt_bb.mdp -c npt_sc_7.gro -r npt_sc_7.gro -p topol.top -o npt_bb.tpr -n index.ndx -maxwarn 2
gmx_mpi mdrun -deffnm npt_bb -v -nb gpu -pme gpu
wait 
#NPT equilibriation [ 8000]
gmx_mpi grompp -f npt_bb_1.mdp -c npt_bb.gro -r npt_bb.gro -p topol.top -o npt_bb_1.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_bb_1 -v -nb gpu -pme gpu
wait 
#NPT equilibriation [ 6000]
gmx_mpi grompp -f npt_bb_2.mdp -c npt_bb_1.gro -r npt_bb_1.gro -p topol.top -o npt_bb_2.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_bb_2 -v -nb gpu -pme gpu
wait 
#NPT equilibriation [ 4000]
gmx_mpi grompp -f npt_bb_3.mdp -c npt_bb_2.gro -r npt_bb_2.gro -p topol.top -o npt_bb_3.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_bb_3 -v -nb gpu -pme gpu
wait 
#NPT equilibriation [ 2000]
gmx_mpi grompp -f npt_bb_4.mdp -c npt_bb_3.gro -r npt_bb_3.gro -p topol.top -o npt_bb_4.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_bb_4 -v -nb gpu -pme gpu
wait 
#NPT equilibriation [ 1000]
gmx_mpi grompp -f npt_bb_5.mdp -c npt_bb_4.gro -r npt_bb_4.gro -p topol.top -o npt_bb_5.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_bb_5 -v -nb gpu -pme gpu
wait 
#NPT equilibriation [ 500]
gmx_mpi grompp -f npt_bb_6.mdp -c npt_bb_5.gro -r npt_bb_5.gro -p topol.top -o npt_bb_6.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_bb_6 -v -nb gpu -pme gpu
wait 
#NPT equilibriation [ 0]
gmx_mpi grompp -f npt_bb_7.mdp -c npt_bb_6.gro -r npt_bb_6.gro -p topol.top -o npt_bb_7.tpr -n index.ndx
gmx_mpi mdrun -deffnm npt_bb_7 -v -nb gpu -pme gpu
wait

# RELAXATION 200ns
gmx_mpi grompp -f short_run.mdp -c npt_bb_7.gro -r npt_bb_7.gro -p topol.top -o unconstrain_200.tpr -n index.ndx
gmx_mpi mdrun -deffnm unconstrain_200 -v -nb gpu -pme gpu

#gmx_mpi grompp -f step7_production.mdp -c short_run.gro -p topol.top -n index.ndx -o mdp_1.tpr
##-machinefile $PBS_NODEFILE -n 40 gmx_mpi mdrun -deffnm mdp_1
gmx_mpi grompp -f step7_production.mdp -c short_run.gro -p topol.top -n index.ndx -o mdp_1.tpr
gmx_mpi mdrun -deffnm mdp_1
