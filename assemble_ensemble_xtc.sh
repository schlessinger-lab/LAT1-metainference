#!/bin/bash
#Load Gromacs-Plumed
ml load gromacs/2020.1-plumed-isdb

# gromacs stride for writing xtc file [ps]
nstx=5

# 1) discard initial 20% of trajectory of each replica
# cycle on all directories
for d in ../3-METAI/rep-*
do
     (cd "$d" && gmx_s_mpi trjcat -f traj_comp.*.xtc -cat -o tmp.xtc) 
     (cd "$d" && nframes=`gmx_s_mpi check -f tmp.xtc |& grep Step | awk '{print $2}'`)
     (cd "$d" && b=`echo ${nframes} ${nstx} | awk '{printf "%d\n",$1*$2*0.2}'`)
     (cd "$d" && gmx_s_mpi trjconv -f tmp.xtc -o traj_rep.xtc)
     (cd "$d" && rm tmp.xtc)
done

# 2) cat together the trajectories of all replicas 
gmx_s_mpi trjcat -f ../3-METAI/rep-*/traj_rep.xtc -cat -o traj_all.xtc
# 3) run driver on global trajectory to fix PBCs
plumed driver --plumed plumed_driver.dat --mf_xtc traj_all.xtc

# 4) convert traj.gro to xtc file
echo 20 | gmx_s_mpi trjconv -f traj.gro -o traj_all_PBC.xtc -s ../3-METAI/gromacs.pdb  -n index2.ndx 

# 5) final clean
rm traj_all.xtc traj.gro                                       
