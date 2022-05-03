#!/bin/bash
################################################################################
#Activate required modules and environments
#modules
ml gromacs/2020.1-plumed-isdb
ml anaconda3/2021.5

conda activate /conda/envs/MDA
###############################################################################
#Run Script

# trajectory file
traj=./traj_all_PBC.xtc

# count number of frames
nframes=`gmx_s_mpi check -f ${traj} |& grep Step | awk '{print $2}'`
# total number of entries in rmsd matrix
nentries=`echo ${nframes} | awk '{printf "%d\n",$1*($1-1)/2}'`
# number of entries per task - rounded up
ntask=`echo ${nentries} | awk '{printf "%d\n",$1/128+1}'`    
# id of the job
id=$((LSB_JOBINDEX - 1))
# range of entries for rmsd calculation
i0=$(( $id * $ntask ))
i1=$(( $i0 + $ntask ))
# check if over limit, and reset
if [ $i1 -gt $nentries ]; then i1=$nentries; fi

# go
python kkh_get_rmsd_matrix.py ./gromacs.pdb $traj 5.0 RMSD.$id $i0 $i1
