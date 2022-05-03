#!/bin/sh
#
# To run the gromos clustering on the pre-calculated RMSD matrices, you need first
# to compile gromos_clustering.cpp. If you have GNU C++ compiler, you can do:

g++ -O3 gromos_clustering.cpp -o gromos_clustering.x

# Load Gromacs
module load load gromacs/2020.1-plumed-isdb

# Set gromacs excecutable
GMX=gmx_s_mpi
# trajectory file
traj=traj_all_PBC.xtc

# count number of frames
nframes=`${GMX} check -f ${traj} |& grep Step | awk '{print $2}'`

#check expected frames
echo $nframes

# rmsd matrix prefix
RMSD_MATRIX=/4-ANALYSIS/4.2-RMSD-MATRIX/FIXED/RMSD
# cutoff [Ang]
CUT=1.5

# do clustering
./gromos_clustering.x $RMSD_MATRIX 128 $CUT $nframes


# Output of gromos_clustering
# NOTE: Numbering scheme for clusters and frame IDs starts from 0!
# You will obtain two output files:
# - log.dat: information about the cluster founds: population, and cluster center.
# - trajectory.dat: cluster assignement for each frame of the original trajectory.
