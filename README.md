# LAT1-metainference

This git contains scricpts and files relating to analysis of metainference molecular dynamics(MD) simulations performed on recently published outward-occluded LAT1 cryo-EM structure. (PDB: 7DSL).

For detailed gude on how to set-up metainference simulations, GROMACS topology and input files as well as PLUMED input files, please visit: https://www.plumed-nest.org/eggs/22/018/

The scripts perform 4 tasks:

1) Assembling and concatenate ensemble simulation

    Scripts: assemble_ensemble_xtc.sh
    
2) Calculation of the RMSD matrix

    Scripts: rmsd_matrix_generation.sh, kkh_get_rmsd_matrix.py
    
3) Clustering with gromos algorithm

    Scripts: gromos_clustering.sh, gromos_clustering.cpp
    
4) Analyzing Ensemble Trajectory

    Scripts: 7dsl_traj_analysis.py, gen_interaction_map.py
    
    
# Dependencies
   
 Gromacs 2020.1
 
 MDAnalysis 2.0
  
 ProLif 3.2
 
 Python 3.6+
 
 RDKit (2020.03+)
 
 Pandas (1.0+)
 
 NumPy
 
 SciPy
 
 tqdm
 
 * Highly recommended to install MDAnalysis and ProLif in an anaconda environment. See the install instructions for those programs.


  
