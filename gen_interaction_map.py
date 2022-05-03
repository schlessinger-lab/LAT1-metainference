#!/usr/bin/env python
import MDAnalysis as mda
import prolif as plf
import sys

#Import structure and trajectory
u = mda.Universe('./7dsl_MetaI.pdb','./traj_all_PBC.xtc', guess_bonds=True)

# create selections for the ligands and protein
prot = u.atoms.select_atoms("protein")
lig = u.atoms.select_atoms("resname HO0")

from rdkit import Chem
from rdkit.Chem import Draw
lmol = plf.Molecule.from_mda(lig)

# set default  interactions
fp = plf.Fingerprint()

# run on entire trajectory *repeat if multiple ligands are present*
fp.run(u.trajectory, lig, prot)

#create dataframe(s)
df = fp.to_dataframe()

# calculate the occurence of each interaction on the trajectory
occ = df.mean()
x1 = occ.loc[occ > 0.3]
x1.to_csv('./interaction_profile.csv')

# regroup all interactions together and do the same
g = (df.groupby(level=["protein"], axis=1)
       .sum()
       .astype(bool)
       .mean())
x2 = g.loc[g > 0.3]
x2.to_csv('./interactions_overall')

# Create ligand interaction profile
from prolif.plotting.network import LigNetwork
df2 = fp.to_dataframe(return_atoms=True)

net1 = LigNetwork.from_ifp(df2, lmol, kind="frame", frame=0, threshold=.3, rotation=270)
net1.save('JX_078_interactions_2.html')

net2 = LigNetwork.from_ifp(df2, lmol, kind="aggregate", threshold=.3, rotation=270)
net2.save('JX_078_interactions.html')
