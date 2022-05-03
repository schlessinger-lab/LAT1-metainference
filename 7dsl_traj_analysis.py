#!/usr/bin/env python

#Import required python modules and dependencies
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib import distances


### Import trajectory and topography files
u = mda.Universe('7dsl_MetaI.gro', 'traj_all_PBC.xtc', guess_bonds=True)
print(len(u.trajectory))


### 1. RMSD analysis
from MDAnalysis.analysis import align, rms

#define variables
Protein_Backbone = 'backbone'
JX_078 = 'resname HO0' 
Binding_Site = 'backbone and (resid 62-68 or resid 140-148 or resid 252-259 or resid 400-406)'

#setup rmsd calculations 
R = rms.RMSD(u,u, select='backbone', groupselections=[Protein_Backbone, JX_078, Binding_Site], ref_frame=0)
R.run()

#create dataframe
df = pd.DataFrame(R.rmsd, columns=['Frame', 'Time (ps)', 'Backbone','Protein_Backbone','JX_078', 'Binding_Site'])

#Transform Time from ps to ns
df['Time (ns)'] = df['Time (ps)']*0.001

#plot RMSD Graphs
plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]

#backbone
ax = df.plot(x='Time (ns)', y=['Protein_Backbone'],kind='line')
ax.set_ylabel(r'RMSD ($\AA$)')
ax.axes.xaxis.label.set_size(18)
ax.axes.yaxis.label.set_size(18)
plt.xticks(np.arange(0, 35, 5))
ax.figure.savefig('./RMSD1.png')

#ligand
ax2 = df.plot(x='Time (ns)', y='JX_078', kind='line')
ax2.set_ylabel('RMSD ($\AA$)')
ax2.axes.xaxis.label.set_size(18)
ax2.axes.yaxis.label.set_size(18)
plt.xticks(np.arange(0, 35, 5))
ax2.figure.savefig('./RMSD2.png')


### 2. RMSF Analysis
#set up RMSF
average = align.AverageStructure(u, u, select='protein and name CA', ref_frame=0).run(verbose=True)
avg_ref = average.universe
aligner = align.AlignTraj(u, avg_ref, select='protein and name CA', in_memory=True).run(verbose=True)
c_alphas = u.select_atoms('protein and name CA')

#run RMSF
R2 = rms.RMSF(c_alphas).run(verbose=True)

#plot RMSF
ax3=plt.plot(c_alphas.resids, R2.rmsf)
plt.xlabel('Residue number')
plt.ylabel('RMSF ($\AA$)')
plt.axvspan(63, 68, zorder=0, alpha=0.2, color='darkorange', label='Binding Region I')
plt.axvspan(140, 148, zorder=0, alpha=0.2, color='darkgreen', label='Binding Region II')
plt.axvspan(251, 259, zorder=0, alpha=0.2, color='purple', label='Binding Region III')
plt.axvspan(400, 408, zorder=0, alpha=0.2, color='darkred', label='Binding Region IV')
plt.ylabel('RMSF ($\AA$)', fontsize=18)
plt.xlabel('Residue Number', fontsize=18)
plt.legend()
plt.savefig('./RMSF.png')


### 3. Distance/Contact Analysis
from MDAnalysis.analysis import contacts

sel_HO0 = "resname HO0"
sel_Key_Res = "resid 62 63 65 66 67 252 255 259 405"
JX_078 = u.select_atoms(sel_HO0)
Key_Res = u.select_atoms(sel_Key_Res)

def fraction_contacts_between(r, r0, radius=3.2, min_radius=2.5):
    is_in_contact = (r < radius) & (r > min_radius)  # array of bools
    fraction = is_in_contact.sum()/r.size
    return fraction

ca = contacts.Contacts(u,
                       select=(sel_HO0, sel_Key_Res),
                       refgroup=(JX_078, Key_Res),
                       method=fraction_contacts_between,
                       radius=5.0,
                       kwargs={'radius': 5.0,
                               'min_radius': 2.4}).run();

ca_df = pd.DataFrame(ca.timeseries,
                    columns=['Time (ns)', 'Contacts with Key Residues'])
ca_df.head()

#Transform Time from ps to ns
ca_df['Time (ns)'] = ca_df['Time (ns)']*0.001
ca_df.head()

#plot
ax4 = ca_df.plot(x='Time (ns)')
plt.ylabel('Fraction of Contacts')
ax4.axes.xaxis.label.set_size(16)
ax4.axes.yaxis.label.set_size(18)
plt.xticks(np.arange(0, 35, 5))
ax4.figure.savefig('./fraction_contacts_overtime.png'
