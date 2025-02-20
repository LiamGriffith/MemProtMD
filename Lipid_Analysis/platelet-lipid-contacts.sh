import matplotlib as mpl
mpl.use('Agg')
import sys
sys.path.append("/storage/chem/lfsmgr/SRG/MemProtMD3/")
import MDAnalysis as mda
import numpy as np
from contacts import LipidContactGenerator
import matplotlib.pyplot as plt
import pandas as pd
from pandas import *
import MDAnalysis.analysis.rms
from matplotlib import rc, rcParams
import pylab

mpl.rc_file_defaults()

# Amino acid mapping
d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# Load the system
U = mda.Universe(sys.argv[1], sys.argv[2])

# Select atoms
A = U.select_atoms("protein and name BB")
All = U.atoms

resids = A.atoms.residues.resids
resnames = A.atoms.residues.resnames

# Create residue dataframe
df = pd.DataFrame(data=A.atoms.residues.resnames, index=A.atoms.residues.resids, columns=['Residue'])

# Generate lipid contacts
generate = LipidContactGenerator(U)

# Include the additional lipids in the ligand selection
ligand_selection = (
    "resname W or resname POPE or resname POPG or resname CARD or "
    "resname POPC or resname SSM or resname DPG3 or resname CHOL or "
    "resname POPS or resname POP6"
)
contacts = generate.build_contacts(
    protein_selection="protein",
    ligand_selection=ligand_selection,
    frameskip=1,
    cutoff=6,
    KDTree=True
)
contacts.aggregate(
    group_protein_by="resid",
    group_ligand_by="resname",
    aggregate_function=lambda x: x.max()
)
data = contacts.time_aggregate(aggregate_function=lambda x: sum(x.values()) / contacts.n_frames)

# Save contact data to a CSV file
data.to_dataframe().to_csv("Lipid-contacts.csv")

# Merge residue info with contact data
oant = data.to_dataframe()
oant = pd.concat([oant, df], axis=1, sort=False, join='inner')
oant['Resid'] = oant.index
oant['combined'] = oant['Residue'] + oant.index.astype(str)

# Add contact information for each lipid and write to separate PDB files
lipids = ["POPE", "POPG", "CARD", "W", "POPC", "SSM", "DPG3", "CHOL", "POPS", "POP6"]
for lipid in lipids:
    U.add_TopologyAttr(mda.core.topologyattrs.Tempfactors(np.zeros(len(All))))
    for i in All:
        if i.residue.resid in oant.index:
            i.tempfactor = oant.loc[i.residue.resid, [lipid]].values[0] if lipid in oant.columns else 0
    U.trajectory[0]
    All.write(f"{lipid}-contacts.pdb")