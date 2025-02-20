import sys
import os
import MDAnalysis
from MDAnalysis.analysis.density import DensityAnalysis

# Ensure correct arguments are provided
if len(sys.argv) < 3:
    print("Usage: python platelet-densities.py <topology file> <trajectory file>")
    sys.exit(1)

# Load trajectory and topology files
try:
    u = MDAnalysis.Universe(sys.argv[1], sys.argv[2])
except Exception as e:
    print(f"Error loading files: {e}")
    sys.exit(1)

# Check if trajectory has frames
if len(u.trajectory) == 0:
    print("Error: Trajectory is empty. Check the input file:", sys.argv[2])
    sys.exit(1)

# Select the protein
Protein = u.select_atoms("protein")
if len(Protein) == 0:
    print("Error: No protein atoms found in the system.")
    sys.exit(1)

# Define lipid selections and analyze density
lipids = {
    "POPC": 40,
    "POPE": 50,
    "SSM": 15,
    "DPG3": 10,
    "CHOL": 50,
    "POPS": 15,
    "POP6": 10,
}

for lipid, count in lipids.items():
    lipid_atoms = u.select_atoms(f"resname {lipid}")
    if len(lipid_atoms) == 0:
        print(f"No atoms found for lipid: {lipid}, skipping.")
        continue
    print(f"Analyzing lipid: {lipid} ({count} molecules expected)")
    try:
        D = DensityAnalysis(lipid_atoms, delta=1.0)
        D.run()
        D.density.export(f"{lipid}.dx", type="double")
    except Exception as e:
        print(f"Error analyzing lipid {lipid}: {e}")

# Analyze water
W = u.select_atoms("resname W")
if len(W) > 0:
    try:
        D = DensityAnalysis(W, delta=1.0)
        D.run()
        D.density.export("water.dx", type="double")
    except Exception as e:
        print(f"Error analyzing water density: {e}")
else:
    print("No water atoms found in the system.")

# Save the protein structure
try:
    u.trajectory[-1]
    Protein.write("density-protein.pdb")
    print("density-protein.pdb successfully written.")
except Exception as e:
    print(f"Error saving protein structure: {e}")
    sys.exit(1)