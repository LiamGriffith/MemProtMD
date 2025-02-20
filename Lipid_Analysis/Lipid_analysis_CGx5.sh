mkdir CAT

for ii in {1..5} # Modify for how many repeats you want
do 
gmx trjconv -f MD$ii-output/MD$ii/md.gro -s MD$ii-output/MD$ii/md.tpr -center -pbc res -conect -o MD$ii-output/MD$ii/md-center.pdb 
done

gmx trjcat -f MD*-output/MD*/md.xtc -o CAT/md.xtc -settime 
gmx trjconv -f CAT/md.xtc -o CAT/md-center.xtc -s MD1-output/MD1/md.tpr -center -pbc res
gmx trjconv -f CAT/md-center.xtc -o CAT/md-fit.xtc -s MD1-output/MD1/md.tpr -center -fit progressive
gmx trjconv -f CAT/md-center.xtc -o CAT/md-rot+trans.xtc -s MD1-output/MD1/md.tpr -center -fit rotxy+transxy
conda activate memprotmd
python platelet-lipid-contacts.py MD1-output/MD1/md-center.pdb CAT/md-center.xtc 
python /storage/chem/lfsmgr/SRG/MemProtMD3/lipid-deformations.py MD1-output/MD1/md-center.pdb CAT/md-rot+trans.xtc 
python platelet-densities.py MD1-output/MD1/md-center.pdb CAT/md-fit.xtc 
gmx trjconv -f MD5-output/MD5/md.gro -o system5.pdb -conect -center -s MD5-output/MD5/md.tpr -pbc res # Modify for how many repeats you want
gmx trjconv -f MD4-output/MD4/md.gro -o system4.pdb -conect -center -s MD5-output/MD5/md.tpr -pbc res # Modify for how many repeats you want
gmx trjconv -f MD3-output/MD3/md.gro -o system3.pdb -conect -center -s MD5-output/MD5/md.tpr -pbc res # Modify for how many repeats you want
gmx trjconv -f MD1-output/MD1/md.gro -o system1.pdb -conect -center -s MD5-output/MD5/md.tpr -pbc res # Modify for how many repeats you want
gmx trjconv -f MD2-output/MD2/md.gro -o system2.pdb -conect -center -s MD5-output/MD5/md.tpr -pbc res # Modify for how many repeats you want
pymol lipid-analysis.pml 

# User inputs:
# If Sum of least squares: 1
# If Centering: 1
# If Output: 0 