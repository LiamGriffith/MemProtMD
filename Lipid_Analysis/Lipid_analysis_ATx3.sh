mkdir CATAT

for ii in {1..3}
do 
gmx trjconv -f CG2AT/FINAL/PR$ii/pr.gro -s CG2AT/FINAL/PR$ii/pr.tpr -center -pbc res -conect -o CG2AT/FINAL/PR$ii/pr-center.pdb -dt 500
done

gmx trjcat -f CG2AT/FINAL/PR*/pr.xtc -o CATAT/pr.xtc -settime 
gmx trjconv -f CATAT/pr.xtc -o CATAT/pr-center.xtc -s CG2AT/FINAL/PR1/pr.tpr -center -pbc res
gmx trjconv -f CATAT/pr-center.xtc -o CATAT/pr-fit.xtc -s CG2AT/FINAL/PR1/pr.tpr -center -fit progressive
gmx trjconv -f CATAT/pr-center.xtc -o CATAT/pr-rot+trans.xtc -s CG2AT/FINAL/PR1/pr.tpr -center -fit rotxy+transxy

conda activate memprotmd

python platelet-lipid-contacts.py CG2AT/FINAL/PR1/pr-center.pdb CATAT/pr-center.xtc
python /storage/chem/lfsmgr/SRG/MemProtMD3/lipid-deformations.py CG2AT/FINAL/PR1/pr-center.pdb CATAT/pr-rot+trans.xtc
python platelet-densities.py CG2AT/FINAL/PR1/pr-center.pdb CATAT/pr-fit.xtc
gmx trjconv -f CG2AT/FINAL/PR3/pr.xtc -o system3.pdb -conect -center -s CG2AT/FINAL/PR3/pr.tpr -pbc res
gmx trjconv -f CG2AT/FINAL/PR2/pr.xtc -o system2.pdb -conect -center -s CG2AT/FINAL/PR2/pr.tpr -pbc res
gmx trjconv -f CG2AT/FINAL/PR1/pr.xtc -o system1.pdb -conect -center -s CG2AT/FINAL/PR1/pr.tpr -pbc res

pymol lipid-analysis.pml