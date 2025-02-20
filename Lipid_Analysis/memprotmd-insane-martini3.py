import os
import sys
import gromacs
import MDAnalysis
import glob
import pandas as pd
import shutil
import subprocess
from gromacs import genion, grompp, editconf, make_ndx, trjconv, confrms, pdb2gmx, mdrun, mdrun_d

def run_insane3(lipid, x, y, z, dm):
    lipid_args = lipid.split()
    cmd = [
        'python3.11', '/storage/chem/lfsmgr/SRG/MemProtMD3/insane3.py',
        *lipid_args, '-salt', '0.15', '-sol', 'W', '-o', 'CG-system.gro',
        '-p', 'topol.top', '-f', 'protein-em.pdb', '-center',
        '-x', f"{x}", '-y', f"{y}", '-z', f"{z}", '-dm', f"{dm}",
    ]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.decode('utf-8'))
    if stderr:
        print(stderr.decode('utf-8'))
        
def update_topology_file(topol, martini):
    replacements = {
        'NA+': 'NA',
        'CL-': 'CL',
        '#include \"'+martini+'.itp\"': (
            '#include "martini_v3.0.0.itp"\n'
            '#include "martini_v3.0.0_ffbonded_v2_openbeta.itp"\n'
            '#include "martini_v3.0.0_ions_v1.itp"\n'
            '#include "martini_v3.0.0_solvents_v1.itp"\n'
            '#include "martini_v3.0.0_phospholipids_v1.itp"\n'
        ),
    }
    lines = []
    with open(topol) as infile:
        for line in infile:
            for src, target in replacements.items():
                line = line.replace(src, target)
            lines.append(line)
    with open(topol, 'w') as outfile:
        for line in lines:
            outfile.write(line)

if len(sys.argv) == 3:
    filename = sys.argv[1]
    lipid = sys.argv[2]
    name = os.path.splitext(filename)[0]
    repeat = 1
    ss = 'alpha'
elif len(sys.argv) == 4:
    filename = sys.argv[1]
    lipid = sys.argv[2]
    name = os.path.splitext(filename)[0]
    repeat = int(sys.argv[3])
    ss = 'alpha'
elif len(sys.argv) == 5:
    filename = sys.argv[1]
    lipid = sys.argv[2]
    name = os.path.splitext(filename)[0]
    repeat = int(sys.argv[3])
    ss = sys.argv[4] 
else:
    print("\nDon't forget an input PDB file and details of the lipids to be included: \n\npython memprotmd-insane-local.py protein.pdb '-l lipids' num-repeats\n\ne.g.\n\npython memprotmd-insane-local.py protein.pdb '-l POPE:4 -l POPG:1 -l CARD:1' 1\n\n")
    sys.exit(1)


if not os.path.exists('martini_v3.0.0.itp'):
	os.system("wget https://github.com/pstansfeld/MemProtMD/raw/main/martini_v300.zip")
	os.system("unzip -o martini_v300.zip")
	os.system("mv martini_v300/* .")
	os.system("cp /storage/chem/lfsmgr/SRG/MemProtMD3/martini300/*phospholipids* .")
	os.system("cp /storage/chem/lfsmgr/SRG/MemProtMD3/martini300/*ffbonded* .")
if not os.path.exists('insane3.py'):
	os.system("wget https://github.com/pstansfeld/MemProtMD/raw/main/insane3.py")

if ss == 'alpha':
	os.system('/storage/chem/lfsmgr/SRG/MemProtMD3/memembed/bin/memembed -o memembed.pdb ' + filename)
elif ss == 'beta':
	os.system('/storage/chem/lfsmgr/SRG/MemProtMD3/memembed/bin/memembed -b -o memembed.pdb ' + filename)
		

make_ndx(f='memembed.pdb', o='index.ndx', input=('del 0', 'del 1-100','rDUM','q'),backup=False)
editconf(f='memembed.pdb',o='centered.pdb',n='index.ndx',c=True,d=3,input=(0,0),backup=False)

u = MDAnalysis.Universe('centered.pdb')

confrms(f2='memembed.pdb', f1='centered.pdb', one=True, o='aligned.gro', input=(3,3),backup=False)
editconf(f='aligned.gro', o='protein.pdb',label='A',resnr=1,n=True,input=(0,0),backup=False)

v = MDAnalysis.Universe('aligned.gro')
dum = v.select_atoms('resname DUM')
box = (u.dimensions/10.0)

dm = (box[2]/2) - ((dum.center_of_mass()[2]/10))

editconf(f='memembed.pdb',o='centered.pdb',label='A',resnr=1,n='index.ndx',c=True,input=(1,0),backup=False,box=[round(box[0]), round(box[1]), round(box[2])])

with open('centered.pdb', 'r') as file :
  filedata = file.read()
filedata = filedata.replace('HSE', 'HIS')
filedata = filedata.replace('HSD', 'HIS')
filedata = filedata.replace('MSE', 'MET')
filedata = filedata.replace(' SE ', ' SD ')
with open('centered.pdb', 'w') as file:
  file.write(filedata)

with open('em.mdp','w') as em:
            em.write('define = -DFLEXIBLE\nintegrator = steep\nnsteps = 25000\nemtol = 100\nemstep = 0.001')

os.system("martinize2 -f centered.pdb  -ff martini3001 -x protein-cg.pdb -o protein-cg.top -dssp /storage/chem/lfsmgr/SRG/miniconda3/bin/mkdssp -scfix -elastic -ef 500 -eu 0.9 -el 0.5 -ea 0 -ep 0 -merge A -maxwarn 1000")
os.system("sed -e 's/^molecule.*/Protein 1/g' molecule*.itp >  protein-cg.itp")

update_topology_file('protein-cg.top', "martini")

editconf(f='protein-cg.pdb',o='protein-cg.pdb',d=3,backup=False)

grompp(f='em.mdp',o='em.tpr',c='protein-cg.pdb',p='protein-cg.top',maxwarn='100',backup=False,v=True)
mdrun_d(deffnm='em', c='protein-em.pdb',backup=False, ntmpi=1)
trjconv(f='protein-em.pdb', o='protein-em.pdb', pbc='res', s='em.tpr', conect=True, input='0',backup=False)

run_insane3(lipid, round(box[0]), round(box[1]), round(box[2]), round(dm))      

update_topology_file('topol.top', "martini_v3")

with open('equil1.mdp','w') as md:
            md.write('integrator = md\ntinit = 0.0\ndt = 0.02\nnsteps = 500000\nnstxout = 0\nnstvout = 0\nnstfout = 0\nnstlog = 50000\nnstenergy = 50000\nnstxout-compressed = 50000\ncompressed-x-precision = 10000\nnstlist  = 10\nns_type  = grid\npbc   = xyz\ncoulombtype  = Reaction_field\nrcoulomb_switch = 0.0\nrcoulomb  = 1.1\nepsilon_r  = 15\nvdw_type  = cutoff\nrvdw_switch  = 0.9\nrvdw   = 1.1\ncutoff-scheme = verlet\ncoulomb-modifier = Potential-shift\nvdw-modifier  = Potential-shift\nepsilon_rf  = 0\nverlet-buffer-tolerance = 0.005\ntcoupl  = v-rescale\ntc-grps  = PROTEIN LIPID SOL_ION\ntau_t  = 1.0 1.0 1.0\nref_t  = 310 310 310\nPcoupl  = no\ngen_vel  = yes\ngen_temp  = 310\ngen_seed  = -1\nconstraints  = none\nconstraint_algorithm = Lincs\ncontinuation  = no\nlincs_order  = 4\nlincs_warnangle = 30\n')
with open('equil2.mdp','w') as md:
            md.write('integrator = md\ntinit = 0.0\ndt = 0.02\nnsteps = 500000\nnstxout = 0\nnstvout = 0\nnstfout = 0\nnstlog = 50000\nnstenergy = 50000\nnstxout-compressed = 50000\ncompressed-x-precision = 10000\nnstlist  = 10\nns_type  = grid\npbc   = xyz\ncoulombtype  = Reaction_field\nrcoulomb_switch = 0.0\nrcoulomb  = 1.1\nepsilon_r  = 15\nvdw_type  = cutoff\nrvdw_switch  = 0.9\nrvdw   = 1.1\ncutoff-scheme = verlet\ncoulomb-modifier = Potential-shift\nvdw-modifier  = Potential-shift\nepsilon_rf  = 0\nverlet-buffer-tolerance = 0.005\ntcoupl  = v-rescale\ntc-grps  = PROTEIN LIPID SOL_ION\ntau_t  = 1.0 1.0 1.0\nref_t  = 310 310 310\nPcoupl  = berendsen\nPcoupltype  = semiisotropic\ntau_p  = 4.0\ncompressibility = 3e-4 3e-4\nref_p  = 1.0 1.0\ngen_vel  = no\nconstraints  = none\nconstraint_algorithm = Lincs\ncontinuation  = yes\nlincs_order  = 4\nlincs_warnangle = 30\n')
with open('cgmd.mdp','w') as md:
            md.write('integrator = md\ntinit = 0.0\ndt = 0.02\nnsteps = 50000000\nnstxout = 0\nnstvout = 0\nnstfout = 0\nnstlog = 50000\nnstenergy = 50000\nnstxout-compressed = 50000\ncompressed-x-precision = 10000\nnstlist  = 10\nns_type  = grid\npbc   = xyz\ncoulombtype  = Reaction_field\nrcoulomb_switch = 0.0\nrcoulomb  = 1.1\nepsilon_r  = 15\nvdw_type  = cutoff\nrvdw_switch  = 0.9\nrvdw   = 1.1\ncutoff-scheme = verlet\ncoulomb-modifier = Potential-shift\nvdw-modifier  = Potential-shift\nepsilon_rf  = 0\nverlet-buffer-tolerance = 0.005\ntcoupl  = v-rescale\ntc-grps  = PROTEIN LIPID SOL_ION\ntau_t  = 1.0 1.0 1.0\nref_t  = 310 310 310\nPcoupl  = c-rescale\nPcoupltype  = semiisotropic\ntau_p  = 4.0\ncompressibility = 3e-4 3e-4\nref_p  = 1.0 1.0\ngen_vel  = yes\ngen_temp  = 310\ngen_seed  = -1\nconstraints  = none\nconstraint_algorithm = Lincs\ncontinuation  = no\nlincs_order  = 4\nlincs_warnangle = 30\n')

grompp(f='em.mdp',o='em.tpr',c='CG-system.gro',maxwarn='100',backup=False,v=True)

genion(s='em.tpr',neutral=True,o='CG-system.gro',p='topol.top',input='W')
grompp(f='em.mdp',o='em.tpr',c='CG-system.gro',maxwarn='100',backup=False,v=True)
mdrun_d(deffnm='em', c='CG-system.pdb',backup=False, ntmpi=1)
trjconv(f='CG-system.pdb', o='CG-system.pdb', pbc='res', s='em.tpr', conect=True, input='0',backup=False)
make_ndx(f='CG-system.pdb', o='index.ndx', input=('del 0', 'del 1-40', '0|rPOP*','1&!0','!1','del 1','name 1 Lipid','name 2 SOL_ION','q'),backup=False)

for rep in range(1, repeat + 1):
	os.system("mkdir MD"+str(rep))
	grompp(f='equil1.mdp',o='equil1-'+str(rep),c='CG-system.pdb',maxwarn=100, n='index.ndx',backup=False)  
	mdrun(deffnm='equil1-'+str(rep),backup=False,v=True,resethway=True,nstlist=100,ntmpi=1)
	grompp(f='equil2.mdp',o='equil2-'+str(rep),c='equil1-'+str(rep)+'.gro',maxwarn=100, n='index.ndx',backup=False)  
	mdrun(deffnm='equil2-'+str(rep),backup=False,v=True,resethway=True,nstlist=100,ntmpi=1)
	grompp(f='cgmd.mdp',o='MD'+str(rep)+'/md',c='equil2-'+str(rep)+'.gro',maxwarn=100, n='index.ndx',backup=False)  
	for file in glob.glob(r'#*'):
		os.remove(file)
	for file in glob.glob(r'#*'):
		os.remove(file)