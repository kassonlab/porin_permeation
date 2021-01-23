#!/usr/bin/python

import os, sys, numpy
from subprocess import call

#### set GROMACS path ################

gromacs = "/opt/gromacs-2016.5"

######################################

'''
python script defining a workflow for calculating permeabilities through OmpF porins.

This script accompanies the following publications:
Simulation-guided engineering of antibiotics for improved bacterial uptake, Ricardo Ferreir
a and Peter Kasson, doi:10.1101/2020.10.08.330332.
Antibiotic Uptake Across Gram-Negative Outer Membranes: Better Predictions Towards Better Antibiotics, Ricardo Ferreira and Peter Kasson, doi:10.1021/acsinfecdis.9b00201.

It is intended with this script to:
	- to manipulate a ligand file (stripping waters/ions/etc);
	- add it to an existing system (above the first monomer);
	- strip overlapping solvent (bulk waters)
	- change the topol file and generate and adequate index;
	- generate the tpr file for the pulling run (thus defining a set of z coordinates for US)

### version 0.0.3 ###
	- add mdp options to python script
	- add command to run GROMACS (optionally)
	- gromacs path set within the script
	- script can now be run from main directory

### version 0.0.2 ###
In this version some modifications were added:
	- VMD is no longer required, being replaced by 'gmx select'
	- For a monoanionic/dianionic compound, it also considers the presence of additional ions and possible overlap with solvent molecules
	- 'os.system' for calling GROMACS was replaced by 'subprocess.call'
	- other minor adjustments
	
### version 0.0.1 ###
In this first version:
	ligand.gro, system.gro --> system_stripped.gro
	system_stripped.gro --> topol.top, index.ndx, run.tpr
	
NOTES:
	-> The ligand file must contain only the equilibrated molecule of interest and counter-ions
	-> Still depends on VMD to strip solvent molecules from the assembled system
'''

#### read ligand and system files ####
ligand = open('ligand.gro', 'r')
lig_lines = ligand.readlines()
system = open('system.gro', 'r')
sys_lines = system.readlines()

#### read informations from files ####
sys_box = sys_lines[-1][2:30]
lig_name = lig_lines[3][5:8]
lig_topol = str(lig_lines[0][:-10] + '.itp')
counter_name = ''
if lig_lines[-2][5:8] != lig_name:
	counter_name = lig_lines[-2][5:7]

#### recenter ligand #################
call(gromacs + "/bin/gmx -quiet editconf -f ligand.gro -o lig_centered.gro -box " + str(sys_box) + " -center 6.200 8.550 8.950", shell=True)

#### check total number of atoms #####
sys_number = int(sys_lines[1].strip(' '))
lig_number = int(lig_lines[1].strip(' '))
sys_lig_number = str(sys_number + lig_number)

#### assemble the first system #######
lig_centered = open('lig_centered.gro', 'r')
lcent_lines = lig_centered.readlines()
assembly = open('system_assembled.gro' , 'w')
assembly.write('Assembled system containing OmpF and ' + lig_name + '\n')		# new title
assembly.write(sys_lig_number + '\n')											# updated atom number
for s in sys_lines[2:-1]:
	assembly.write(s)
for f in lcent_lines[2:-1]:
	assembly.write(f)
assembly.write(sys_lines[-1])
assembly.close()
system.close()
ligand.close()

### FOR VMD USAGE (DEPRECATED) ##################################

# strip waters with VMD and compile the final system
#strip_file = open ('strip.tcl', 'w')
#if lcent_lines[-2][5:8] != lig_name:
#	strip_file.write('set strip [atomselect top "not same residue as (water and within 2 of (resname ' + str(lig_name) + ' ' + str(lcent_lines[-2][5:7]) + '))"]\n')
#else:
#	strip_file.write('set strip [atomselect top "not same residue as (water and within 2 of resname ' + str(lig_name) + ')"]\n')
#lig_centered.close()
#strip_file.write('$strip writepdb system_stripped.pdb\n')
#strip_file.write('exit\n')
#strip_file.close()
#os.popen('vmd -dispdev text system_assembled.gro < strip.tcl')
#os.popen('@gromacs2016 gmx editconf -quiet -f system_stripped.pdb -o system_stripped.gro', shell=True)

#################################################################

#### strip files (gmx select) ########
strip_file = open('strip.txt', 'w')
if lcent_lines[-2][5:8] != lig_name:
	strip_file.write('group "System" and not same residue as (group "SOL" and within 0.2 of resname ' + str(lig_name) + ' ' + str(lcent_lines[-2][5:7]) + ')\n')
else:
	strip_file.write('group "System" and not same residue as (group "SOL" and within 0.2 of resname ' + str(lig_name) + ')\n')
lig_centered.close()
strip_file.close()
call(gromacs + "/bin/gmx -quiet select -s system_assembled.gro -on index_new.ndx -sf strip.txt", shell=True)
call(gromacs + "/bin/gmx -quiet editconf -f system_assembled.gro -n index_new.ndx -o system_stripped.gro", shell=True)

#### get number of waters ############
stripped = open('system_stripped.gro')
strip_lines = stripped.readlines()
waters = 0
for st in strip_lines:
	if st[5:8] == 'SOL':
		waters = waters + 1
total_waters = str(waters / 3)

#### create topol.top and index.ndx ##
top_file = open('topol.top', 'r')
top_lines = top_file.readlines()
topol_new = open('topol_new.top', 'w')
for t in top_lines[1:45]:
	topol_new.write(t)
topol_new.write('\n; Include ligand topology\n#include "' + str.lower(lig_topol) + '"\n')		# add ligand topology
for t in top_lines[45:60]:
	topol_new.write(t)
topol_new.write('SOL\t\t\t\t' + str(total_waters) + '\n')										# add new solvent
for t in top_lines[61:]:
	topol_new.write(t)
topol_new.write(str(lig_name) + '\t\t\t1\n')													# write ligand name
if strip_lines[-2][5:8] != lig_name:															# check for counter-ions at the bottom of the file
	topol_new.write(str(strip_lines[-2][5:7]) + '\t\t\t1\n')

#### create new index ################
index_auto = open('index_auto.txt', 'w')
index_auto.write('26|27|28|29|30\n')
index_auto.write('name 36 OM\n')
index_auto.write('1|36\n')
index_auto.write('1 & ri 1-340\n')
index_auto.write('name 38 monomer_A\n')
index_auto.write('q\n')
index_auto.close()
call(gromacs + "/bin/gmx -quiet make_ndx -quiet -f system_stripped.gro -o index.ndx < index_auto.txt", shell=True)
stripped.close()
top_file.close()
topol_new.close()

#### delete unneeded files ###########
os.system('rm system_stripped.pdb lig_centered.gro index_auto.txt strip.txt')

#### generate mdp and tpr files ######
mdp = open('npt_pull_vec.mdp', 'w')
mdp.write('define = \n')
mdp.write('integrator = md \n')
mdp.write('nsteps = 50000000 \n')
mdp.write('dt = 0.002 \n')
mdp.write('nstxout = 500000 \n')
mdp.write('nstvout = 500000 \n')
mdp.write('nstenergy = 50000 \n')
mdp.write('nstlog = 50000 \n')
mdp.write('nstxout-compressed = 50000 \n')
mdp.write('continuation = yes \n')
mdp.write('constraint_algorithm = lincs \n')
mdp.write('constraints = all-bonds \n')
mdp.write('lincs_iter = 1 \n')
mdp.write('lincs_order = 4 \n')
mdp.write('ns_type = grid \n')
mdp.write('rlist = 1.2 \n')
mdp.write('rcoulomb = 1.2 \n')
mdp.write('rvdw = 1.2 \n')
mdp.write('vdwtype = Cut-off \n')
mdp.write('vdw-modifier = Force-switch \n')
mdp.write('rvdw_switch = 1.0 \n')
mdp.write('coulombtype = PME \n')
mdp.write('pme_order = 4 \n')
mdp.write('fourierspacing = 0.16 \n')
mdp.write('cutoff-scheme = Verlet \n')
mdp.write('nstlist = 40 \n')
mdp.write('tcoupl = Nose-Hoover \n')
mdp.write('tc-grps = Protein_OM ' + str(lig_name) + ' Water_and_ions \n')
mdp.write('tau_t = 1.0 1.0 1.0 \n')
mdp.write('ref_t = 310 310 310 \n')
mdp.write('nsttcouple = 1 \n')
mdp.write('nhchainlength = 1 \n')
mdp.write('pcoupl = Parrinello-Rahman \n')
mdp.write('pcoupltype = semiisotropic \n')
mdp.write('tau_p = 5.0 \n')
mdp.write('ref_p = 1.0 1.0 \n')
mdp.write('compressibility = 4.5e-5 0 \n')
mdp.write('pbc = xyz \n')
mdp.write('DispCorr = EnerPres \n')
mdp.write('gen_vel = no \n')
mdp.write('nstcomm = 1 \n')
mdp.write('nstcalcenergy = 1 \n')
mdp.write('comm-mode = Linear \n')
mdp.write('comm-grps = Protein_OM ' + str(lig_name) + ' Water_and_ions \n')
mdp.write('pull = yes \n')
mdp.write('pull_ngroups = 1 \n')
mdp.write('pull_ncoords = 1 \n')
mdp.write('pull_group1_name = ' + str(lig_name) + '\n')
mdp.write('pull_nstxout = 10 \n')
mdp.write('pull_nstfout = 10 \n')
mdp.write('pull_coord1_type = umbrella \n')
mdp.write('pull_coord1_geometry = direction-periodic \n')
mdp.write('pull_coord1_origin = 6.200 8.550 8.950 \n')
mdp.write('pull_coord1_groups = 0 1 \n')
mdp.write('pull_coord1_vec = 0 0 -1 \n')
mdp.write('pull_coord1_rate = 0.0000875 \n')
mdp.write('pull_coord1_k = 1000 \n')
mdp.write('pull_coord1_start = no\n')
mdp.close()

if not os.path.exists('PULL'):
	os.makedirs('PULL')
call(gromacs + "/bin/gmx -quiet grompp -f npt_pull_vec.mdp -c system_stripped.gro -p topol_new.top -n index.ndx -o ./PULL/OmpF_" + str(lig_name) + ".tpr -maxwarn 1", shell=True)

#### optionally run GROMACS job ######
start = raw_input("Do you want to start GROMACS mdrun (y/n)? ")
if start == "y":
	os.chdir('./PULL')
	call(gromacs + "/bin/gmx -quiet mdrun -s *.tpr -nb gpu", shell=True)

#### exit the script #################
quit()
