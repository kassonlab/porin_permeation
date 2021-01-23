#!/usr/bin/python

import os, sys, numpy, glob, time
from decimal import Decimal
import matplotlib.pyplot as plt

'''
Python script to process umbrella sampling runs to obtain diffusion coefficients and PMF profiles

This script accompanies the following publications:
Simulation-guided engineering of antibiotics for improved bacterial uptake, Ricardo Ferreir
a and Peter Kasson, doi:10.1101/2020.10.08.330332.
Antibiotic Uptake Across Gram-Negative Outer Membranes: Better Predictions Towards Better Antibiotics, Ricardo Ferreira and Peter Kasson, doi:10.1021/acsinfecdis.9b00201.

The program takes two arguments, initial and final sampling times, from sys.argv[1] and sys.argv[2].
In the current version, it only takes sampling each 10ns, which was proved to be sufficient 
to calculate permeabilities. Thus:
	
	>>> python 3_analyze_runs_vXx.py

With this script we intent to:
	- calculate the variance at each bin (from pullx files);
	- generate autocorrelation times (IACT) and, with variances, diffusions in each bin (Marrink & Berendsem;
	- calculate resistance to permeation in each bin acording to the homogeneous solubility-diffusion model;
	- generate a permeability coefficient (P) as the inverse of permeation resistance;
	- estimate a standard deviation for P, by bootstrapping PMFs

### version 0.1.2 (v5) ###########
	- initial and final time are now inserted interactively;
	- resistance is now calculated from force integration within each bin (and not from the PMFs);
	- add resistance plotting (now plots diffusions, PMF bootstrap and resistance);

### version 0.1.1 (v4 and v4a) ###
	- autocorrelation times (IACT) are used in bootstrapping PMFs (to certify that random PMFs are uncorrelated);
	- add bootstrap function to generate 100 PMF profiles;
	- resistance is now calculated as the average of the bootstrapped PMFs;
	- resistance scalling factor was removed from calculations (not needed);

### version 0.1.0 (v3) ###
	- remove dependence from GROMACS g_wham built-in module;
	- calculate PMF from average forces in pullf files;
	- calculate integrated autocorrelation times (IAC) through the approach by Zhu & Hummer (DOI: 10.1021/ct2009279)
	  from block averaging method by Flyvbjerg & Petersen (DOI:10.1063/1.457480);
	- get bin centers from "output_bins.txt", generated with "2_get_us_windows.py" script;
	- print bin center, diffusions and pmf's to text file, along with diff and pmf plots (png);

### version 0.0.2 (v2) ###
	- added mdp options to script (generate GROMACS 4.6.x tpr files)
	- script can now be run from main directory

### version 0.0.1 (v1) ###	
	- must use GROMACS 4.6.7 version due to 'position' keyword in g_wham
	- assume that the general name for pull files is 'pullx_US_*.xvg' and 'pullf_US_*.xvg'
	- assumes that all necessary files are in the working directory
	- assumes that pull-nstxout and pull-nstfout (in mdp file) is set to 10
'''
########### define functions ##########
def format_e(n):
	# round the final P value
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]
    
def blockAverage(datastream, maxBlockSize=0):									
	# block averaging function (needed for IACT calculation)
	Nobs = len(datastream)
	minBlockSize = 1;
	if maxBlockSize == 0:
		maxBlockSize = int(Nobs/5); 
	NumBlocks = maxBlockSize - minBlockSize
	blockMean = numpy.zeros(NumBlocks)
	blockVar = numpy.zeros(NumBlocks)
	blockCtr = 0
	for blockSize in range(minBlockSize, maxBlockSize):
		Nblock = int(Nobs/blockSize)
		obsProp = numpy.zeros(Nblock)
		for i in range(1,Nblock+1):			
			ibeg = (i-1) * blockSize
			iend =  ibeg + blockSize
			obsProp[i-1] = numpy.mean(datastream[ibeg:iend])
		# blockMean[blockCtr] = numpy.mean(obsProp) --> only used if we need averages
		blockVar[blockCtr] = numpy.var(obsProp)/(Nblock - 1)
		blockCtr += 1
	return blockVar[-1]
	
def BootStrap (x,step=1000):
	if(numpy.mean(x)) == 0:
		return 0.000, 0.000
	else:
		avg =[]
		x = numpy.array(x)
		n = len(x)
		idx = numpy.random.randint(0,n,(step,n))
		sample_x = x[idx]
		avg = numpy.sort(numpy.mean(sample_x,1))
		return avg
		
### initialize script header ##########
os.system('clear')
print('*********** 3_analyze_us_runs.py script ***********')
print('                                                   ')
print('          PYTHON SCRIPT TO CALCULATE PORIN         ')
print('         PERMEABILITIES FROM US SIMULATIONS        ')
print('                                                   ')
print('***************************************************\n')

### initialize lists and variables ###
positions = []						 				# for storing bin centers (from output_bins.txt)
variances = []						 				# for storing variances (from pullx files)
iact = []						     				# for storing autocorrelations times
diffusions = []						 				# for storing diffusions (variances * iact)

tau = 0.0							 				#
init1 = int(raw_input("Enter starting time (ns): "))#		
init2 = int(raw_input("Enter ending time (ns)  : "))# 
init3 = int(raw_input("How many monomers  : "))     # 
diffusion = 0.0						 				# 
exponential = 0.0					 				#
integral = 0.0		 				 				#

meanff = []							 				# forces read from pullf.xvg files
forces_array = []					 				# for storing arrays of forces from bootstrap 
pmf_array = []						 				# for storing bin energies for resistance calculations
pmf_plot = []										# for plotting final PMFs
resistance_array = [0]				 				# for storing final resistance (R) sums

#### define starting points ##########
# defines start and end, accounting for file header (18 lines)
start = int((init1 / 0.00002) + 19)
end = int((init2 / 0.00002) + 19)

#### create directories ##############
if not os.path.exists('./US/WHAM/PCALC'):
    os.makedirs('./US/WHAM/PCALC')
os.chdir('./US')

# get average positions from output_bins.txt file:
fx = open("output_bins.txt","r")
fxlines = fx.readlines()
for l in fxlines[1:]:
		time, position = l.split()
		positions.append(float(position))
fx.close()

### go to target directory ###########
os.chdir('./WHAM/PCALC')

#### open and sort files by name #####
pullx = sorted(glob.glob('../*pullx*.xvg'))
pullx.sort(key=lambda f: int(filter(str.isdigit, f)))
pullf = sorted(glob.glob('../*pullf*.xvg'))
pullf.sort(key=lambda f: int(filter(str.isdigit, f)))

#### verify if all files are present ##
if len(positions) != (len(pullx) or len(pullf)):
	os.system('clear')
	print('******************** WARNING **********************')
	print('                                                   ')
	print('        The number of bins does not match          ')
	print('      the number of US files. Please verify.       ')
	print('                                                   ')
	print('***************************************************\n')
	sys.exit()
	
### get variances and IACTs ###########
p = 0
for xx in pullx:
	print("Calculating variance and IACT in bin " + str(p + 1) + ", z = " + str(round(positions[p],3)))
	x = []
	pullfile = open(xx, 'r')
	pulllines = pullfile.readlines()
	# 1 ps for each 50 lines
	for l in pulllines[start:end:50]:
		time, z = l.split('\t')
		# z in nm
		x.append(float(z))
	# calculate and append variances
	variances.append(numpy.var(x))													
	# calculate and append autocorrelation times
	tau = 0.0
	dt = 1.0
	# tau in ps
	tau = (((len(x)* blockAverage(x, maxBlockSize=0)) / numpy.var(x))-1) * (dt / 2)
	iact.append(tau)
	pullfile.close()
	p = p + 1
	
### get PMF from average forces #######
b = 0
for f in pullf:
	print("Calculating and bootstrapping PMFs, bin " + str(b + 1) + ", IACT = " + str(int(iact[b])))
	meanff = []
	ff = open(f,"r")
	fflines = ff.readlines()
	# decorrelate forces for generating PMFs
	# 1 ps each 50 lines
	i = int(round(iact[b]) * 50)
	# to prevent when IACT = 0 ps
	if i == 0:
		i = 1
	for lines in fflines[start:end:i]:
		time, energy = lines.split('\t')
		# forces in kJ/mol/nm
		meanff.append(float(energy))
	forces_array.append(BootStrap(meanff,step=100))
	b = b + 1

b = 0
e = len(forces_array)
while b < e:
	point = 0.0
	pointp = 0.0
	pmf = []
	pmff = []
	new_list = zip(*forces_array)[b]
	bb = 1
	ee = len(new_list)
	pmf.append(new_list[0])
	pmff.append(new_list[0])
	point = new_list[0]
	pointp = new_list[0]
	while bb < ee:
		# bin centers in nm, forces in kJ/mol/nm
		# to use in resistance calculations
		point = (positions[bb] - positions[bb - 1]) * new_list[bb]
		pmf.append(point)
		# to construct the final PMFs
		pointp += (positions[bb] - positions[bb - 1]) * new_list[bb]
		pmff.append(pointp)
		bb = bb + 1
	# add to corresponding lists	
	pmf_array.append(pmf)
	pmf_plot.append(pmff)
	b = b + 1

### calculate D from var and IACT #####
b = 0
e = len(positions)
while b < e:
	# diffusion in cm^2/s
	diffusion = ((float(variances[b]) * (10 ** (-14))) / (float(iact[b]) * (10 ** (-12))))
	diffusions.append(float(diffusion))
	b = b + 1
	
### calculate R from Diff and PMFs ####
b = 0
e = len(positions) - 1
while b < e:
	resistance = []
	for item in pmf_array:
		# gas constant R (in kJ/mol.K) and T (in Kelvin)
		integral = ((positions[b] - positions[b + 1]) * (10 ** (-7))) * ((numpy.exp((item[b]) / (0.0083144598 * 310)) / diffusions[b]))
		# in cm/s (average PMF)
		resistance.append(integral)
	# append the resistance sum (to allow averages and stdev)
	resistance_array.append(numpy.sum(resistance))
	b = b + 1
	
### estimate mean and stdev ###########
mean = (init3 / numpy.mean(resistance_array))
stdev = (1 / numpy.std(resistance_array))

### write results/data to file ########
data_file = open('output_results_' + str(init1) + '-' + str(init2) + '.dat', 'w')
data_file.write('Bin center\tDiff\t\t\tPMF\n')
b1 = 0
e1 = len(positions)
while b1 < e1:
	data_file.write(str(round(positions[b1],3)) + \
	'\t\t' + str(format_e((float(variances[b1]) * (10 ** (-14))) / (float(iact[b1]) * (10 ** (-12))))) + \
	'\t\t' + str(round(pmf[b1],3)) + '\n')
	b1 = b1 + 1
data_file.write('\nCalculated permeation coefficient is ' + str(format_e(Decimal(mean))) + ' +/- ' + str(format_e(Decimal(stdev))) + ' cm/s\n')
data_file.close()

### print out D, PMF and R curves #####
fig1 = plt.figure()
ax1 = fig1.add_subplot(311)
ax1.set_title(str(format_e(Decimal(mean))) + ' +/- ' + str(format_e(Decimal(stdev))) + ' cm/s')
ax1.plot(positions, diffusions)
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('z (nm)')
plt.ylabel('D (cm^2/s)')
ax2 = fig1.add_subplot(312)
for i in pmf_plot:
	ax2.plot(positions, i)
plt.xlabel('z (nm)')
plt.ylabel('PMF (kJ/mol)')
ax3 = fig1.add_subplot(313)
ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax3.plot(positions, resistance_array)
plt.xlabel('z (nm)')
plt.ylabel('R (cm^(-1)/s)')

print('\nCalculated permeation coefficient is ' + str(format_e(Decimal(mean))) + ' +/- ' + str(format_e(Decimal(stdev))) + ' cm/s')

plt.tight_layout()
plt.savefig('output_plots_' + str(init1) + '-' + str(init2) + '.svg', format="svg")
plt.show()

### exit the script ###################
quit()
