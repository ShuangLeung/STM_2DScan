'''
                                  =============
                                  STM-2DScan.py
                                  =============

---------------------------------------------------------------------------------
Introduction:

STM-2DScan.py is a postprocessing script for VASP code to generate STM images 
based on DFT-calculations. It firstly imports volumetric data from a file with 
CHGCAR format (e.g. CHGCAR, PARCHG and CURRENT (output data file of BSKAN)...) 
and then interpolates the data onto specified two-dimensional (2D) slice in an 
arbitrary xy-plane, which will be saved as a '.png' format file.

This script is a free software, which means that you can redistribute it or 
modify it under the terms of the GNU General Public License. Of course, any 
feedback of modifications, bug fixes and improvements is appreciate.

---------------------------------------------------------------------------------
Usage:

Step1. Setting the height of tip position. For convenience, an alternative 
       option is to set a range of heights in order to generate a series of 
       corresponding 2D slices at a time. 

Step2. Selecting a scan mode to obtain the 2D slice for visualization. There 
       are three types optional scan modes, including two common types: the 
       constant-height mode and the constant-current mode. The working 
       mechanism of these two modes are similar to the ways of stm.scan and 
       stm.scan2 functions in ASE (You can visit the website of ASE to learn 
       about more details). More simply, you also can select the third mode 
       to make the 2D slice by cutting a xy-plane at a given height from 3D 
       charge density. Although this mode produces almost the same result as 
       that of the constant-height mode, allowing for more general applica-
       tions, e.g. display a 2D image of CHGCAR or CURRENT (Of course, VESTA 
       or other software can do it better).

Step3. Optionally, you can replicate the unit-cell by inputting the cell ranges
       along a and b vectors, respectively.

Step4. Plot the 2D contour by using matplotlib and save as a .png format image.
       Multifarious colormaps included with matplotlib are available, allowing 
       you to choose your most favorite one.

---------------------------------------------------------------------------------
Important Notes: 

1. This script runs in the Python3 environment.

2. The script depends on the following python modules and packages:
    numpy - http://numpy.scipy.org/                                  
    ase - https://wiki.fysik.dtu.dk/ase/            
    matplotlib - http://matplotlib.sourceforge.net/
    
   If your OS doesn’t have above packages installed, you can install 
   them using a command like:
    pip install --upgrade numpy ase matplotlib

3. Multifarious colormaps included with matplotlib are alternative. 
   See https://matplotlib.org/examples/color/colormaps_reference.html
   for more details.

4. Reading the charge density file is a time-consuming progress, which may 
   cost a few seconds (even reach to the order of minutes).

---------------------------------------------------------------------------------
Author: ShuangLeung (sleung1924@gmail.com)                         
Date of last version: 2019/08/30
---------------------------------------------------------------------------------
'''

import sys, os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from ase.calculators.vasp import VaspChargeDensity

print("-"*86)
starttime = time.perf_counter()
print("Starting the program at")
print(time.strftime("%Y-%m-%d %H:%M:%S"))
print("-"*86)

# Check if the PARCHG exists and import charge density data from it.
# If not, enter the filename of CHGCAR/PARCHG format file containing charge density data.
if os.path.exists("PARCHG"):
	print("PARCHG file already exists, reading data now...\n") 
	vasp_charge = VaspChargeDensity('PARCHG')
	density = vasp_charge.chg[-1]
	atoms = vasp_charge.atoms[-1]
	del vasp_charge
	print("Done!\n")

else:
	active = True
	while active:
		inputstr1 = input("Enter the filename of CHGCAR/PARCHG format file containing data:") 
		if os.path.exists(inputstr1):
			print("\nReading charge density data from %s...\n" %inputstr1)
			vasp_charge = VaspChargeDensity(filename=inputstr1)
			density = vasp_charge.chg[-1]
			atoms = vasp_charge.atoms[-1]
			del vasp_charge
			print("Done!\n")
			break
		else:
			print("\n*** ERROR: The input filename doesn't exist!!! Please retype the filename again. ***\n")
			active = True

# Read size of XYZ grids.
ngridpoints = np.array(density.shape)

# Read scaling factor and unit cell of the crystal structure.
unit_cell = atoms.get_cell()
cell_lengths = np.sqrt(np.dot(unit_cell,unit_cell.transpose()).diagonal())

# Select the way of setting the height of tip position to make the 2D slice.
inputstr2 = input("Select the way of setting the height of tip position:\n %s %s %s" % ("1: Setting a specified height;[Default option]\n",
												"2: Specifying a range of heights.\n",
												"Enter your option [1/2]:"))
height_option = inputstr2.strip()
# Set the default option.
height_options = ["1","2"]
if height_option not in height_options:
	print("\n*** WARNING: the input has syntax errors. Select the defaul option automatically. ***\n")
	height_option = "1"
# Initialize the input height of tip position.
heights = []

# Select the scan mode to get the visualization of 2D slice
inputstr3 = input("\nSelect the STM scan mode to obtain 2D slice for visualization:\n %s %s %s %s" % ("1: Constant-height mode;[Default option]\n",
												"2: Constant-current mode;\n",
												"3: 2D-slice at a specified height.\n",
												"Enter your option [1/2/3]:"))
scan_mode_option = inputstr3.strip()
# Set the default option.
scan_mode_options = ["1","2","3"]
if scan_mode_option not in scan_mode_options:
	print("\n*** WARNING: the input has syntax errors. Select the constant-height mode automatically. ***\n")
	scan_mode_option = "1"

# Input the height(s) (in Angstroms) of tip position to make the 2D slice(s).
# Find the integer plane(s) (scan_mode_option==1/2) or 
# the closest plane(s) (scan_mode_option==3) corresponding to the input height(s).
if height_option == "1" :
	while True:
		try:
			inputstr4 = input("\nEnter the height in Angstroms along c vector to make 2D slice：")
			height = float(inputstr4.strip())
			heights.append(height)
			break
		except:
			print("\n*** ERROR: The input is non-digit!!! Please retype again. ***\n")

	print("\nAttempting to find a plane/grid-point at a height of %4f Angs.\n" % height)

else:
	while True:
		try:
			inputstr5 = input("\nEnter the initial height in Angstroms along c vector to make 2D slice:")
			height1 = float(inputstr5.strip())
			break
		except:
			print("\n*** ERROR: The input is incorrect!!! Please retype again. ***\n")
	while True:
		try:
			inputstr6 = input("Enter the final height in Angstroms along c vector to make 2D slice:")
			height2 = float(inputstr6.strip())
			break
		except:
			print("\n*** ERROR: The input is incorrect!!! Please retype again. ***\n")

	# Find the closest planes corresponding to the input heights.
	if scan_mode_option == "3":
		plane_index1 = int(round(height1/cell_lengths[2]*ngridpoints[2]))%ngridpoints[2]
		plane_index2 = int(round(height2/cell_lengths[2]*ngridpoints[2]))%ngridpoints[2]

	# Find the integer planes (grid points) corresponding to the input heights.
	else:
		plane_index1 = int(height1/cell_lengths[2]*ngridpoints[2])%ngridpoints[2]
		plane_index2 = int(height2/cell_lengths[2]*ngridpoints[2])%ngridpoints[2]

	# Determine the maximal linear insertion points.
	max_inspoint = abs(plane_index2-plane_index1)+1

	while True:
		try:
			inputstr7 = input("Enter the amount of linear insertion points(max_inspoints <= %d!) to divide the height:" %max_inspoint)
			point = int(inputstr7)
			if point > max_inspoint:
				point = max_inspoint #Set the default amount of linear insertion points
				for i in range(0,point+2):
					h = abs(height2-height1)*i/(point+1)+min(height1,height2)
					heights.append(h)
			else:
				for i in range(0,point+2):
					h = abs(height2-height1)*i/(point+1)+min(height1,height2)
					heights.append(h)		
			break
		except:
			print("\n*** ERROR: The input is incorrect!!! Please retype again. ***\n")

	print("\nAttempting to find a serials of planes/grid-points at the heights ranging from %4f to %4f Angs.\n" % (height1,height2))

# Make a supercell.
while True:
	try:
		inputstr8 = input("Enter the supercell range along a vector:")
		m = int(inputstr8)
		break
	except:
		print("\n*** ERROR: The input is not an integer!!! Please input again. ***\n")

while True:
	try:
		inputstr9 = input("Enter the supercell range along b vector:")
		n = int(inputstr9)
		break
	except:
		print("\n*** ERROR: The input is not an integer!!! Please input again. ***\n")

# The size of grids along a/b vector.
supercell_xngridpoints = (ngridpoints[0]-1)*m+1
supercell_yngridpoints = (ngridpoints[1]-1)*n+1

# Make arrays of x and y values.
supercell_xarray = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float)
supercell_yarray = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float)

# Make arrays of supercell_density2D/I/H with the same dimensions as x/y arrays.
supercell_density2D = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float)
I = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float)
H = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float)

#Find projection of b vector onto a vector.
ytox = np.dot(unit_cell[0],unit_cell[1].T)/cell_lengths[0]
#Find component of b vector perpendicular to a vector.
ynormal = np.cross(unit_cell[0],unit_cell[1].T)/cell_lengths[0]
ynormal = np.sqrt(np.dot(ynormal,ynormal.T))

# Plot the STM-2Dscan images in matplotlib.
# Multifarious colormaps included with matplotlib are alternative.
colormaps = ['afmhot','gray','bone','gist_heat','gnuplot','gnuplot2','CMRmap','jet','viridis','plasma','inferno','magma','cividis']
print("\nSelect the colormap to plot the STM-2Dscan images:\n",
	" 0: afmhot\n",
	" 1: gray\n",
	" 2: bone\n",
	" 3: gist_heat\n",
	" 4: gnuplot\n",
	" 5: gnuplot2\n",
	" 6: CMRmap\n",
	" 7: jet\n",
	" 8: viridis\n",
	" 9: plasma\n",
	"10: inferno\n",
	"11: magma\n",
	"12: cividis")

inputstr10 = input('Enter your option [0-12]:')
colormaps_option = inputstr10.strip()
colormaps_options = []
for i in range(0,13):
	colormaps_options.append(str(i))

if colormaps_option in colormaps_options:
	cmap_No = int(colormaps_option)
else:
	print("\n*** WARNING: the input has syntax errors. Select the default option automatically. ***\n")
	cmap_No = 0

print("\nGenerating STM-2Dslice images...\n")
for h in heights:
	# Contant height mode.
	if scan_mode_option == "1":
		n1 = h/cell_lengths[2]*ngridpoints[2]
		dn1 = n1-np.floor(n1)
		n1 = int(n1)%ngridpoints[2]
		ldos = (1-dn1)*density[:,:,n1]+dn1*density[:,:,(n1+1)%ngridpoints[2]]

	# Contant current mode.
	elif scan_mode_option == "2":
		n2 = h/cell_lengths[2]*ngridpoints[2]
		dn2 = n2-np.floor(n2)
		n2 = int(n2)%ngridpoints[2]
		# Get the averaged current.
		averaged_current = ((1-dn2)*density[:,:,n2].mean()+dn2*density[:,:,(n2+1)%ngridpoints[2]].mean())
		c1 = density[:,:,n2]
		c2 = density[:,:,(n2+1)%ngridpoints[2]]
	
	# 2D-slice at a specified height.
	else:
		plane_index = int(round(h/cell_lengths[2]*ngridpoints[2]))%ngridpoints[2]
		density2D = density[:,:,plane_index]

	for i in range(supercell_xngridpoints):
		for j in range(supercell_yngridpoints):
			supercell_xarray[i][j] = float(i)/float(supercell_xngridpoints)*cell_lengths[0]*m+float(j)/float(supercell_yngridpnoints)*ytox*
			supercell_yarray[i][j] = float(j)/float(supercell_yngridpoints)*ynormal*n
			mi = i%(ngridpoints[0]-1)
			nj = j%(ngridpoints[1]-1)

			if scan_mode_option == "1":
				I[i][j] = ldos[mi][nj]

			elif scan_mode_option == "2":
				if c2[mi][nj]-c1[mi][nj] == 0:
					H[i][j] = n1*cell_lengths[2]/ngridpoints[2]

				else:
					H[i][j] = (n2+(averaged_current-c1[mi][nj])/(c2[mi][nj]-c1[mi][nj]))*cell_lengths[2]/ngridpoints[2]
				
			else:
				supercell_density2D[i][j] = density2D[mi][nj]

	# Plot the 2D contour in matplotlib. 
	# The newly generated images will be named after the letter "H/C/S"(constant-height mode/constant-current mode/2D-slice)
	# and height value (of tip position).
	if scan_mode_option == "1":
		P = I
		mode_label = "H" 
	elif scan_mode_option == "2":
		P = H
		mode_label = "C"
	else:
		P = supercell_density2D
		mode_label = "S"

	plt.figure()
	plt.rcParams['figure.max_open_warning'] = 50
	plt.gca(aspect='equal')
	plt.axis('off')
	plt.xticks(())
	plt.yticks(())
	cm = plt.cm.get_cmap('%s' %colormaps[cmap_No])
	plt.contourf(supercell_xarray,supercell_yarray,P, 900, cmap=cm)
	plt.colorbar()
	plt.savefig(mode_label+'_'+str(round(h,3))+'.png', dpi=300, bbox_inches='tight')

print("Done!\n")

print("-"*86)
starttime = time.perf_counter()
print("Ending program at")
print(time.strftime("%Y-%m-%d %H:%M:%S"))
print("-"*86)
