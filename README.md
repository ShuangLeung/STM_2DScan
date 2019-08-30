# STM_2DScan

----------------------------------------------------------------------------------------------------------------------------------------
Introduction:

STM-2DScan.py is a postprocessing script for VASP code to generate STM images based on DFT-calculations. It firstly imports volumetric data from a file with CHGCAR format (e.g. CHGCAR, PARCHG and CURRENT (output data file of BSKAN)...) and then interpolates the data onto specified two-dimensional (2D) slice in an arbitrary xy-plane, which will be saved as a '.png' format file.

This script is a free software, which means that you can redistribute it or modify it under the terms of the GNU General Public License. Of course, any feedback of modifications, bug fixes and improvements is appreciate.

----------------------------------------------------------------------------------------------------------------------------------------
Usage:

Step1. Setting the height of tip position. For convenience, an alternative option is to set a range of heights in order to generate a          series of corresponding 2D slices at a time. 

Step2. Selecting a scan mode to obtain the 2D slice for visualization. There are three types optional scan modes, including two common          types: the constant-height mode and the constant-current mode. The working mechanism of these two modes are similar to the ways          of stm.scan and stm.scan2 functions in ASE (You can visit the website of ASE to learn about more details). More simply, you also        can select the third mode to make the 2D slice by cutting a xy-plane at a given height from 3D charge density. Although this mode        produces almost the same result as that of the constant-height mode, allowing for more general applications, e.g. display a 2D          image of CHGCAR or CURRENT (Of course, VESTA or other software can do it better).

Step3. Optionally, you can replicate the unit-cell by inputting the cell ranges along a and b vectors, respectively.

Step4. Plot the 2D contour by using matplotlib and save as a .png format image. Multifarious colormaps included with matplotlib are            available, allowing you to choose your most favorite one.

----------------------------------------------------------------------------------------------------------------------------------------
Important Notes: 

1. This script runs in the Python3 environment.

2. The script depends on the following python modules and packages:                                                                     
    numpy - http://numpy.scipy.org/                                  
    ase - https://wiki.fysik.dtu.dk/ase/            
    matplotlib - http://matplotlib.sourceforge.net/
    
   If your OS doesn’t have above packages installed, you can install them using a command like:
    pip install --upgrade numpy ase matplotlib

3. Multifarious colormaps included with matplotlib are alternative. See https://matplotlib.org/examples/color/colormaps_reference.html      for more details.

4. Reading the charge density file is a time-consuming progress, which may cost a few seconds (even reach to the order of minutes).

----------------------------------------------------------------------------------------------------------------------------------------
Authors: ShuangLeung                                                                                           Email:sleung1924@gmail.com                         
Date of last version: 2019/08/30

----------------------------------------------------------------------------------------------------------------------------------------
