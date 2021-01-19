Experimental data:

1) The program 'processing_raw_data.m' was run first. This program performs the image analysis, 
cell tracking, and some other calculations, e.g. msd  and the cell positions as a function of radius.
It required only one input parameter, which is the directory with the parameters file.
This is a single .mat file that contains all the parameter needed for the program to run.
All the details for the parameter file can be found in the header of the 'processing_raw_data.m'.

2) The 'density_vs_velcoty.m' is the post-processing program, that calculates the local velocity and density as a function of radius and time.
It only requires one input with the directory where the proccessed saved files from step 1 are located.


Simulations:
This code assumes that the raw data `results_rho0_*.txt` is contained in a subfolder named `data`.