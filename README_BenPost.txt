README_BenPost: BenPost.c is a rudimentary post processor for CU-BEN.

Usage:

Keep in same folder as CU-BEN's results*.txt output files

After creating and running BenPost executable, a "master.pvd" file will output for ParaView (http://www.paraview.org/) accessibility

The BenPost executable also outputs a folder named "vtu_files" which contains a "slave_*.vtu" file for each time step of a given CUBEN simulation

This folder can be opened in ParaView and may subsequently be saved as a .vtk file
