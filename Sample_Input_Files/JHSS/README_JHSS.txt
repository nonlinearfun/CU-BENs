README_model: This folder contains the input file of a ship hull (JHSS). Detailed description of the model can be found in JHSS model test information - misc.pdf. The ship hull is modeled with 91,776 DKT shell elements. The model is constrained in rotations as well as y- and z- direction. Nonlinear pressure pulse loading is imposed on the stern to simulate power plant vibration. The total time of this analysis is 0.012 sec, divided into 30 equal time steps. Displacements in the longitudinal axis are computed. This is a first order elastic analysis that uses the Newmark implicit time integration scheme. 

Usage:

Move input file to same folder as all other source code files. 

Rename input file to “model_def.txt”. 