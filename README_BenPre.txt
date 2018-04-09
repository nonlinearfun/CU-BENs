README_BenPre: BenPre.py is a rudimentary pre-processing routine for CU-BEN. This routine was written for Python 2.7. It operated on files that conform to the abaqus.inp formatting, along with supplemental command line prompts, to create an input file (i.e. model_def.txt) for CU-BEN. 

This routine is ONLY intended to build models containing three node triangular SHELL elements (to be used in FE context). This code will not prepare a CU-BEN input file for other elements, or for a CU-BEN acoustic FSI analyses.

While this code is helpful in constructing CU-BEN input files, it is HIGHLY recommended that the user review the input files generated: in order to ensure that she/he obtains the desired model specification. 

Usage:

Execute BenPre.py in the same directory as the CU-BEN.inp file, created with ABAQUS (R). The input file generated from ABAQUS MUST be named CU-BEN.inp. 

To execute BenPre.py, issue the command "python BenPre.py" (without quotes) within the command line.

The user will be prompted to enter the analysis type, solution algorithm type, solver algorithm type, and to choose to “execute node-renumbering” in the command line. Options for each of the analysis (etc) types will be provided and the user is to enter the numerical flag associated with their choice. The code then will access the ABAQUS input file to transcribe the geometry established in ABAQUS to be read by CU-BEN. Depending on the choice of analysis/algorithm, the user will be prompted to enter additional parameters for the solution. If the prompt requires more than one value, a guide will be printed below the prompt to help the user input the values correctly. 

A tutorial on how to correctly set up the ABAQUS file, to be read seamlessly by BenPre, can be found in the "Introduction to CU-BEN" document. The primary concerns to consider when constructing an ABAQUS model are to use only element type S3 shell elements and to title the .inp file "CU-BEN.inp". If running a dynamic analysis, it is important to explicitly define a time step and enter the load history in a tabular form. The ABAQUS input deck card, *step, created for the time steps must be called "TimeSteps". If the load on any nodes is not zero at time t=0, or if there are any initial prescribed displacement, velocity, or accelerations on any of the nodes, then this must be defined in a Step called "Time0". More details can be found in the "Introduction to CU-BEN" document found at http://earls.cee.cornell.edu/software/. 

The code will output a file named "model_def.txt". This can be used as an input file in CU-BEN. 