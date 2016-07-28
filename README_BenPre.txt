README_BenPre: BenPre.py is a rudimentary pre-processing code for CU-BEN. This code was written in Python 2.7. It works with abaqus .inp files and command line prompts to create an input file for CU-BEN. 

This code is ONLY to build models containing ONLY SHELL elements to be used in FE context. This code will not prepare a CU-BEN input file for other elements or for FSI analysis.

While this code is helpful in constructing CU-BEN input files, it is HIGHLY recommended that the user reviews the input files generated. 

Usage:

Keep in same folder as CU-BEN.inp file created in Abaqus. The input file generated from Abaqus MUST be named CU-BEN.inp. 

Run by entering "python BenPre.py" in the command line.

The user will be prompted to enter the analysis type, solution algorithm type, solver algorithm type, and to choose to execute node-renumbering in the command line. Options for each of the analysis (etc) types will be provided and the user is to enter the numerical flag associated with their choice. The code then will access the Abaqus input file to transcribe the geometry established in Abaqus to be read by CU-BEN. Depending on the choice of analysis/algorithm, the user will be prompted to enter additional parameters for the solution. If the prompt requires more than one value, a guide will be printed below the prompt to help the user input the values correctly. 

A tutorial on how to correctly set up the Abaqus file to be read seemlessly by BenPre can be found in the "Introduction to CU-BEN" document. The primary concerns to consider when constructing an Abaqus model are to use only element type S3 shell elements and to title the .inp file "CU-BEN.inp". If running a dynamic analysis, it is important to define an explicit time step and enter the load pattern in a tabular form. The Step created for the time steps must be called "TimeSteps". If the load on any nodes in not zero at time t=0 or is there are any initial perscribed displacement, velocity, or accerlation on any of the nodes, this must be defined in a Step called "Time0". More details can be found in "Introduction to CU-BEN". 

The code will output a file named "model_def.txt". This can be used as an input file in CU-BEN. 