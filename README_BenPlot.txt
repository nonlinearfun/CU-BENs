README_BenPlot: BenPlot.py is a rudimentary plotting routine for CU-BEN. BenPlot.py is ported to Python 2.7.

Usage:

Keep in same folder as CU-BEN's results1.txt and results2.txt files. 

Run by entering "python BenPlot.py" in the command line.

The user will be prompted to input the xyz-coordinates pertaining to the joint of interest. If the combination of coordinate values is not found, the routine will output an error and terminate. If the coordinates do not have active degrees of freedom (DOFs), the routine will output an error and terminate. 

For each DOF associated with the specified joint, the routine will create a plot. If the analysis is dynamic, the load response will be plotted against the solution time. If the analysis is static and nonlinear, the load proportionality factor will be plotted against the load response. The plots are outputted based on their ordering. If the user would like to save the plot, they will be prompted about the path for the save. When the user closes the window containing the plot, the next DOF associated with the given joint will be plotted. 