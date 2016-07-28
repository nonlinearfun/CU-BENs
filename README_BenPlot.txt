README_BenPlot: BenPlot.py is a rudimentary plotting code for CU-BEN. This code was written in Python 2.7.

Usage:

Keep in same folder as CU-BEN's results1.txt and results2.txt files. 

Run by entering "python BenPlot.py" in the command line.

The user will be prompted to input the xyz-coordinates of the joint of interest. If the combination of coordinate values is not found, the code will output an error and terminate. If the coordinates do not have active degrees of freedom, the code output an error and terminate. 

For each DOF associated with the specified joint, the code will create a plot. If the analysis is dynamic, the load response will be plotted against the time. If the analysis is static and nonlinear, the the load proportionality factor will be plotted against the load response. The plots will be output based on their numbering. If the user would like to save the plot, save the figure to the user's location of choice. When the user closes the window containing the plot, the next DOF on the joint will be plotted. 