#Beta version, Justyna Kosianka 1/04/18

import csv
import matplotlib.pyplot as plt

#Input coordinates
x_coord,y_coord,z_coord = map(float, raw_input("Enter Joint Coordinates to Plot [x,y,z]: ").split(','))

#Search results1.txt for joint number associated with coordinates
results1 = list(csv.reader(open('results1.txt', 'rb'), delimiter='\t'))

#Convert
x_coord = str('%.6f' % x_coord)
y_coord = str('%.6f' % y_coord)
z_coord = str('%.6f' % z_coord)

#Find index of x_coords
x_ind = [(index, row.index(x_coord)) for index, row in enumerate(results1) if x_coord in row]

x_col = [x[1] for x in x_ind]
x_row = [x[0] for x in x_ind]
indices = [i for i, x in enumerate(x_col) if x == 3]
x_lines = [x_row[i] for i in indices]
for i in x_lines:
    results1[i][3] = "None"

#Then scan for y_coords
y_ind = [(index, row.index(y_coord)) for index, row in enumerate(results1) if y_coord in row]

y_col = [y[1] for y in y_ind]
y_row = [y[0] for y in y_ind]
indices = [i for i, y in enumerate(y_col) if y == 4]
y_lines = [y_row[i] for i in indices]
for i in y_lines:
    results1[i][4] = "None"

#Then scan those indeces for z_coords to find final line index
z_ind = [(index, row.index(z_coord)) for index, row in enumerate(results1) if z_coord in row]

z_col = [z[1] for z in z_ind]
z_row = [z[0] for z in z_ind]
indices = [i for i, z in enumerate(z_col) if z == 5]
z_lines = [z_row[i] for i in indices]

#If coordinates have no joint number, output error and terminate
if not x_lines:
    print "Coordinates invalid"
    exit()
if not y_lines:
    print "Coordinates invalid"
    exit()
if not z_lines:
    print "Coordinates invalid"
    exit()

#Then get joint number
joint_line = list(set(x_lines) & set(y_lines) & set(z_lines))
if not joint_line:
    print "Coordinates invalid"
    exit()
joint = results1[joint_line[0]][1]

#Scan to find indices that match joint number and choose the index that is one beyond the one established before
joint_ind = [(index, row.index(joint)) for index, row in enumerate(results1) if joint in row]

j_col = [j[1] for j in joint_ind]
j_row = [j[0] for j in joint_ind]
indices = [i for i, j in enumerate(j_col) if j == 1]
j_line = [j_row[i] for i in indices]
indices = [i for i in j_line if i > joint_line[0]]
j_line = indices[0]

#Go to that row and collect DOFs, track the label
DOF = [0]*7
DOF_label = [0]*7

DOF[0] = int(results1[j_line][3])
if DOF[0] != 0:
    DOF_label[0] = "DOF %d: X-Translation" % DOF[0]

DOF[1] = int(results1[j_line][5])
if DOF[1] != 0:
    DOF_label[1] = "DOF %d: Y-Translation" % DOF[1]

DOF[2] = int(results1[j_line][7])
if DOF[2] != 0:
    DOF_label[2] = "DOF %d: Z-Translation" % DOF[2]

DOF[3] = int(results1[j_line][9])
if DOF[3] != 0:
    DOF_label[3] = "DOF %d: X-Rotation" % DOF[3]

DOF[4] = int(results1[j_line][11])
if DOF[4] != 0:
    DOF_label[4] = "DOF %d: Y-Rotation" % DOF[4]

DOF[5] = int(results1[j_line][13])
if DOF[5] != 0:
    DOF_label[5] = "DOF %d: Z-Rotation" % DOF[5]

DOF[6] = int(results1[j_line][15])
if DOF[6] != 0:
    DOF_label[6] = "DOF %d: Warping" % DOF[6]

#If joint does not have DOFs, output error and terminate
if not [i for i in DOF if i > 0]:
    print "Chosen Coordinate does not have active DOFs"
    exit()

#Go to results2 and track the column associated with DOF
results2 = list(csv.reader(open('results2.txt', 'rb'), delimiter='\t'))
results2 = results2[2:(len(results2)-1)][:]

constant_axis = map(float,[d[1] for d in results2])

#Plot seperate figure for each DOF
fig_count = 1

for i in range(0,6):
    if DOF[i] != 0:
        plt.figure(fig_count)
        fig_count = fig_count+1
        DOF_plot = map(float,[d[DOF[i]+3] for d in results2])
        #If dynamic, plot time on x and displacement on y
        if ((results1[4][1] == 'Dynamic (Newmark)') or (results1[4][1] == 'Dynamic (Newmark with geometric nonlinearity)')):
            plt.plot(constant_axis,DOF_plot)
            plt.xlabel('Time (sec)')
            if 0 <= i <= 2:
                plt.ylabel('Displacement')
            elif 3 <= i <= 5:
                plt.ylabel('Rotation')
            elif i == 6:
                plt.ylabel('Warping')
            plt.title(DOF_label[i])
        #If nonlinear, plot LPF on y and displacement on x
        elif ((results1[4][1] == 'Newton Raphson') or (results1[4][1] == 'Modified Newton Raphson') or (results1[4][1] == 'Modified Spherical Arc Length')):
            plt.plot(DOF_plot,constant_axis)
            plt.ylabel('Load Proportionality Factor (LPF)')
            if 0 <= i <= 2:
                plt.xlabel('Displacement')
            elif 3 <= i <= 5:
                plt.xlabel('Rotation')
            elif i == 6:
                plt.xlabel('Warping')
            plt.title(DOF_label[i])
        else:
            print "This type of analysis does not produce plot."
            exit()
        plt.show()
