#BenPre - a automated assistant for building CU-BEN input files from Abaqus *.inp files

#Beta version, Justyna Kosianka, 1/04/18

import csv

#User passes through file name for Abaqus input
file_name = raw_input("Enter Abaqus *.inp File Name:\n>> ")

#Collect geometry information from Abaqus input file
#Input CU-BEN.inp
CUBEN = list(csv.reader(open(file_name, 'rb'), delimiter=','))

#Determine number of joints from .inp file
var = "*Node"
line_node = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
if not line_node:
    var = "*NODE"
    line_node = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
if not line_node:
    var = "*Node                                                                                                                   "
    line_node = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
line_node = line_node[0][0]
var = "*Element"
line_elem = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
if not line_elem:
    var = var = "*ELEMENT"
    line_elem = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
line_elem = line_elem[0][0]
NJ = (line_elem - line_node - 1)

#Determine number of elements from .inp file for SHELL ONLY
var = "*Nset"
line_Nset = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
if not line_Nset:
    var = "*NSET"
    line_Nset = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
if not line_Nset:
    print "Model was not built in accordance to file creation protocol. This will not successfully build a CU-BEN input file."
    exit()
line_Nset1 = line_Nset[0][0]
line_Nset = [y[0] for y in line_Nset]
var = "*Elset"
line_Elset = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
if not line_Elset:
    var = "*ELSET"
    line_Elset = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
line_Elset1 = line_Elset[0][0]
line_Elset = [y[0] for y in line_Elset]

if line_Nset1 < line_Elset1:
    NSH = (line_Nset1 - line_elem - 1)
else:
    NSH = (line_Elset1 - line_elem - 1)

#Collect SHELL member incidence from *Element, type=S3 in .inp file
minc = CUBEN[(line_elem + 1):(line_elem + NSH + 1)]

#Check if model is meshed correctly to build model_def.txt
if len(minc[0]) > 4:
    print "Model was not meshed with S3 shell elements in Abaqus. Please re-mesh and try again."
    exit()

var = "** STEP: TimeSteps"
line_TimeSteps = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
if not line_TimeSteps:
    var = "** STEP: TimeStep"
    line_TimeSteps = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
if not line_TimeSteps:
    var = "** STEP: Step-1"
    line_TimeSteps = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
if not line_TimeSteps:
    print "Model was not loaded in accordance to file creation protocol. Please re-apply load and try again."
    exit()
line_TimeSteps = line_TimeSteps[0][0]
var = "*End Step"
line_EndStep = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
line_EndStep = [y[0] for y in line_EndStep]
end_TimeSteps = [j for j in line_EndStep if j > line_TimeSteps]
end_TimeSteps = end_TimeSteps[0]

#Create a file called model_def.txt
model_def = open('model_def.txt','w')

#User inputs analysis type
ANAFLAG = input("Enter analysis type:\n[1 - 1st order elastic / 2 - 2nd order elastic / 3 - 2nd order inelastic]\n>> ")

#User inputs algorithm type
ALGFLAG = input("Enter solution algorithm type:\n[0 - (Static) Direct Stiffness / 1 - (Static) Newton Raphson / 2 - (Static) Modified Newton Raphson / 3 - (Static) Modified Spherical Arc Length Method / 4 - (Dynamic) Newmark Implicit Integration Method / 5 - (Dynamic) Nonlinear Newmark Implicit Integration Method]\n>> ")

#Check for invalid combination of analysis and algorithm
if ((ALGFLAG == 0) and (ANAFLAG > 1)):
    print "Invalid combination of analysis and solution algorithm type: Static Direct Stiffness method can only analyze 1st order elatic models."
    ANAFLAG = input("Enter analysis type:\n[1 - 1st order elastic / 2 - 2nd order elastic / 3 - 2nd order inelastic]\n>> ")
    ALGFLAG = input("Enter solution algorithm type:\n[0 - (Static) Direct Stiffness / 1 - (Static) Newton Raphson / 2 - (Static) Modified Newton Raphson / 3 - (Static) Modified Spherical Arc Length Method / 4 - (Dynamic) Newmark Implicit Integration Method / 5 - (Dynamic) Nonlinear Newmark Implicit Integration Method]\n>> ")
elif ((ALGFLAG == 4) and (ANAFLAG > 1)):
    print "Invalid combination of analysis and solution algorithm type: Dynamic Newmark Implicit Integration Method can only analyze 1st order linear elastic models."
    ANAFLAG = input("Enter analysis type:\n[1 - 1st order elastic / 2 - 2nd order elastic / 3 - 2nd order inelastic]\n>> ")
    ALGFLAG = input("Enter solution algorithm type:\n[0 - (Static) Direct Stiffness / 1 - (Static) Newton Raphson / 2 - (Static) Modified Newton Raphson / 3 - (Static) Modified Spherical Arc Length Method / 4 - (Dynamic) Newmark Implicit Integration Method / 5 - (Dynamic) Nonlinear Newmark Implicit Integration Method]\n>> ")
elif ((ALGFLAG == 5) and (ANAFLAG == 1)):
    print "Invalid combination of analysis and solution algorithm type: Dynamic Nonlinear Newmark Implicit Integration Method cannot analyze 1st order linear elastic models."
    ANAFLAG = input("Enter analysis type:\n[1 - 1st order elastic / 2 - 2nd order elastic / 3 - 2nd order inelastic]\n>> ")
    ALGFLAG = input("Enter solution algorithm type:\n[0 - (Static) Direct Stiffness / 1 - (Static) Newton Raphson / 2 - (Static) Modified Newton Raphson / 3 - (Static) Modified Spherical Arc Length Method / 4 - (Dynamic) Newmark Implicit Integration Method / 5 - (Dynamic) Nonlinear Newmark Implicit Integration Method]\n>> ")
elif ((ALGFLAG == 1) and (ANAFLAG == 1)):
    print "Invalid combination of analysis and solution algorithm type: Static Newton Raphson Method cannot analyze 1st order linear elastic models."
    ANAFLAG = input("Enter analysis type:\n[1 - 1st order elastic / 2 - 2nd order elastic / 3 - 2nd order inelastic]\n>> ")
    ALGFLAG = input("Enter solution algorithm type:\n[0 - (Static) Direct Stiffness / 1 - (Static) Newton Raphson / 2 - (Static) Modified Newton Raphson / 3 - (Static) Modified Spherical Arc Length Method / 4 - (Dynamic) Newmark Implicit Integration Method / 5 - (Dynamic) Nonlinear Newmark Implicit Integration Method]\n>> ")
elif ((ALGFLAG == 2) and (ANAFLAG == 1)):
    print "Invalid combination of analysis and solution algorithm type: Static Modified Newton Raphson Method cannot analyze 1st order linear elastic models."
    ANAFLAG = input("Enter analysis type:\n[1 - 1st order elastic / 2 - 2nd order elastic / 3 - 2nd order inelastic]\n>> ")
    ALGFLAG = input("Enter solution algorithm type:\n[0 - (Static) Direct Stiffness / 1 - (Static) Newton Raphson / 2 - (Static) Modified Newton Raphson / 3 - (Static) Modified Spherical Arc Length Method / 4 - (Dynamic) Newmark Implicit Integration Method / 5 - (Dynamic) Nonlinear Newmark Implicit Integration Method]\n>> ")
elif ((ALGFLAG == 3) and (ANAFLAG == 1)):
    print "Invalid combination of analysis and solution algorithm type: Static Modified Spherical Arc Length Method cannot analyze 1st order linear elastic models."
    ANAFLAG = input("Enter analysis type:\n[1 - 1st order elastic / 2 - 2nd order elastic / 3 - 2nd order inelastic]\n>> ")
    ALGFLAG = input("Enter solution algorithm type:\n[0 - (Static) Direct Stiffness / 1 - (Static) Newton Raphson / 2 - (Static) Modified Newton Raphson / 3 - (Static) Modified Spherical Arc Length Method / 4 - (Dynamic) Newmark Implicit Integration Method / 5 - (Dynamic) Nonlinear Newmark Implicit Integration Method]\n>> ")

model_def.write('%d\n' % ANAFLAG)
model_def.write('%d\n' % ALGFLAG)

#User inputs solver algorithm type
SLVFLAG = input("Enter solver algorithm type:\n[0 - CU-BEN for symmetric matrices / 1 - CLAPACK solver for symmetric and non-symmetric matrices]\n>> ")
model_def.write('%d\n' % SLVFLAG)

#User inputs node-numbering algorithm
OPTFLAG = input("Enter flag for execution of node-renumbering algorithm:\n[1 - No / 2 - Yes]\n>> ")
model_def.write('%d\n' % OPTFLAG)

#Write number of joints to model_def.txt
model_def.write('%d\n' % NJ)

#Write nuber of shell elements to model_def.txt
model_def.write('0,0,%d,0,0\n' % NSH)

#Write member incidence to
for i in range(0, NSH):
    model_def.write('%d,%d,%d\n' % (int(minc[i][1]), int(minc[i][2]), int(minc[i][3])))

#Collect joint constraints from BC-set and *Boundary in .inp file
var = "*Boundary"
line_Boundary = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
line_Boundary = [y[0] for y in line_Boundary]
line_Boundary = [j for j in line_Boundary if ((j < end_TimeSteps) and (j > line_TimeSteps))]
num_Boundary = len(line_Boundary[:])

for i in range(0,num_Boundary):
    BC_name = CUBEN[line_Boundary[i]+1][0]
    var = " nset=" + BC_name
    line_BCset = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
    line_BCset = line_BCset[0][0]
    var = "*Step"
    line_Step = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
    line_Step = [y[0] for y in line_Step]
    var = "*End Assembly"
    line_Assembly = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
    line_Assembly = [y[0] for y in line_Assembly]

    test = [j for j in line_Nset if j > line_BCset]
    test = test + [j for j in line_Elset if j > line_BCset]
    test = test+ [j for j in line_Step if j > line_BCset]
    test = test+ [j for j in line_Assembly if j > line_BCset]
    end_BCset = min(test)

    BCset_length = end_BCset - line_BCset - 1

    var = "** "
    line_star = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
    line_star = [y[0] for y in line_star]
    test = [j for j in line_star if j > line_Boundary[i]]
    name_test = [j for j in line_Boundary if j > line_Boundary[i]]
    for j in range(0,len(name_test)):
        if CUBEN[name_test[j]-1][0] == BC_name:
            test = test + [name_test[j]]
        elif CUBEN[name_test[j]-2][0] == BC_name:
            test = test + [name_test[j]-1]
    end_Boundary = min(test)
    Boundary_length = end_Boundary - line_Boundary[i] - 1

    BCset = CUBEN[(line_BCset + 1):(end_BCset)]
    
    BCset_array = None

    for j in range(0,BCset_length):
        for k in range(1,len(BCset[:][j])):
            if BCset[j][k] < BCset[j][k-1]:
                BCset_array = range(int(BCset[0][0]),int(BCset[j][k-1])+1,int(BCset[j][k]))

    if BCset_array:
        for j in range(0,len(BCset_array)):
            for l in range(0,Boundary_length):
                model_def.write('%d,%d\n' % (int(BCset_array[j]), int(CUBEN[(line_Boundary[i]+1+l)][1])))
    else:
        for j in range(0,BCset_length):
            for k in range(0,len(BCset[:][j])):
                for l in range(0,Boundary_length):
                    model_def.write('%d,%d\n' % (int(BCset[j][k]), int(CUBEN[(line_Boundary[i]+1+l)][1])))

    BCset_array = None

#Terminate list in model_def.txt using 0,0
model_def.write('0,0\n')

#Collect joint coordinates from *Node
coords = CUBEN[(line_node + 1):(line_elem)]
for i in range(0, NJ):
    model_def.write('%g,%g,%g\n' % (float(coords[i][1]), float(coords[i][2]), float(coords[i][3])))


#Collect shell element properties (elastic modulus, Poisson's ratio, thickness, density, maximum allowable yield stress)
#Go through Sections and apply shell section and material to elements defined in elset
var = "*Shell Section"
line_Sections = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
line_Sections = [y[0] for y in line_Sections]
num_Sections = len(line_Sections[:])

elemProp = [[0 for x in range(5)] for y in range(NSH)]

for i in range(0,num_Sections):
    var = CUBEN[line_Sections[i]][2]
    var = " name=" + var[10:]
    line_Material = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
    line_Material = line_Material[0][0]
    var = CUBEN[line_Sections[i]][1]
    if (var == " elset=Eall"):
        for j in range(0,NSH):
            #Elastic modulus - first entry under *Elastic in *Material
            elemProp[j][0] = float(CUBEN[line_Material+4][0])
            #Poisson's ratio - second entry under *Elastic in *Material
            elemProp[j][1] = float(CUBEN[line_Material+4][1])
            #Thickness - defined in *Shell Section, first entry
            elemProp[j][2] = float(CUBEN[line_Sections[i]+1][0])
            #Density - defined under *Density in *Material
            elemProp[j][3] = float(CUBEN[line_Material+2][0])
            #Maximum allowable yield stress - defined under *Plastic in Material, first entry
            elemProp[j][4] = float(CUBEN[line_Material+6][0])
    else:
        elset_Section = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
        elset_Section = elset_Section[0][0]+1
        if (len(CUBEN[elset_Section-1]) > 2):
            if (CUBEN[elset_Section-1][2] == " generate"):
                start = int(CUBEN[elset_Section][0])
                end = int(CUBEN[elset_Section][1])
                increment = int(CUBEN[elset_Section][2])
                for j in range(start-1,end,increment):
                    #Elastic modulus - first entry under *Elastic in *Material
                    elemProp[j][0] = float(CUBEN[line_Material+4][0])
                    #Poisson's ratio - second entry under *Elastic in *Material
                    elemProp[j][1] = float(CUBEN[line_Material+4][1])
                    #Thickness - defined in *Shell Section, first entry
                    elemProp[j][2] = float(CUBEN[line_Sections[i]+1][0])
                    #Density - defined under *Density in *Material
                    elemProp[j][3] = float(CUBEN[line_Material+2][0])
                    #Maximum allowable yield stress - defined under *Plastic in Material, first entry
                    elemProp[j][4] = float(CUBEN[line_Material+6][0])
        else:
            test = [j for j in line_Nset if j > elset_Section]
            test = test + [j for j in line_Elset if j > elset_Section]
            test = test + [j for j in line_Assembly if j > elset_Section]
            name_test = [j for j in line_Sections if j > elset_Section]
            for j in range(0,len(name_test)):
                name_test[j] = name_test[j] -1
            test = test + name_test
            end_elset = min(test)
            elset = CUBEN[(elset_Section):(end_elset)]
            for j in range(0,len(elset)):
                for k in range(0,len(elset[:][j])):
                    l = int(elset[j][k])-1
                    #Elastic modulus - first entry under *Elastic in *Material
                    elemProp[l][0] = float(CUBEN[line_Material+4][0])
                    #Poisson's ratio - second entry under *Elastic in *Material
                    elemProp[l][1] = float(CUBEN[line_Material+4][1])
                    #Thickness - defined in *Shell Section, first entry
                    elemProp[l][2] = float(CUBEN[line_Sections[i]+1][0])
                    #Density - defined under *Density in *Material
                    elemProp[l][3] = float(CUBEN[line_Material+2][0])
                    #Maximum allowable yield stress - defined under *Plastic in Material, first entry
                    elemProp[l][4] = float(CUBEN[line_Material+6][0])

for i in range(0, NSH):
    model_def.write('%g,%g,%g,%g,%g\n' % (float(elemProp[i][0]), float(elemProp[i][1]), float(elemProp[i][2]), float(elemProp[i][3]), float(elemProp[i][4])))

#Start collecting algorithm specific input variables
if (ALGFLAG < 4): #Static
    
    #Collect load (joint, direction, force) from .inp file
    var = "*Cload"
    line_Load = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
    line_Load = [y[0] for y in line_Load]
    num_Load = len(line_Load[:])

    for i in range(0,num_Load):
        Load_name = CUBEN[line_Load[i]+1][0]
        var = " nset=" + Load_name
        line_Loadset = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
        line_Loadset = line_Loadset[0][0]
    
        test = [j for j in line_Nset if j > line_Loadset]
        test = test + [j for j in line_Elset if j > line_Loadset]
        test = test+ [j for j in line_Assembly if j > line_Loadset]
        end_Loadset = min(test)
    
        Loadset_length = end_Loadset - line_Loadset - 1
    
        test = [j for j in line_star if j > line_Load[i]]
        name_test = [j for j in line_Load if j > line_Load[i]]
        for j in range(0,len(name_test)):
            if CUBEN[name_test[j]-1][0] == Load_name:
                test = test + [name_test[j]]
            elif CUBEN[name_test[j]-2][0] == Load_name:
                test = test + [name_test[j]-1]
        end_Load = min(test)
        Load_length = end_Load - line_Load[i] - 1
    
        Loadset = CUBEN[(line_Loadset + 1):(end_Loadset)]
    
        for j in range(0,Loadset_length):
            for k in range(0,len(Loadset[:][j])):
                for l in range(0,(Load_length)):
                    model_def.write('%d,%d,%f\n' % (int(Loadset[j][k]), int(CUBEN[(line_Load[i]+1+l)][1]), float(CUBEN[(line_Load[i]+1+l)][2])))
    #End of loads
    model_def.write('0,0,0\n')

    #User inputs algorithm parameters
    if (ALGFLAG == 0) or ((ALGFLAG == 1) and (ANAFLAG == 1)): #(Direct Stiffness)
        model_def.write('1\n')
    elif ALGFLAG == 1: #(Newton Raphson)
        LPFmax = input("Enter maximum load proportionality factor:\n>> ")
        model_def.write('%f\n' % LPFmax)
    elif ALGFLAG == 2: #(Modified Newton Raphson)
        print "Enter load proportionality factor parameters:"
        LPFmax = input("Maximum lambda:\n>> ")
        LPF = input("Initial lambda:\n>> ")
        dLPF = input("Increment of lambda:\n>> ")
        dLPFmax = input("Maximum increment of lambda:\n>> ")
        dLPFmin = input("Minimum increment of lambda:\n>> ")
        model_def.write('%g,%g,%g,%g,%g\n' % (LPFmax, LPF, dLPF, dLPFmax, dLPFmin))
        print "Enter maximums and minimums on counters:"
        itemax = input("Maximum number of iterations within load increment:\n>> ")
        submax = input("Maximum number of times to step back load due to unconverged solution:\n>> ")
        solmin = input("Minimum number of converged solutions before increasing increment of lambda:\n>> ")
        model_def.write('%d,%d,%d\n' % (itemax, submax, solmin))
    elif ALGFLAG == 3: #(Modified Spherical Arc Length Method)
        print "Enter MSAL parameters:"
        jnum,jdir,dk = map(float, raw_input("Initial prescribed displacement at DOF k:\n[displacement, direction, DOF k]\n>> ").split(','))
        model_def.write('%g,%g,%g\n' % (jnum,jdir,dk))
        alpha = input("Factor limiting size of load increment:\n>> ")
        model_def.write('%g\n' % alpha)
        psi_thresh = input("Threshold to eliminate load control in the arc length criterion in the vicinity of a critical point (A minimum value of 0.1 is recommended):\n>> ")
        model_def.write('%g\n' % psi_thresh)
        iteopt = input("Optimum number of iterations (a value of 6 is recommended):\n>> ")
        model_def.write('%g\n' % iteopt)
        lpfmax,dkimax = map(float, raw_input("Maximum lambda and allowable displacment:\n[maximum lambda, allowable displacement]\n>> "))
        model_def.write('%g,%g\n' % (lpfmax, dkimax))
        print "Enter maximums / minimums on counters:"
        itemax = input("Maximum number of iterations within load increment:\n>> ")
        submax = input("Maximum number of times to step back load due to unconverged solution:\n>> ")
        imagmax = input("Maximum number of times to step back load due to arc length criterion producing imaginary roots:\n>> ")
        negmax = input("maximum number of times to step back load due to arc length criterion producing two negative roots:\n>> ")
        model_def.write('%g,%g,%g,%g\n' % (itemax, submax, imagmax, negmax))

else: #Dynamic

    #Enter initial number of time steps (in load) and total time for analysis (s);
    var = "*Dynamic"
    line_Dynamic = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
    line_Dynamic = [y[0] for y in line_Dynamic]
    line_Dynamic = [j for j in line_Dynamic if ((j < end_TimeSteps) and (j > line_TimeSteps))]
    line_Dynamic = line_Dynamic[0] + 1
    ttot = float(CUBEN[line_Dynamic][1])
    ntstpsinpt = int(round(ttot/float(CUBEN[line_Dynamic][0])))
    model_def.write('%d,%g\n' % (ntstpsinpt, ttot))
            
    #enter reference concentrated load(s) on joints for each time step (in load); i = 0:ntstpsinpt
    #Collect load (joint, direction, force) from .inp file
    var = "*Cload"
    line_Load = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
    line_Load = [y[0] for y in line_Load]
    line_Load = [j for j in line_Load if ((j < end_TimeSteps) and (j > line_TimeSteps))]
    num_Load = len(line_Load[:])
    
    #Determine if Abaqus File contains load history
    var = "*Amplitude"
    line_Amp = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
    #Determine if Abaqus File contains non-zero initial conditions
    var = "** STEP: Time0"
    line_Time0 = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
    if line_Amp:
        line_Amp = [y[0] for y in line_Amp]
        num_Amp = len(line_Amp)
        if (num_Amp == 1):
            line_Amp = line_Amp[0]
            test = [j for j in line_star if j > line_Amp]
            end_Amp = min(test)
            Amp_line_length = end_Amp - line_Amp - 1
        else:
            for i in range(0,num_Amp):
                test = [j for j in line_star if j > line_Amp[i]]
                test = test + [j for j in line_Amp if j > line_Amp[i]]
                end_Amp = min(test)
                Amp_line_length = end_Amp - line_Amp[i] - 1
                if (Amp_line_length > 1):
                    line_Amp = line_Amp[i]
                    break
    
        Load_pattern = [0 for x in range(ntstpsinpt+1)]
    
        for i in range(0,Amp_line_length):
            for j in range(0, len(CUBEN[:][line_Amp+1+i]),2):
                k = int(round(float(CUBEN[line_Amp+1+i][j])/(ttot/ntstpsinpt)))
                Load_pattern[k] = float(CUBEN[line_Amp+1+i][j+1])

        for i in range(0,num_Load):
            Load_name = CUBEN[line_Load[i]+1][0]
            var = " nset=" + Load_name
            line_Loadset = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
            line_Loadset = line_Loadset[0][0]
        
            test = [j for j in line_Nset if j > line_Loadset]
            test = test + [j for j in line_Elset if j > line_Loadset]
            test = test+ [j for j in line_Assembly if j > line_Loadset]
            end_Loadset = min(test)
        
            Loadset_length = end_Loadset - line_Loadset - 1
        
            test = [j for j in line_star if j > line_Load[i]]
            name_test = [j for j in line_Load if j > line_Load[i]]
            for j in range(0,len(name_test)):
                if CUBEN[name_test[j]-1][0] == Load_name:
                    test = test + [name_test[j]]
                elif CUBEN[name_test[j]-2][0] == Load_name:
                    test = test + [name_test[j]-1]
            end_Load = min(test)
            Load_length = end_Load - line_Load[i] - 1
        
            Loadset = CUBEN[(line_Loadset + 1):(end_Loadset)]
        
            #Collect initial values for load if not equal to 0 at time t=0
            if line_Time0:
                line_Time0 = line_Time0[0][0]
                end_Time0 = [j for j in line_EndStep if j > line_Time0]
                end_Time0 = end_Time0[0]
                #Collect load sets, directions, and magnitude for initial loads
                var = "*Cload"
                line_LoadInit = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
                line_LoadInit = [y[0] for y in line_LoadInit]
                line_LoadInit = [j for j in line_LoadInit if ((j < end_Time0) and (j > line_Time0))]
                num_LoadInit = len(line_LoadInit[:])

                for n in range(0,num_LoadInit):
                    LoadInit_name = CUBEN[line_LoadInit[n]+1][0]
                    var = " nset=" + LoadInit_name
                    line_LoadInitset = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
                    line_LoadInitset = line_LoadInitset[0][0]
        
                    test = [j for j in line_Nset if j > line_LoadInitset]
                    test = test + [j for j in line_Elset if j > line_LoadInitset]
                    test = test+ [j for j in line_Assembly if j > line_LoadInitset]
                    end_LoadInitset = min(test)
        
                    LoadInitset_length = end_LoadInitset - line_LoadInitset - 1
        
                    test = [j for j in line_star if j > line_LoadInit[n]]
                    name_test = [j for j in line_Load if j > line_LoadInit[n]]
                    for k in range(0,len(name_test)):
                        if CUBEN[name_test[k]-1][0] == Load_name:
                            test = test + [name_test[j]]
                        elif CUBEN[name_test[k]-2][0] == Load_name:
                            test = test + [name_test[j]-1]
                    end_LoadInit = min(test)
                    LoadInit_length = end_LoadInit - line_LoadInit[n] - 1
        
                    LoadInitset = CUBEN[(line_LoadInitset + 1):(end_LoadInitset)]
                
                    LoadInit = [[0 for x in range(3)] for y in range((LoadInitset_length)*(len(LoadInitset[:][0])))]
                    d = 0
                    for a in range(0,LoadInitset_length):
                        for b in range(0,len(LoadInitset[:][a])):
                            for c in range(0,LoadInit_length):
                                LoadInit[d][0] = LoadInitset[a][b]
                                LoadInit[d][1] = CUBEN[(line_LoadInit[n]+1+c)][1]
                                LoadInit[d][2] = CUBEN[(line_LoadInit[n]+1+c)][2]
                                d = d + 1
                    if n == 0:
                        LoadInit_Total = LoadInit
                    else:
                        LoadInit_Total = LoadInit_Total + LoadInit

            for j in range(0,Loadset_length):
                for k in range(0,len(Loadset[:][j])):
                    for l in range(0,Load_length):
                        for m in range(0,ntstpsinpt+1):
                            if (m == 0):
                                if not line_Time0:
                                    model_def.write('%d,%d,%g\n' % (int(Loadset[j][k]), int(CUBEN[(line_Load[i]+1+l)][1]), float(CUBEN[(line_Load[i]+1+l)][2])*Load_pattern[m]))
                                else:
                                    if Loadset[j][k] in [y[0] for y in LoadInit_Total]:
                                        LoadInit_ind = [(index, row.index(Loadset[j][k])) for index, row in enumerate([y[0] for y in LoadInit_Total]) if Loadset[j][k] in row]
                                        if len(LoadInit_ind) > 1:
                                            for n in range(0,len(LoadInit_ind)):
                                                model_def.write('%d,%d,%g\n' % (int(Loadset[j][k]), int(LoadInit_Total[(LoadInit_ind[n][0])][1]), float(LoadInit_Total[(LoadInit_ind[n][0])][2])))
                                        else:
                                            model_def.write('%d,%d,%g\n' % (int(Loadset[j][k]), int(LoadInit_Total[(LoadInit_ind[0][0])][1]), float(LoadInit_Total[(LoadInit_ind[0][0])][2])))
                                    else:
                                        model_def.write('%d,%d,%g\n' % (int(Loadset[j][k]), int(CUBEN[(line_Load[i]+1+l)][1]), float(CUBEN[(line_Load[i]+1+l)][2])*Load_pattern[m]))
                            else:
                                model_def.write('%d,%d,%g\n' % (int(Loadset[j][k]), int(CUBEN[(line_Load[i]+1+l)][1]), float(CUBEN[(line_Load[i]+1+l)][2])*Load_pattern[m]))
    
    #Terminate list in model_def.txt using 0,0,0
    model_def.write('0,0,0\n')
    
    #enter all non-zero initial conditions for node displacement, velocity, and acceleration (in load)
    #joint,dir,disp,vel,acc;
    if line_Time0:
        var = "*Boundary"
        line_Boundary = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
        line_Boundary = [y[0] for y in line_Boundary]
        line_Boundary = [j for j in line_Boundary if ((j < end_Time0) and (j > line_Time0))]
        num_Boundary = len(line_Boundary[:])

        for i in range(0,num_Boundary):
            BC_name = CUBEN[line_Boundary[i]+1][0]
            var = " nset=" + BC_name
            line_BCset = [(index, row.index(var)) for index, row in enumerate(CUBEN) if var in row]
            line_BCset = line_BCset[0][0]
    
            test = [j for j in line_Nset if j > line_BCset]
            test = test + [j for j in line_Elset if j > line_BCset]
            test = test+ [j for j in line_Step if j > line_BCset]
            test = test+ [j for j in line_Assembly if j > line_BCset]
            end_BCset = min(test)
    
            BCset_length = end_BCset - line_BCset - 1
    
            test = [j for j in line_star if j > line_Boundary[i]]
            name_test = [j for j in line_Boundary if j > line_Boundary[i]]
            for j in range(0,len(name_test)):
                if CUBEN[name_test[j]-1][0] == BC_name:
                    test = test + [name_test[j]]
                elif CUBEN[name_test[j]-2][0] == BC_name:
                    test = test + [name_test[j]-1]
            end_Boundary = min(test)
            Boundary_length = end_Boundary - line_Boundary[i] - 1
    
            BCset = CUBEN[(line_BCset + 1):(end_BCset)]
    
            for j in range(0,BCset_length):
                for k in range(0,len(BCset[:][j])):
                    for l in range(0,(Boundary_length)):
                        if (len(CUBEN[:][line_Boundary[i]+1+l]) == 4):
                            if not ('BoundaryInit' in locals() or 'BoundaryInit' in globals()):
                                if (len(CUBEN[:][line_Boundary[i]]) == 2):
                                    BoundaryInit = [int(BCset[j][k]), int(CUBEN[(line_Boundary[i]+1+l)][1]), float(CUBEN[(line_Boundary[i]+1+l)][3]), 0, 0]
                                elif ((len(CUBEN[:][line_Boundary[i]]) == 3) and (CUBEN[line_Boundary[i]][2] == " type=VELOCITY")):
                                    BoundaryInit = [int(BCset[j][k]), int(CUBEN[(line_Boundary[i]+1+l)][1]), 0, float(CUBEN[(line_Boundary[i]+1+l)][3]), 0]
                                elif ((len(CUBEN[:][line_Boundary[i]]) == 3) and (CUBEN[line_Boundary[i]][2] == " type=ACCELERATION")):
                                    BoundaryInit = [int(BCset[j][k]), int(CUBEN[(line_Boundary[i]+1+l)][1]), 0, 0, float(CUBEN[(line_Boundary[i]+1+l)][3])]
                            else:
                                if (len(CUBEN[:][line_Boundary[i]]) == 2):
                                    BoundaryInit = BoundaryInit + [int(BCset[j][k]), int(CUBEN[(line_Boundary[i]+1+l)][1]), float(CUBEN[(line_Boundary[i]+1+l)][3]), 0, 0]
                                elif ((len(CUBEN[:][line_Boundary[i]]) == 3) and (CUBEN[line_Boundary[i]][2] == " type=VELOCITY")):
                                    BoundaryInit = BoundaryInit + [int(BCset[j][k]), int(CUBEN[(line_Boundary[i]+1+l)][1]), 0, float(CUBEN[(line_Boundary[i]+1+l)][3]), 0]
                                elif ((len(CUBEN[:][line_Boundary[i]]) == 3) and (CUBEN[line_Boundary[i]][2] == " type=ACCELERATION")):
                                    BoundaryInit = BoundaryInit + [int(BCset[j][k]), int(CUBEN[(line_Boundary[i]+1+l)][1]), 0, 0, float(CUBEN[(line_Boundary[i]+1+l)][3])]


        BoundaryInit_Array = [[0 for x in range(5)] for y in range((len(BoundaryInit)/5))]
        
        j = 0
        for i in range(0,len(BoundaryInit),5):
            BoundaryInit_Array[j] = [BoundaryInit[i], BoundaryInit[i+1], BoundaryInit[i+2], BoundaryInit[i+3], BoundaryInit[i+4]]
            j = j + 1

        BoundaryInit_Final = [[0 for x in range(5)] for y in range((len(BoundaryInit)/5))]

        for i in range(0,len(BoundaryInit_Array)):
            if (BoundaryInit_Array[i][0] != 0):
                joint_ind = [(index, row.index(BoundaryInit_Array[i][0])) for index, row in enumerate(BoundaryInit_Array) if BoundaryInit_Array[i][0] in row]
                dir_ind = [(index, row.index(BoundaryInit_Array[i][1])) for index, row in enumerate(BoundaryInit_Array) if BoundaryInit_Array[i][1] in row]
                joint_col = [x[1] for x in joint_ind]
                joint_row = [x[0] for x in joint_ind]
                indices = [y for y, x in enumerate(joint_col) if x == 0]
                joint_lines = [joint_row[y] for y in indices]
                
                dir_col = [x[1] for x in dir_ind]
                dir_row = [x[0] for x in dir_ind]
                indices = [y for y, x in enumerate(dir_col) if x == 1]
                dir_lines = [dir_row[y] for y in indices]
                
                set_lines = list(set(joint_lines) & set(dir_lines))
                BoundaryInit_Final[i][0] = BoundaryInit_Array[i][0]
                BoundaryInit_Final[i][1] = BoundaryInit_Array[i][1]
                for j in set_lines:
                    if BoundaryInit_Array[j][2] != 0:
                        BoundaryInit_Final[i][2] = BoundaryInit_Array[j][2]
                    if BoundaryInit_Array[j][3] != 0:
                        BoundaryInit_Final[i][3] = BoundaryInit_Array[j][3]
                    if BoundaryInit_Array[j][4] != 0:
                        BoundaryInit_Final[i][4] = BoundaryInit_Array[j][4]
                    BoundaryInit_Array[j][0] = 0

        for i in range(0,len(BoundaryInit_Final)):
            if (BoundaryInit_Final[i][0] != 0):
                model_def.write('%d,%d,%g,%g,%g\n' % (BoundaryInit_Final[i][0], BoundaryInit_Final[i][1], BoundaryInit_Final[i][2], BoundaryInit_Final[i][3], BoundaryInit_Final[i][4]))

    #Terminate list in model_def.txt using 0,0,0,0,0
    model_def.write('0,0,0,0,0\n')

    numopt = input("Enter numerical damping scheme:\n[0 - none / 1 - Generalized-alpha / 2 - HHT method / 3 - WBZ method]\n>> ")
    if numopt > 0:
        spectrds = input("Spectral Radius:\n>> ")
    else:
        spectrds = 1

    if (spectrds < 0):
        print "Invalid value for spectral radius: cannot be less than 0."
        spectrds = input("Enter spectral radius:\n>> ")
    if (spectrds > 1):
        print "Invalid value for spectral radius: cannot be greater than 1."
        spectrds = input("Enter spectral radius:\n>> ")
    model_def.write('%d,%g\n' % (numopt,spectrds))

    if ALGFLAG == 5: #(Nonlinear Newmark Implicit Integration Method)
        print "Enter load proportionality factor parameters:"
        LPFmax = input("Maximum lambda:\n>> ")
        LPF = input("Initial lambda:\n>> ")
        dLPF = input("Increment of lambda:\n>> ")
        dLPFmax = input("Maximum increment of lambda:\n>> ")
        dLPFmin = input("Minimum increment of lambda:\n>> ")
        model_def.write('%g,%g,%g,%g,%g\n' % (LPFmax, LPF, dLPF, dLPFmax, dLPFmin))
        print "Enter maximums and minimums on counters:"
        itemax = input("Maximum number of iterations within load increment:\n>> ")
        submax = input ("Maximum number of times to step back load due to unconverged solution:\n>> ")
        solmin = input("Minimum number of converged solutions before increasing increment of lambda:\n>> ")
        model_def.write('%d,%d,%d\n' % (itemax,submax,solmin))
                       
#Enter tolerances on out-of-balance displacements, forces, and energy
if ANAFLAG == 1:
    model_def.write('0.001,0.001,0.001')
else:
    dtol,ftol,etol = map(float, raw_input("Enter Tolerances on Out-of-Balance Displacements, Forces, and Energy:\n[displacement tolerance, force tolerance, energy tolerance]\n>> ").split(','))
    model_def.write('%g,%g,%g' % (dtol, ftol, etol))
