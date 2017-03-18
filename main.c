//********************************************************************************
//**																			**
//**  Pertains to CU-BEN ver 3.14												**
//**																			**
//**  Copyright (c) 2016 C. J. Earls                                            **
//**  Developed by C. J. Earls, Cornell University                              **
//**  All rights reserved.														**
//**                                                                            **
//**  Contributors:                                                             **
//**    Christopher Stull                                                       **
//**    Heather Reed                                                            **
//**    Justyna Kosianka                                                        **
//**    Wensi Wu                                                                **
//**                                                                            **
//**  Redistribution and use in source and binary forms, with or without		**
//**  modification, are permitted provided that the following conditions are	**
//**  met:																		**
//**																			**
//**  - Redistributions of source code must retain the above copyright			**
//**    notice, this list of conditions and the following disclaimer.			**
//**																			**
//**  - Redistributions in binary form must reproduce the above copyright		**
//**    notice, this list of conditions and the following disclaimer listed		**
//**    in this license in the documentation and/or other materials				**
//**    provided with the distribution.											**
//**																			**
//**  - Neither the name of the copyright holders nor the name of Cornell		**
//**    University may be used to endorse or promote products derived from		**
//**    this software without specific prior written permission.				**
//**																			**
//**  Private, research, and institutional usage is without charge.				**
//**  Distribution of modified versions of this soure code is admissible		**
//**  UNDER THE CONDITION THAT THIS SOURCE CODE REMAINS UNDER COPYRIGHT OF		**
//**  THE ORIGINAL DEVELOPERS, BOTH SOURCE AND OBJECT CODE ARE MADE FREELY      **
//**  AVAILABLE WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE MODIFICATIONS.	**
//**  Distribution of this code as part of a commercial system is permissible	**
//**  ONLY BY DIRECT ARRANGEMENT WITH THE DEVELOPERS.							**
//**																			**
//**																			**
//**  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS		**
//**  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT			**
//**  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR		**
//**  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT		**
//**  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,		**
//**  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT			**
//**  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,		**
//**  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY		**
//**  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT		**
//**  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE		**
//**  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.		**
//**																			**
//********************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prototypes.h"


#include <Accelerate/Accelerate.h>
//#include "f2c.h"

/*
 CU-BEN Serial Version 3.14 (October 26, 2016)

 Analysis Options:
    1st order elastic, i.e. "linear elastic"
    2nd order elastic, i.e. "geometrically nonlinear"
    2nd order inelastic, i.e. "geometrically and materially nonlinear"
    Fluid-structure interaction (FSI)
 Nonlinear Solution Algorithm Options:
    (Static) Newton Raphson (NR)
    (Static) Modified Newton Raphson (MNR)
    (Static) Modified Spherical Arc Length (MSAL)
    (Dynamic) Newmark Implicit Integration Method
    (Dynamic) Nonlinear Newmark Implicit Integration Method
 Element Options (Updated Lagrangian):
    6-dof space trusses
    14-dof space frames
    18-dof triangular shells
    24-dof continuum brick (solid and fluid)

 Input data:
    enter flag for analysis type (in main) - ANAFLAG
        1 - 1st order elastic
        2 - 2nd order elastic
        3 - 2nd order inelastic
        *** for 1st order elastic analysis, all 2nd order / inelastic analysis related variables must be entered (for consistency in input files), but will be ignored during solution
        4 - Fluid-structure interaction
        *** for FSI analysis, also enter (on same line), flag for fsi incidence array - FSIINCFLAG tracks which elements have FSI nodes, then which face of the element (if the element is a brick) is on the FSI interface, and then the (global) nodes that are on the interface.  This is used later in the assembly if the global A and G matrices that translate the normal pressures into x, y, and z, displacements (i.e. the off-diagonal matrices in the monolithic K and M matrices.)
            0 - Input fsi incidence array in fsiinc.txt
            1 - Allow BEN to calculate fsi incidence array
    enter flag for solution algorithm type (in main) - ALGFLAG
        1 - (Static) Newton Raphson
        2 - (Static) Modified Newton Raphson
        3 - (Static) Modified Spherical Arc Length
        4 - (Dynamic) Newmark Implicit Integration Method
        5 - (Dynamic) Nonlinear Newmark Implicit Integration Method
        *** for 1st order elastic analysis, ALGFLAG is automatically set to 0
    enter flag for solver algorithm type (in main) - SLVFLAG
        0 - CU_BEN for symmetric matrices
        1 - CLAPACK solver for symmetric and non-symmetric matrices
    enter flag for execution of node-renumbering algorithm (in main) - optflag
        1 - no
        2 - yes
    enter number of joints (in main) - NJ
    enter number of elements (in main):
        truss - NE_TR
        frame - NE_FR
        shell - NE_SH
        solid brick - NE_SBR
        fluid brick - NE_FBR
        *** enter on single line as: NE_TR, NE_FR, NE_SH, NE_SBR, NE_FBR
    enter truss element member incidences (in struc) - minc[i,1],minc[i,2];
        i = 1 to NE_TR
    enter frame element member incidences (in struc) - minc[i,1],minc[i,2];
        i = 1 to NE_FR
    enter shell element member incidences (in struc) - minc[i,1],minc[i,2],minc[i,3];
        i = 1 to NE_SH
    enter solid brick element member incidences (in struc) - minc[i,1],minc[i,2],minc[i,3],minc[i,4],minc[i,5],minc[i,6],minc[i,7],minc[i,8];
        i = 1 to NE_SBR
    enter fluid brick element member incidences (in struc) - minc[i,1],minc[i,2],minc[i,3],minc[i,4],minc[i,5],minc[i,6],minc[i,7],minc[i,8];
        i = 1 to NE_FBR
    enter joint constraint(s) (in struc) - jnum,jdir; end = 0,0
        *** For acoustic fluid elements, constraint direction 7 corresponds with pressure
        *** for warping DOFs - jnum,jdir,restrnt
        0 - fixed
        1 - free
    enter joint coordinates (in prop) - x[i,1],x[i,2],x[i,3]; i = 1 to NJ
    if (NE_TR > 0) {
        enter truss element properties (in prop_tr); i = 1 to NE_TR:
            elastic modulus - emod[i]
            cross-sectional area - carea[i]
            density - dens[i]
            maximum allowable yield stress - yield[i]
            *** enter on single line as: emod[i],carea[i],dens[i],yield[i]
    }
    if (NE_FR > 0) {
        enter frame element material and geometric properties (in prop_fr); i = 1 to NE_FR:
            material properties
            elastic modulus - emod[i]
            shear modulus - gmod[i]
            density - den[i]
            geometric properties:
            cross-sectional area - carea[i]
            moment of inertia, strong-axis - istrong[i]
            moment of inertia, weak-axis - iweak[i]
            moment of inertia, polar - ipolar[i]
            moment of inertia, warping - iwarp[i]
            *** enter on single line as: emod[i],gmod[i],...,ipolar[i],iwarp[i]
        enter frame element member end offsets (in prop_fr) -
            member - osflag[i]
            End-1, Direction 1 - offset[i,1]
            End-1, Direction 2 - offset[i,2]
            End-1, Direction 3 - offset[i,3]
            End-2, Direction 1 - offset[i,4]
            End-2, Direction 2 - offset[i,5]
            End-2, Direction 3 - offset[i,6]
            *** enter on single line as:
                osflag[i], offset[i,1],offset[i,2],offset[i,3],offset[i,4],offset[i,5], offset[i,6]
                end = 0,0,0,0,0,0,0
        enter coordinates of frame element auxiliary point; plane formed by member end-points and the auxiliary point orients the strong-axis of the frame element (in prop_fr) - auxpt[i,1],auxpt[i,2],auxpt[i,3]; i = 1 to NE_FR
        enter element yield criteria (in prop_fr); i = 1 to NE_FR:
            maximum allowable yield stress - yield[i]
            plastic section modulus, strong-axis - zstrong[i]
            plastic section modulus, weak-axis - zweak[i]
            *** enter on single line as: yield[i],zstrong[i],zweak[i]
        enter frame element member end bending releases (in prop_fr):
            *** this only applies to frame elements for which member ends are released, i.e. do not enter frame elements which contain no member end bending releases
                0 - rigid
                1 - released
            member - mendrel[i,1]
            End-1, strong-axis - mendrel[i,2]
            End-1, weak-axis - mendrel[i,3]
            End-2, strong-axis - mendrel[i,4]
            End-2, weak-axis - mendrel[i,5]
            *** enter on single line as:
                mendrel[i,1],mendrel[i,2],mendrel[i,3],mendrel[i,4],mendrel[i,5];
                end = 0,0,0,0,0
    }
    if (NE_SH > 0) {
        enter shell element properties (in prop_sh); i = 1 to NE_SH:
            elastic modulus - emod[i]
            Poisson's Ratio - nu[i]
            thickness - thick[i]
            density - dens[i]
            maximum allowable yield stress - yield[i]
            *** enter on single line as: emod[i],nu[i],thick[i],dens[i],yield[i]
    }
    if (NE_SBR > 0) {
        enter brick element properties (in prop_fsi); i = 1 to NE_SBR:
            elastic modulus - emod[i]
            Poisson's Ratio - nu[i]
            density (solid) - dens[i]
            maximum allowable yield stress - yield[i]
            *** enter on single line as: emod[i],nu[i],dens[i],yield[i]
    }
    if (NE_FBR > 0) {
        enter fluid brick element properties (in prop_fsi); ONLY ONCE:
            density (fluid) - fdens
            bulk modulus - bmod
    }
    if (ANAFLAG == 4) { // For FSI analysis
        enter joint absorbtion areas (enter zero if joint is non-absorbing) (in prop_fsi); i = 1 to NJ
        enter joint normal orientation points to orient the direction of the outward facing normal vector of the interface (in prop_fsi); i = 1 to NJ
            *** enter on single line as: norpt[i,1],norpt[i,2],norpt[i,3]
        enter number of time steps (in main) and total time for analysis (s);
            *** enter on single line as: ntstpsinpt, ttot
        enter concentrated load, nodal acceleration(s) and fluid incident pressure(s) applied during time stepon joints for each time step (in load_fsi) - i = 0:ntstpsinpt
            joint,dir,force,fpress,facc;
            end = 0,0,0,0,0
        enter initial conditions for node displacement or pressure and 1st or 2nd derivatives.  If dir = 4, initial condition refers to a fluid DOF (in load_fsi)
            joint,dir,disp,vel,acc;
            end = 0,0,0,0,0
        enter alpha and beta parameters for Newmark time integration scheme
            alpha, delta
    }
    else { // non-FSI analysis
        if (ALGFLAG < 4) { //Static
            enter reference concentrated load(s) on joints (in load) - joint,dir,force;
            end = 0,0,0
            if (NE_FR > 0) {
                enter reference uniformly distributed load(s) on frame elements (in load) -
                    frame,dir,force; end = 0,0,0
            }
        }
        else { // Dynamic Analysis
            if (ALGFLAG == 4){
                enter initial number of time steps (in load) and total time for analysis (s);
                    *** enter on single line as: ntstpsinpt, ttot
                enter reference concentrated load(s) on joints for each time step (in load); i = 0:ntstpsinpt
                    joint,dir,force;
                    end = 0,0,0
                enter all non-zero initial conditions for node displacement, velocity, and acceleration (in load)
                    joint,dir,disp,vel,acc;
                    end = 0,0,0,0,0
                enter Newmark time integration scheme option and spectral radius
                    numopt,spectrds
                    0, 1 - standard Newmark
                    1, 1 - standard Newmark
                    1, 0.9 - Generalized-alpha method with spectral radius of 0.9
                    2, 1 - standard Newmark
                    2, 0.9 - HHT method with spectral radius of 0.9
                    3, 1 - standard Newmark
                    3, 0.9 - WBZ method with spectral radius of 0.9
            }
            if (ALGFLAG == 5){
                The load history is under linear interpolation assumption in the case when damping scheme is applied.
                Limit the size of delta T that may be used to maintain fidelity with the loading history.
                enter initial number of time steps (in load) and total time for analysis (s);
                    *** enter on single line as: ntstpsinpt, ttot
                enter reference concentrated load(s) on joints for each time step (in load); i = 0:ntstpsinpt
                    joint,dir,force;
                    end = 0,0,0
                enter all non-zero initial conditions for node displacement, velocity, and acceleration (in load)
                    joint,dir,disp,vel,acc;
                    end = 0,0,0,0,0
                enter Newmark time integration scheme option and spectral radius
                    numopt,spectrds
                    0, 1 - standard Newmark
                    1, 1 - standard Newmark
                    1, 0.9 - Generalized-alpha method with spectral radius of 0.9
                    2, 1 - standard Newmark
                    2, 0.9 - HHT method with spectral radius of 0.9
                    3, 1 - standard Newmark
                    3, 0.9 - WBZ method with spectral radius of 0.9
                enter load proportionality factor parameters (in main):
                    maximum lambda - lpfmax
                    initial lambda - lpf
                    increment of lambda - dlpf
                    maximum increment of lambda - dlpfmax
                    minimum increment of lambda - dlpfmin
                    *** enter on single line as: lpfmax,lpf,dlpf,dlpfmax,dlpfmin
                enter maximums / minimums on counters (in main):
                    maximum number of iterations within load increment - itemax
                    maximum number of times to step back load due to unconverged solution - submax
                    minimum number of converged solutions before increasing increment of lambda - solmin
                    *** enter on single line as: itemax,submax,solmin
            }
        }
        if (ALGFLAG == 1) {
            enter maximum load proportionality factor (in main) - lpfmax
        }
        else if (ALGFLAG == 2) {
            enter load proportionality factor parameters (in main):
                maximum lambda - lpfmax
                initial lambda - lpf
                increment of lambda - dlpf
                maximum increment of lambda - dlpfmax
                minimum increment of lambda - dlpfmin
                *** enter on single line as: lpfmax,lpf,dlpf,dlpfmax,dlpfmin
            enter maximums / minimums on counters (in main):
                maximum number of iterations within load increment - itemax
                maximum number of times to step back load due to unconverged solution - submax
                minimum number of converged solutions before increasing increment of lambda - solmin
                *** enter on single line as: itemax,submax,solmin
        }
        else { //MSAL
            enter MSAL parameters:
                initial prescribed displacement at DOF "k" (in msal) - jnum,jdir,dk
                factor limiting size of load increment (in main) - alpha
                threshold to eliminate load control in the arc length criterion in the vicinity of a critical point; a minimum value of 0.1 is recommended (in main) - psi_thresh
                optimum number of iterations; a value of 6 is recommended (in main) - iteopt
                maximum lambda and allowable displacment (in main) - lpfmax,dkimax
            enter maximums / minimums on counters (in main):
                maximum number of iterations within load increment - itemax
                maximum number of times to step back load due to unconverged solution - submax
                maximum number of times to step back load due to arc length criterion producing imaginary roots - imagmax
                maximum number of times to step back load due to arc length criterion producing two negative roots - negmax
                *** enter on single line as: itemax,submax,imagmax,negmax
        }
    }
    enter tolerances on out-of-balance displacements, forces, and energy (in main) -
        toldisp,tolforc,tolener
 */

// Number of joints, number of truss, frame, and shell elements, and number of equations
long NJ, NE_TR, NE_FR, NE_SH, NE_SBR, NE_FBR, NE_BR, NEQ, SNDOF, FNDOF, NTSTPS, ntstpsinpt;
double dt, ttot;
// "666" is an unlikely mistake; initialization allows for assumption of empty input file
int ANAFLAG = 666, ALGFLAG, OPTFLAG, SLVFLAG, FSIFLAG, FSIINCFLAG, brFSI_FLAG, shFSI_FLAG;
FILE *IFP[3], *OFP[7]; // Pointers to input and output file

int main (int argc, char **argv)
{
    int i, j; // Counter variables

    // Open I/O for business!
    do {
        IFP[0] = fopen("model_def.txt", "r"); // Open input file for reading
    } while (IFP[0] == 0);

    char file[20];
    for (i = 0; i < 7; ++i) {
        sprintf(file, "results%d.txt", i + 1);
        do {
            OFP[i] = fopen(file, "w"); // Open output file for writing
        } while (OFP[i] == 0);
    }

    // Read in analysis / algorithm /solver type from input file
    fscanf(IFP[0], "%d,", &ANAFLAG);
    if (ANAFLAG == 4) {
        fscanf(IFP[0], "%d\n", &FSIINCFLAG);
    }
    fscanf(IFP[0], "%d\n", &ALGFLAG);
    fscanf(IFP[0], "%d\n", &SLVFLAG);
    if (ANAFLAG == 1 && ALGFLAG != 4) {
        fprintf(OFP[0], "Analysis Type:\n\t1st Order Elastic\n");
        ALGFLAG = 0;
    } else if (ANAFLAG == 2) {
        fprintf(OFP[0], "Analysis Type:\n\t2nd Order Elastic\n");
    } else if (ANAFLAG == 3) {
        fprintf(OFP[0], "Analysis Type:\n\t2nd Order Inelastic\n");
    } else if (ANAFLAG == 1 && ALGFLAG == 4) {
        fprintf(OFP[0], "Analysis Type:\n\t1stOrder Elastic Dynamic\n");
    } else if (ANAFLAG == 2 && ALGFLAG == 5) {
        fprintf(OFP[0], "Analysis Type:\n\t2nd Order Elastic Dynamic\n");
    } else if (ANAFLAG == 3 && ALGFLAG == 5) {
        fprintf(OFP[0], "Analysis Type:\n\t2nd Order Inelastic Dynamic\n");
    }else if (ANAFLAG == 4 && ALGFLAG != 5) {
        fprintf(OFP[0], "Analysis Type:\n\tFluid Structure Interaction\n");
    } else if (ANAFLAG == 666) {
        fprintf(OFP[0], "***ERROR*** Input file is empty\n");
        goto EXIT1;
    } else {
        fprintf(OFP[0], "***ERROR*** Invalid entry for analysis type\n");
        goto EXIT1;
    }

    if (ALGFLAG == 0) {
        fprintf(OFP[0], "\nAlgorithm Type:\n\tDirect Stiffness\n");
    } else if (ALGFLAG == 1) {
        fprintf(OFP[0], "\nAlgorithm Type:\n\tNewton Raphson\n");
    } else if (ALGFLAG == 2) {
        fprintf(OFP[0], "\nAlgorithm Type:\n\tModified Newton Raphson\n");
    } else if (ALGFLAG == 3) {
        fprintf(OFP[0], "\nAlgorithm Type:\n\tModified Spherical Arc Length\n");
    } else if (ALGFLAG == 4 || ALGFLAG == 5)   {
        fprintf(OFP[0], "\nAlgorithm Type:\n\tDynamic (Newmark)\n");
    } else {
        fprintf(OFP[0], "\n***ERROR*** Invalid entry for algorithm type\n");
        goto EXIT1;
    }

    // Read in optimization flag, number of joints and elements from input file
    fscanf(IFP[0], "%d\n", &OPTFLAG);
    fscanf(IFP[0], "%ld\n", &NJ);
    fscanf(IFP[0], "%ld,%ld,%ld,%ld,%ld\n", &NE_TR, &NE_FR, &NE_SH, &NE_SBR, &NE_FBR);
    fprintf(OFP[0], "\nControl Variables:\n\tNumber of Joints: %ld\n", NJ);
    fprintf(OFP[0], "\tNumber of Truss Elements: %ld\n", NE_TR);
    fprintf(OFP[0], "\tNumber of Frame Elements: %ld\n", NE_FR);
    fprintf(OFP[0], "\tNumber of Shell Elements: %ld\n", NE_SH);
    fprintf(OFP[0], "\tNumber of Solid Brick Elements: %ld\n", NE_SBR);
    fprintf(OFP[0], "\tNumber of Fluid Brick Elements: %ld\n", NE_FBR);

    // Total number of bricks
    NE_BR = NE_SBR + NE_FBR;

    // Set up flags to determine whether the solid elements are shells or bricks
    if (NE_SH > NE_SBR) {
        shFSI_FLAG = 1;
    }
    else {
        brFSI_FLAG = 1;
    }

    if (OPTFLAG == 2) {
        // Open I/O for business!
        do {
            // Open optimized input file for writing
            IFP[1] = fopen("model_def_OPT.txt", "w");
        } while (IFP[1] == 0);

        // Write control variables to optimized input file
        fprintf(IFP[1], "%d\n%d\n%d\n1\n%ld\n%ld,%ld,%ld,%ld,%ld\n", ANAFLAG, ALGFLAG, SLVFLAG, NJ,
                NE_TR, NE_FR, NE_SH, NE_BR, NE_FBR);
    }



    // Memory management variables
    /* Pointer-to-pointer-to-int array (5 arrays of type int are defined during program
     execution) */
    int *p2p2i[5];
    // Counter to track number of arrays of type int for which memory is allocated
    int ni = 0;
    /* Pointer-to-pointer-to-long array (7 arrays of type long are defined during program
     execution) */
    long *p2p2l[7];
    // Counter to track number of arrays of type long for which memory is allocated
    int nl = 0;
    /* Pointer-to-pointer-to-double array (101 arrays of type double are defined during
     program execution) */
    double *p2p2d[101];
    // Counter to track number of arrays of type double for which memory is allocated
    int nd = 0;

    /*
     Define secondary variables which DO NOT depend upon NEQ, common to both NR and MSAL
     algorithms
     Note: variables for which the reference configuration changes during the analysis
     take on the following notation:
     xyz = configuration at last successful load increment
     xyz_i = current configuration
     xyz_ip = configuration immediately preceding current configuration
     Note: variables with "_temp" are temporary variables which may be reverted back
     to permanent counterparts upon unsuccessful / invalid load increment
     */

    /*
     Variables related to joint coordinates
     */
    // Joint coordinates
    double *x = alloc_dbl (NJ*3);
    if (x == NULL) {
        goto EXIT1;
    }
    p2p2d[nd] = x;
    nd++;
    double *x_temp = alloc_dbl (NJ*3);
    if (x_temp == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = x_temp;
    nd++;
    double *x_ip = alloc_dbl (NJ*3);
    if (x_ip == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = x_ip;
    nd++;

    // Auxiliary points (frame element only)
    double *auxpt = alloc_dbl (NE_FR*3);
    if (auxpt == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = auxpt;
    nd++;

    // Member end offsets (frame element only)
    double *offset = alloc_dbl (NE_FR*6);
    if (offset == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = offset;
    nd++;
    // Flag for member end offsets (frame elements only)
    int *osflag = alloc_int (NE_FR);
    if (osflag == NULL) {
        goto EXIT2;
    }
    p2p2i[ni] = osflag;
    ni++;

    // Frame element member end coordinates (frame elements only)
    double *xfr = alloc_dbl (NE_FR*6);
    if (xfr == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = xfr;
    nd++;
    double *xfr_temp = alloc_dbl (NE_FR*6);
    if (xfr_temp == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = xfr_temp;
    nd++;

    // Non-zero local coordinates (shell element only)
    double *xlocal = alloc_dbl (NE_SH*3);
    if (xlocal == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = xlocal;
    nd++;

    /*
     Variables related to element material properties
     */
    // Elastic modulus
    double *emod = alloc_dbl (NE_TR+NE_FR+NE_SH+NE_BR);
    if (emod == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = emod;
    nd++;
    // Solid density
    double *dens = alloc_dbl (NE_TR+NE_FR+NE_SH+NE_BR);
    if (dens == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = dens;
    nd++;
    // Fluid density
    double *fdens = alloc_dbl (1);
    if (fdens == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = fdens;
    nd++;
    // Fluid bulk modulus
    double *bmod = alloc_dbl (1);
    if (bmod == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = bmod;
    nd++;
    // Shear modulus (frame element only)
    double *gmod = alloc_dbl (NE_FR+NE_BR);
    if (gmod == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = gmod;
    nd++;
    // Poisson's Ratio (shell element only)
    double *nu = alloc_dbl (NE_SH+NE_BR);
    if (nu == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = nu;
    nd++;

    /*
     Variables related to element geometrical properties
     */
    // Cross-sectional area (truss and frame elements only)
    double *carea = alloc_dbl (NE_TR+NE_FR+NE_BR);
    if (carea == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = carea;
    nd++;

    // Thickness (shell elements only)
    double *thick = alloc_dbl (NE_SH);
    if (thick == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = thick;
    nd++;

    // Face area and deformed face area (shell elements only)
    double *farea = alloc_dbl (NE_SH);
    if (farea == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = farea;
    nd++;
    double *deffarea = alloc_dbl (NE_SH);
    if (deffarea == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = deffarea;
    nd++;
    double *deffarea_i = alloc_dbl (NE_SH);
    if (deffarea_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = deffarea_i;
    nd++;
    double *deffarea_ip = alloc_dbl (NE_SH);
    if (deffarea_ip == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = deffarea_ip;
    nd++;

    // Longitudinal length and deformed longitudinal length (truss and frame elements only)
    double *llength = alloc_dbl (NE_TR+NE_FR);
    if (llength == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = llength;
    nd++;
    double *llength_temp = alloc_dbl (NE_TR+NE_FR);
    if (llength_temp == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = llength_temp;
    nd++;
    double *defllen = alloc_dbl (NE_TR+NE_FR);
    if (defllen == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = defllen;
    nd++;
    double *defllen_i = alloc_dbl (NE_TR+NE_FR);
    if (defllen_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = defllen_i;
    nd++;
    double *defllen_ip = alloc_dbl (NE_TR+NE_FR);
    if (defllen_ip == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = defllen_ip;
    nd++;

    // Side-length and deformed side-length (shell elements only)
    double *slength = alloc_dbl (NE_SH*3);
    if (slength == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = slength;
    nd++;
    double *defslen = alloc_dbl (NE_SH*3);
    if (defslen == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = defslen;
    nd++;
    double *defslen_i = alloc_dbl (NE_SH*3);
    if (defslen_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = defslen_i;
    nd++;
    double *defslen_ip = alloc_dbl (NE_SH*3);
    if (defslen_ip == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = defslen_ip;
    nd++;

    // Strong-axis moment of inertia (frame elements only)
    double *istrong = alloc_dbl (NE_FR);
    if (istrong == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = istrong;
    nd++;
    // Weak-axis moment of inertia (frame elements only)
    double *iweak = alloc_dbl (NE_FR);
    if (iweak == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = iweak;
    nd++;
    // Polar moment of inertia (frame elements only)
    double *ipolar = alloc_dbl (NE_FR);
    if (ipolar == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = ipolar;
    nd++;
    // Warping moment of inertia (frame elements only)
    double *iwarp = alloc_dbl (NE_FR);
    if (iwarp == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = iwarp;
    nd++;

    // Direction cosines
    double *c1 = alloc_dbl (NE_TR+NE_FR*3+NE_SH*3);
    if (c1 == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = c1;
    nd++;
    double *c2 = alloc_dbl (NE_TR+NE_FR*3+NE_SH*3);
    if (c2 == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = c2;
    nd++;
    double *c3 = alloc_dbl (NE_TR+NE_FR*3+NE_SH*3);
    if (c3 == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = c3;
    nd++;
    double *c1_i = alloc_dbl (NE_TR+NE_FR*3+NE_SH*3);
    if (c1_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = c1_i;
    nd++;
    double *c2_i = alloc_dbl (NE_TR+NE_FR*3+NE_SH*3);
    if (c2_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = c2_i;
    nd++;
    double *c3_i = alloc_dbl (NE_TR+NE_FR*3+NE_SH*3);
    if (c3_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = c3_i;
    nd++;
    double *c1_ip = alloc_dbl (NE_TR+NE_FR*3+NE_SH*3);
    if (c1_ip == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = c1_ip;
    nd++;
    double *c2_ip = alloc_dbl (NE_TR+NE_FR*3+NE_SH*3);
    if (c2_ip == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = c2_ip;
    nd++;
    double *c3_ip = alloc_dbl (NE_TR+NE_FR*3+NE_SH*3);
    if (c3_ip == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = c3_ip;
    nd++;


    /*
     Variables related to element yield criteria
     */
    // Maximum allowable yield stress
    double *yield = alloc_dbl (NE_TR+NE_FR+NE_SH+NE_BR);
    if (yield == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = yield;
    nd++;

    // Yield flag (frame elements only)
    int *yldflag = alloc_int (NE_FR*2);
    if (yldflag == NULL) {
        goto EXIT2;
    }
    p2p2i[ni] = yldflag;
    ni++;

    // Strong-axis section modulus (frame elements only)
    double *zstrong = alloc_dbl (NE_FR);
    if (zstrong == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = zstrong;
    nd++;
    // Weak-axis section modulus (frame elements only)
    double *zweak = alloc_dbl (NE_FR);
    if (zweak == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = zweak;
    nd++;

    // Equivalent plastic curvatures (shell elements only)
    double *chi = alloc_dbl (NE_SH*3);
    if (chi == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = chi;
    nd++;
    double *chi_temp = alloc_dbl (NE_SH*3);
    if (chi_temp == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = chi_temp;
    nd++;

    /*
     Variables related to element forces
     */
    // Internal force vectors
    double *ef = alloc_dbl (NE_TR*2+NE_FR*14+NE_SH*18);
    if (ef == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = ef;
    nd++;
    double *ef_i = alloc_dbl (NE_TR*2+NE_FR*14+NE_SH*18);
    if (ef_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = ef_i;
    nd++;
    double *ef_ip = alloc_dbl (NE_TR*2+NE_FR*14+NE_SH*18);
    if (ef_ip == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = ef_ip;
    nd++;

    // Reference fixed-end force vectors (frame elements only)
    double *efFE_ref = alloc_dbl (NE_FR*14);
    if (efFE_ref == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = efFE_ref;
    nd++;
    // Cumulative fixed-end force vectors (frame elements only)
    double *efFE = alloc_dbl (NE_FR*14);
    if (efFE == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = efFE;
    nd++;
    double *efFE_i = alloc_dbl (NE_FR*14);
    if (efFE_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = efFE_i;
    nd++;
    double *efFE_ip = alloc_dbl (NE_FR*14);
    if (efFE_ip == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = efFE_ip;
    nd++;

    // Internal membrane forces (shell elements only)
    double *efN = alloc_dbl (NE_SH*9);
    if (efN == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = efN;
    nd++;
    double *efN_temp = alloc_dbl (NE_SH*9);
    if (efN_temp == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = efN_temp;
    nd++;
    // Internal bending moments (shell elements only)
    double *efM = alloc_dbl (NE_SH*9);
    if (efM == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = efM;
    nd++;
    double *efM_temp = alloc_dbl (NE_SH*9);
    if (efM_temp == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = efM_temp;
    nd++;
    double *Jinv = alloc_dbl (3*3*8*NE_BR); // 3x3 matrix for each integration point of each brick
    if (Jinv == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = Jinv;
    nd++;


    /*
     Miscellaneous variables
     */
    // Member end bending releases (frame elements only)
    int *mendrel = alloc_int (NE_FR*5);
    if (mendrel == NULL) {
        goto EXIT2;
    }
    p2p2i[ni] = mendrel;
    ni++;

    /*
     Finite element model related variables
     */
    // Optimized node numbering scheme
    long *jnt = alloc_long (NJ);
    if (jnt == NULL) {
        goto EXIT2;
    }
    p2p2l[nl] = jnt;
    nl++;
    // Member global DOF code
    long *mcode = alloc_long (NE_TR*6+NE_FR*14+NE_SH*18+NE_SBR*24+NE_FBR*24);
    if (mcode == NULL) {
        goto EXIT2;
    }
    p2p2l[nl] = mcode;
    nl++;
    // Joint constraint code
    long *jcode = alloc_long (NJ*7);
    if (jcode == NULL) {
        goto EXIT2;
    }
    p2p2l[nl] = jcode;
    nl++;
    // Member incidences
    long *minc = alloc_long (NE_TR*2+NE_FR*2+NE_SH*3+NE_BR*8);
    if (minc == NULL) {
        goto EXIT2;
    }
    p2p2l[nl] = minc;
    nl++;
    // Warping restraint flags
    int *wrpres = alloc_int (NJ*3);
    if (wrpres == NULL) {
        goto EXIT2;
    }
    p2p2i[ni] = wrpres;
    ni++;

    /*
     Fluid-structure interaction related variables
     */
    // FSI incidences
    long *fsiinc = alloc_long (NE_SBR*6*4+NE_SH*3*1e5);
    if (fsiinc == NULL) {
        goto EXIT2;
    }
    p2p2l[nl] = fsiinc;
    nl++;
    // Number of FSI faces per solid fsi element
    long *elface = alloc_long (NE_SBR+NE_SH);
    if (elface == NULL) {
        goto EXIT2;
    }
    p2p2l[nl] = elface;
    nl++;
    // Absorbing joint areas
    double *abspt = alloc_dbl (NJ);
    if (abspt == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = abspt;
    nd++;
    // Normal orientation points
    double *norpt = alloc_dbl (NJ*3); // points that orient the f-s normals otuward from the structure
    if (norpt == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = norpt;
    nd++;

    double *nnorm = alloc_dbl (NJ*3); // global normal vecs for f-s joint
    if (nnorm == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = nnorm;
    nd++;
    double *jac = alloc_dbl (9); // global normal vecs for f-s joint
    if (jac == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = jac;
    nd++;
    double *tarea = alloc_dbl (NJ); // tributary areas for f-s interface
    if (tarea == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = tarea;
    nd++;

    // Pass control to struc function
    int errchk = struc (jcode, minc, wrpres, jnt);

    // Terminate program if errors encountered
    if (errchk == 1) {
        goto EXIT2;
    }



    // Pass control to codes function
    codes (mcode, jcode, minc, wrpres);

    // Print number of equations
    fprintf(OFP[0], "\nNumber of equations (system DOFs): %ld\n", NEQ);

    // Pass control to fsi function
    if (ANAFLAG == 4) {
        fsi (mcode, jcode, minc, elface, fsiinc);
    }

    /*
     Define secondary variables which DO depend upon NEQ, common to both NR and MSAL
     algorithms
     Note: variables for which the reference configuration changes during the analysis
     take on the following notation:
     xyz = configuration at last successful load increment
     xyz_ip = configuration immediately preceding current configuration
     Note: variables with "_temp" are temporary variables which may be reverted back
     to permanent counterparts upon occurence of unsuccessful load increment
     */

    /* Only allocate memory to variables that will be used
     - depends on analysis/algorithm */
    long NEQ_nonlin = 0;
    long NEQ_dyn = 0;
    long NEQ_FSI = 0;

    if (ANAFLAG == 2 || ANAFLAG == 3 || ALGFLAG == 3) {NEQ_nonlin = NEQ;}
    if (ALGFLAG == 4 || ALGFLAG == 5) {NEQ_dyn = NEQ;}
    if (ANAFLAG == 4) {NEQ_FSI = NEQ;}

    double *q = alloc_dbl (NEQ); // Generalized joint reference load vector
    if (q == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = q;
    nd++;
    // Generalized joint total load vector, i.e. lpf * q[NEQ]
    double *qtot = alloc_dbl (NEQ);
    if (qtot == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = qtot;
    nd++;

    //Generalized dynamic equivalent external load vector
    double *dyn = alloc_dbl (NEQ);
    if (dyn == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = dyn;
    nd++;

    // Total and incremental generalized nodal displacement vectors
    double *d = alloc_dbl (NEQ);
    if (d == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = d;
    nd++;
    double *dd = alloc_dbl (NEQ);
    if (dd == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = dd;
    nd++;
    double *d_temp = alloc_dbl (NEQ);
    if (d_temp == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = d_temp;
    nd++;
    // Generalized internal force vector
    double *f = alloc_dbl (NEQ);
    if (f == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = f;
    nd++;
    double *f_temp = alloc_dbl (NEQ_nonlin);
    if (f_temp == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = f_temp;
    nd++;
    // Generalized internal force vector from previous load increment
    double *fp = alloc_dbl (NEQ_nonlin);
    if (fp == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = fp;
    nd++;
    double *f_ip = alloc_dbl (NEQ_nonlin);
    if (f_ip == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = f_ip;
    nd++;
    double *r = alloc_dbl (NEQ); // Residual force vector, i.e. qtot - f
    if (r == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = r;
    nd++;
    // Array for CLAPACK
    int *ipiv = alloc_int (NEQ);
    if (ipiv == NULL) {
        goto EXIT2;
    }
    p2p2i[ni] = ipiv;
    ni++;

    // Skyline storage parameters for stiffness matrix
    long *maxa = alloc_long (NEQ+1);
    if (maxa == NULL) {
        goto EXIT2;
    }
    p2p2l[nl] = maxa;
    nl++;
    long lss;
    // Full system stiffness matrix
    double *ss_fsi = alloc_dbl (NEQ_FSI*NEQ_FSI);
    if (ss_fsi == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = ss_fsi;
    nd++;
    // Full system mass matrix
    double *sm_fsi = alloc_dbl (NEQ_FSI*NEQ_FSI);
    if (sm_fsi == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = sm_fsi;
    nd++;
    // System daming array
    double *sd_fsi = alloc_dbl (NEQ_FSI);
    if (sd_fsi == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = sd_fsi;
    nd++;
    // Array of input fluid pressures
    double *fpres = alloc_dbl (FNDOF);
    if (fpres == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = fpres;
    nd++;
    // Array of input fluid accelerations
    double *facc = alloc_dbl (FNDOF);
    if (facc == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = facc;
    nd++;
    // Array of previous displacements
    double *um = alloc_dbl (NEQ);
    if (um == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = um;
    nd++;
    // Array of previous velocities
    double *vm = alloc_dbl (NEQ);
    if (vm == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = vm;
    nd++;
    // Array of previous accelerations
    double *am = alloc_dbl (NEQ);
    if (am == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = am;
    nd++;
    // Array of next displacements
    double *uc = alloc_dbl (NEQ);
    if (uc == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = uc;
    nd++;
    // Array of next velocities
    double *vc = alloc_dbl (NEQ);
    if (vc == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = vc;
    nd++;
    // Array of next accelerations
    double *ac = alloc_dbl (NEQ);
    if (ac == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = ac;
    nd++;
    // Array of displacements at current iteration
    double *uc_i = alloc_dbl (NEQ);
    if (uc_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = uc_i;
    nd++;
    // Array of velocities at current iteration
    double *vc_i = alloc_dbl (NEQ);
    if (vc_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = vc_i;
    nd++;
    // Array of accelerations at current iteration
    double *ac_i = alloc_dbl (NEQ);
    if (ac_i == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = ac_i;
    nd++;
    // L matrix
    double *L = alloc_dbl (SNDOF*FNDOF); // L matrix = G*A
    if (L == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = L;
    nd++;
    double *LT = alloc_dbl (FNDOF*SNDOF); // L transpose
    if (LT == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = LT;
    nd++;
    double *G = alloc_dbl (SNDOF*FNDOF); // Matrix of direction cosines
    if (G == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = G;
    nd++;
    // Diagonal area matrix
    double *A = alloc_dbl (FNDOF*FNDOF); //
    if (A == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = A;
    nd++;
    // Effective load vector (for use in dynamic analysis)
    double *Reff = alloc_dbl (NEQ); //
    if (Reff == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = Reff;
    nd++;
    // Effective mass vector (for use in dynamic analysis)
    double *Meff = alloc_dbl (NEQ); //
    if (Meff == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = Meff;
    nd++;

    // Newmark integration constants
    double alpha, delta;

    // Pass control to skylin function
    errchk = skylin (maxa, mcode, &lss);

    // Print length of stiffness array
    fprintf(OFP[0], "\nLength of stiffness array: %ld\n", lss);

    // Terminate program if errors encountered
    if (errchk == 1) {
        goto EXIT2;
    }

    //Define variable which depends upon lss
    // Effective stiffness matrix (for use in dynamic analysis)
    double *Keff = alloc_dbl (lss); //
    if (Keff == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = Keff;
    nd++;

    // Define secondary variable which depends upon lss
    // Generalized stiffness array
    double *ss = alloc_dbl (lss);
    if (ss == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = ss;
    nd++;

    // Define secondary variable which depends upon lss
    // Generalized mass array
    double *sm = alloc_dbl (lss);
    if (sm == NULL) {
        goto EXIT2;
    }
    p2p2d[nd] = sm;
    nd++;

    // Read in joint coordinates from input file
    fprintf(OFP[0], "\nJoint Coordinates:\n\tJoint\t\tDirection-1\tDirection-2\t");
    fprintf(OFP[0], "Direction-3\n");
    for (i = 0; i < NJ; ++i) {
        fscanf(IFP[0], "%lf,%lf,%lf\n", &x[jnt[i]*3], &x[jnt[i]*3+1], &x[jnt[i]*3+2]);
    }
    for (i = 0; i < NJ; ++i) {
        fprintf(OFP[0], "\t%d\t\t%lf\t%lf\t%lf\n", i + 1, x[i*3], x[i*3+1], x[i*3+2]);
    }
    if (OPTFLAG == 2) {
        for (i = 0; i < NJ; ++i) {
            fprintf(IFP[1], "%lf,%lf,%lf\n", x[i*3], x[i*3+1], x[i*3+2]);
        }
    }

    fprintf(OFP[0], "\nJoint Degrees of Freedom:\n\tJoint\t\tX-Translation\t");
    fprintf(OFP[0], "Y-Translation\tZ-Translation\tX-Rotation\tY-Rotation\t");
    fprintf(OFP[0], "Z-Rotation\tWarping\n");
    for (i = 0; i < NJ; ++i) {
        fprintf(OFP[0], "\t%d\t\t", i + 1);
        for (j = 0; j < 7; ++j) {
            fprintf(OFP[0], "%ld\t\t", jcode[i*7+j]);
        }
        fprintf(OFP[0], "\n");
    }

    if (NE_TR > 0) {
        // Pass control to prop_tr function
        prop_tr (x, emod, carea, dens, llength, yield, c1, c2, c3, minc);

    }

    if (NE_FR > 0) {
        // Pass control to prop_fr function
        prop_fr (x, xfr, emod, gmod, dens, offset, osflag, auxpt, carea, llength, istrong,
                 iweak, ipolar, iwarp, yield, zstrong, zweak, c1, c2, c3, mendrel, minc);
    }

    if (NE_SH > 0) {
        // Pass control to prop_sh function
        prop_sh (x, emod, nu, xlocal, thick, dens, farea, slength, yield, c1, c2, c3, minc);
    }

    if (NE_BR > 0) {
        // Pass control to prop_br function
        prop_br (emod, nu, yield, dens, fdens, bmod, farea);
    }

    // If fluid-structure interaction analysis
    if (ANAFLAG == 4){
        // Pass control to prop_fsi function
        prop_fsi (x, emod, nu, dens, fdens, bmod, farea, slength, yield, minc, elface, fsiinc, nnorm, tarea,
                  ss, ss_fsi, sd_fsi, abspt, norpt, mcode, jcode, L);

        // Scan in the user desired number of time steps and total analysis
        fscanf(IFP[0], "%ld,%lf\n", &ntstpsinpt, &ttot);

        // Calculate dt
        dt = ttot/ntstpsinpt;
        ntstpsinpt += 1;

        // Allocate memory to arrays of input times, applied forces, pressures and accelerations
        double *tinpt = alloc_dbl (ntstpsinpt); // Time array
        if (tinpt == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = tinpt;
        nd++;

        double *pinpt = alloc_dbl (SNDOF*ntstpsinpt); // Applied mechanical forces acting on solid nodes
        if (pinpt == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = pinpt;
        nd++;

        double *presinpt = alloc_dbl (FNDOF*ntstpsinpt); // Applied fluid pressures acting on fluid nodes
        if (presinpt == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = presinpt;
        nd++;

        double *accinpt = alloc_dbl (FNDOF*ntstpsinpt); // Applied normal, incident fluid pressures acting on fluid nodes
        if (accinpt == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = accinpt;
        nd++;

        // Pass control to the L_br function
        L_br(minc, mcode, jcode, jcode, nnorm, tarea, L, A, G);

        // Initialize previous displacements, velocities, and accelerations
        for (i = 0; i < NEQ; ++i) {
            um[i] = 0; vm[i] = 0; am[i] = 0;
        }

        // Pass control to load_fsi function
        load_fsi (jcode, tinpt, pinpt, presinpt, accinpt, fdens, um, vm, am);

        /* Allocate memory to arrays of linearly interpolated loads, pressures and accelerations
         based on the actual time step for transient analysis */
        double *tstps = alloc_dbl (NTSTPS); // Time array based on actual dt
        if (tstps == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = tstps;
        nd++;

        double *qdyn = alloc_dbl (NEQ*NTSTPS); // Array of external agencies acting on the structure/fluid
        if (qdyn == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = qdyn;
        nd++;

        double *apload = alloc_dbl (SNDOF*NTSTPS); // Array of linearly interpolated appled mechanical forces
        if (apload == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = apload;
        nd++;

        double *pres = alloc_dbl (FNDOF*NTSTPS); // Array of linearly interpolated fluid pressures
        if (pres == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = pres;
        nd++;

        double *acc = alloc_dbl (FNDOF*NTSTPS); // Array of linearly interpolated fluid incident accelerations
        if (acc == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = acc;
        nd++;

        double *Lp = alloc_dbl (SNDOF*NTSTPS); // Array of linearly interpolated fluid incident accelerations
        if (Lp == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = Lp;
        nd++;

        double *Au = alloc_dbl (FNDOF*NTSTPS); // Array of linearly interpolated fluid incident accelerations
        if (Au == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = Au;
        nd++;

        // Pass control to q_fsi function
        q_fsi (jcode, qdyn, tstps, apload, pres, acc, L, A, Lp, Au, fdens, tinpt, pinpt, presinpt, accinpt);

        fscanf(IFP[0], "%lf,%lf\n", &alpha, &delta);

        // Pass control to stiff_fsi and mass_fsi functions
        stiff_fsi(minc, mcode, jcode, nnorm, tarea,farea, thick, deffarea, slength, defslen, L, A, ss, ss_fsi,
                  x, xlocal, emod, nu, Jinv, jac, yield, c1, c2, c3, ef, d, chi, efN, efM, maxa);
        mass_fsi (minc, mcode, jcode, nnorm, tarea, carea, farea, thick, slength, L, LT, sm, sm_fsi, x, dens, fdens, Jinv, jac);

        double ssd; // Dummy variable for solve function
        int det; // Flag for sign of determinant of tangent stiffness matrix

        // Variables our putput function
        double time = *(tstps);
        int dum = 0;

        // Pass control to output function
        output (&time, &dum, uc, ef, 0);

        // Pass control to solve function
        errchk = solve (jcode, ss, ss_fsi, sm, sm_fsi, sd_fsi, r, dd, maxa, &ssd, &det, um, vm, am, uc, vc, ac, qdyn, tstps,
                        Keff, Reff, Meff, alpha, delta, ipiv, 0, 1);
    }

    // Analysis for non-FSI
    if (ANAFLAG != 4) {

        if (ALGFLAG > 3) { // Dynamic analysis

            // Scan in the user desired number of time steps and total analysis
            fscanf(IFP[0], "%ld,%lf\n", &ntstpsinpt, &ttot);

            // Calculate dt
            dt = ttot/ntstpsinpt; ntstpsinpt += 1;
        }
        else {
            ntstpsinpt = 0;
        }

        // Allocate memory to arrays of input times and applied forces
        double *tinpt = alloc_dbl (ntstpsinpt); // Time array
        if (tinpt == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = tinpt;
        nd++;

        double *pinpt = alloc_dbl (NEQ*ntstpsinpt); // Applied mechanical forces acting on nodes
        if (pinpt == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = pinpt;
        nd++;

        double *dinpt = alloc_dbl (NEQ*ntstpsinpt); // Applied distributed loads acting on frame elements
        if (dinpt == NULL) {
            goto EXIT2;
        }
        p2p2d[nd] = dinpt;
        nd++;

        // Pass control to load function
        errchk = load (q, efFE_ref, x, llength, offset, osflag, c1, c2, c3, jnt, mcode, jcode, minc, tinpt, pinpt, dinpt, um, vm, am);

        // Terminate program if errors encountered
        if (errchk == 1) {
            goto EXIT2;
        }

        /*
         Define secondary non-array variables, common to both NR and MSAL algorithms
         */

        // Load proportionality parameters
        int det; // Flag for sign of determinant of tangent stiffness matrix
        long ptr; // Points to correct location in ef arrays
        double lpf, dlpf; // Current and incremental
        double dlpfp; // Incremental from previous load increment
        double dlpfmax, dlpfmin; // Maximum / minimum incremental
        double lpfmax; // Maximum allowable
        // Tolerances on displacement, force, and energy for convergence test
        double toldisp, tolforc, tolener;
        double intener1; // Internal energy from first equilibrium iteration
        int convchk; // Flag for convergence test
        int itecnt, itemax; // Iteration counter and maximum number of iterations
        int subcnt, submax; // Subdivision counter and maximum number of subdivisions
        int frcchk_fr, frcchk_sh; // Error check variable on forces functions


        if (ALGFLAG < 3) { // Static analysis
            // Initialize generalized total nodal displacement and internal force vectors
            for (i = 0; i < NEQ; ++i) {
                d[i] = 0;
                f[i] = 0;
            }

            // Initialize element force vectors
            for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                ef[i] = 0;
            }
            // Initialize truss deformed length variables
            for (i = 0; i < NE_TR; ++i) {
                defllen[i] = llength[i];
            }
            // Initialize frame element variables
            for (i = 0; i < NE_FR; ++i) {
                yldflag[i*2] = yldflag[i*2+1] = 0;
                defllen[NE_TR+i] = llength[NE_TR+i];
                for (j = 0; j < 14; ++j) {
                    efFE[i*14+j] = 0;
                }
            }
            // Initialize shell element variables
            for (i = 0; i < NE_SH; ++i) {
                deffarea[i] = farea[i];
                for (j = 0; j < 3; ++j) {
                    defslen[i*3+j] = slength[i*3+j];
                    chi[i*3+j] = 0;
                }
                for (j = 0; j < 9; ++j) {
                    efN[i*9+j] = 0;
                    efM[i*9+j] = 0;
                }
            }

            // Pass control to output function
            output (&lpf, &itecnt, d, ef, 0);

            if (ANAFLAG == 1) {

                // Read in solver parameters from input file
                fscanf(IFP[0], "%lf\n", &lpfmax);
                if (OPTFLAG == 2) {
                    fprintf(IFP[1], "%le\n", lpfmax);
                }

                /* Compute generalized total external load vector, accounting for
                 generalized fixed-end load vector */
                for (i = 0; i < NEQ; ++i) {
                    qtot[i] = q[i] * lpfmax;
                }

                // Initialize tangent stiffness matrix to zero
                for (i = 0; i < lss; ++i) {
                    ss[i] = 0;
                }

                if (NE_TR > 0) {
                    // Pass control to stiff_tr function
                    stiff_tr (ss, emod, carea, llength, defllen, yield, c1, c2, c3, ef, maxa,
                              mcode);
                }

                if (NE_FR > 0) {
                    // Pass control to stiff_fr function
                    stiff_fr (ss, emod, gmod, carea, offset, osflag, llength, defllen,
                              istrong, iweak, ipolar, iwarp, yldflag, yield, zstrong, zweak, c1,
                              c2, c3, ef, efFE, mendrel, maxa, mcode);
                }
                if (NE_SH > 0) {
                    // Pass control to stiff_sh function
                    stiff_sh (ss, emod, nu, x, xlocal, thick, farea, deffarea, slength,
                              defslen, yield, c1, c2, c3, ef, d, chi, efN, efM, maxa, minc, mcode);
                }

                if (NE_BR > 0) {
                    // Pass control to stiff_sh function
                    stiff_br (ss, x, emod, nu, minc, mcode, jcode, Jinv, jac);
                }

                double ssd;
                // Solve the system for incremental displacements
                //double ssd; // Dummy variable for solve function
                if (lss == 1) {
                    /* Carry out computation of incremental displacement directly for
                     lss = 1 */
                    d[0] = qtot[0] / ss[0];
                } else {

                    // Pass control to solve function
                    errchk = solve (jcode, ss, ss_fsi, sm, sm_fsi, sd_fsi, qtot, d, maxa, &ssd, &det, um, vm, am, uc, vc, ac, pinpt, tinpt, Keff, Reff, Meff, alpha, delta, ipiv, 0, 1);

                    // Terminate program if errors encountered
                    if (errchk == 1) {
                        goto EXIT2;
                    }
                }

                if (NE_TR > 0) {
                    // Pass control to forces_tr function
                    forces_tr (f, ef, d, emod, carea, llength, defllen, yield, c1, c2, c3,
                               mcode);
                }

                if (NE_FR > 0) {
                    // Pass control to forces_fr function
                    forces_fr (f, ef, ef, efFE_ref, efFE, efFE, yldflag, d, emod, gmod,
                               carea, offset, osflag, llength, defllen, istrong, iweak, ipolar,
                               iwarp, yield, zstrong, zweak, c1, c2, c3, c1, c2, c3, mendrel, mcode,
                               &dlpf, &itecnt);
                }

                if (NE_SH > 0) {
                    // Pass control to forces_sh function
                    forces_sh (f, ef, ef, efN, efM, d, d, chi, x, x, emod, nu, xlocal, thick,
                               farea, deffarea, slength, defslen, yield, c1, c2, c3, c1, c2, c3,
                               minc, mcode, jcode);
                }

                // Pass control to output function
                output (&lpfmax, &itecnt, d, ef, 1);

                fprintf(OFP[0], "\nSolution successful\n");

            } else {  // Nonliner analysis
                /*
                 Define secondary non-array variables, specific to NR / MNR algorithm
                 */
                double ssd; // Dummy variable for solve function
                int inccnt; // Load increment counter
                int solcnt, solmin; // Solution counter and minimum number of solutions

                // Read in solver parameters from input file
                fscanf(IFP[0], "%lf,%lf,%lf,%lf,%lf\n", &lpfmax, &lpf, &dlpf, &dlpfmax,
                       &dlpfmin);
                fscanf(IFP[0], "%d,%d,%d\n", &itemax, &submax, &solmin);
                fscanf(IFP[0], "%lf,%lf,%lf\n", &toldisp, &tolforc, &tolener);
                if (OPTFLAG == 2) {
                    fprintf(IFP[1], "%le,%le,%le,%le,%le\n", lpfmax,
                            lpf, dlpf, dlpfmax, dlpfmin);
                    fprintf(IFP[1], "%d,%d,%d\n", itemax, submax, solmin);
                    fprintf(IFP[1], "%lf,%lf,%lf\n", toldisp, tolforc, tolener);
                }

                // Initialize load step, converged solution, and subdivision counters
                inccnt = solcnt = subcnt = 0;
                /* Begin load incrementation; load will be incremented until load
                 proportionality factor is equal to user specified maximum */
                do {
                    // If load proportionality factor exceeds maximum, set equal to maximum
                    if (lpf > lpfmax) {
                        lpf = lpfmax;
                    }
                    /* Set all temporary variables and variables which refer to the structure
                     in its current configuration to values obtained at last successful
                     load increment; this step is required so as not to overwrite structure
                     properties prematurely if load increment is unsuccessful / invalid */
                    for (i = 0; i < NEQ; ++i) {
                        /* Compute generalized total external load vector, accounting for
                         generalized fixed-end load vector */
                        qtot[i] = q[i] * lpf;
                        /* Store generalized internal force vector from previous
                         configuration */
                        fp[i] = f[i];
                        d_temp[i] = d[i];
                        f_temp[i] = f[i];
                    }
                    dlpfp = dlpf;

                    // General
                    for (i = 0; i < NJ*3; ++i) {
                        x_temp[i] = x[i];
                    }
                    for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                        ef_i[i] = ef_ip[i] = ef[i];
                    }
                    for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                        c1_i[i] = c1_ip[i] = c1[i];
                        c2_i[i] = c2_ip[i] = c2[i];
                        c3_i[i] = c3_ip[i] = c3[i];
                    }
                    // Truss
                    for (i = 0; i < NE_TR; ++i) {
                        defllen_i[i] = defllen_ip[i] = defllen[i];
                    }
                    // Frame
                    for (i = 0; i < NE_FR; ++i) {
                        defllen_i[NE_TR+i] = defllen_ip[NE_TR+i] = defllen[NE_TR+i];
                        for (j = 0; j < 6; ++j) {
                            xfr_temp[i*6+j] = xfr[i*6+j];
                        }
                        for (j = 0; j < 14; ++j) {
                            efFE_i[i*14+j] = efFE_ip[i*14+j] = efFE[i*14+j];
                        }
                    }
                    // Shell
                    for (i = 0; i < NE_SH; ++i) {
                        deffarea_i[i] = deffarea_ip[i] = deffarea[i];
                        for (j = 0; j < 3; ++j) {
                            defslen_i[i*3+j] = defslen_ip[i*3+j] = defslen[i*3+j];
                            chi_temp[i*3+j] = chi[i*3+j];
                        }
                        for (j = 0; j < 9; ++j) {
                            efN_temp[i*9+j] = efN[i*9+j];
                            efM_temp[i*9+j] = efM[i*9+j];
                        }
                    }

                    // Re-initialize iteration counter at the start of each increment
                    itecnt = 0;

                    /* Start of each equilibrium iteration within load increment; iterations
                     will continue until convergence is reached or iteration count exceeds
                     user specified maximum */
                    frcchk_fr = frcchk_sh = 0;
                    do {
                        // Compute residual force vector
                        for (i = 0; i < NEQ; ++i) {
                            r[i] = qtot[i] - f_temp[i];
                        }

                        if (ALGFLAG == 1 || (ALGFLAG == 2 && itecnt == 0)) {
                            // Initialize tangent stiffness matrix to zero
                            for (i = 0; i < lss; ++i) {
                                ss[i] = 0;
                            }

                            if (NE_TR > 0) {
                                // Pass control to stiff_tr function
                                stiff_tr (ss, emod, carea, llength, defllen_ip, yield, c1_ip,
                                          c2_ip, c3_ip, ef_ip, maxa, mcode);
                            }
                            if (NE_FR > 0) {
                                // Pass control to stiff_fr function
                                stiff_fr (ss, emod, gmod, carea, offset, osflag, llength,
                                          defllen_ip, istrong, iweak, ipolar, iwarp, yldflag,
                                          yield, zstrong, zweak, c1_ip, c2_ip, c3_ip, ef_ip,
                                          efFE_ip, mendrel, maxa, mcode);
                            }
                            if (NE_SH > 0) {
                                // Pass control to stiff_sh function
                                stiff_sh (ss, emod, nu, x_temp, xlocal, thick, farea,
                                          deffarea_ip, slength, defslen_ip, yield, c1_ip, c2_ip,
                                          c3_ip, ef_ip, d_temp, chi_temp, efN_temp, efM_temp, maxa,
                                          minc, mcode);
                            }
                        }

                        // Solve the system for incremental displacements
                        if (lss == 1) {
                            /* Carry out computation of incremental displacement directly for
                             lss = 1 */
                            dd[0] = r[0] / ss[0];
                        } else {
                            if (ALGFLAG == 1 || (ALGFLAG == 2 && itecnt == 0)) {
                                // Pass control to solve function
                                errchk = solve (jcode, ss, ss_fsi, sm, sm_fsi, sd_fsi, r, dd, maxa, &ssd, &det, um, vm, am, uc, vc, ac, pinpt, tinpt,
                                                Keff, Reff, Meff, alpha, delta, ipiv, 0, 1);
                            } else {
                                // Pass control to solve function
                                errchk = solve (jcode, ss, ss_fsi, sm, sm_fsi, sd_fsi, r, dd, maxa, &ssd, &det, um, vm, am, uc, vc, ac, pinpt, tinpt,
                                                Keff, Reff, Meff, alpha, delta, ipiv, 1, 1);
                            }
                            // Terminate program if errors encountered
                            if (errchk == 1) {
                                goto EXIT2;
                            }
                        }

                        /* Update generalized total nodal displacement vector, store
                         generalized internal force vector from previous iteration, and
                         re-initialize generalized internal force vector */
                        for (i = 0; i < NEQ; ++i) {
                            d_temp[i] += dd[i];
                            f_ip[i] = f_temp[i];
                            f_temp[i] = 0;
                        }

                        // Pass control to updatc function
                        updatc (x_temp, x_ip, xfr_temp, dd, defllen_i, deffarea_i, defslen_i,
                                offset, osflag, auxpt, c1_i, c2_i, c3_i, minc, jcode);

                        if (NE_TR > 0) {
                            // Pass control to forces_tr function
                            forces_tr (f_temp, ef_i, d, emod, carea, llength, defllen_i,
                                       yield, c1_i, c2_i, c3_i, mcode);
                        }

                        if (NE_FR > 0) {
                            // Pass control to forces_fr function
                            frcchk_fr = forces_fr (f_temp, ef_ip, ef_i, efFE_ref, efFE_ip,
                                                   efFE_i, yldflag, dd, emod, gmod, carea, offset, osflag,
                                                   llength, defllen_ip, istrong, iweak, ipolar, iwarp, yield,
                                                   zstrong, zweak, c1_ip, c2_ip, c3_ip, c1_i, c2_i, c3_i,
                                                   mendrel, mcode, &dlpf, &itecnt);
                        }

                        if (NE_SH > 0) {
                            // Pass control to forces_sh function
                            frcchk_sh = forces_sh (f_temp, ef_ip, ef_i, efN_temp, efM_temp,
                                                   dd, d_temp, chi_temp, x_temp, x_ip, emod, nu, xlocal, thick,
                                                   farea, deffarea_ip, slength, defslen_ip, yield, c1_ip, c2_ip,
                                                   c3_ip, c1_i, c2_i, c3_i, minc, mcode, jcode);
                        }

                        // Update element internal forces from previous iteration
                        for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                            ef_ip[i] = ef_i[i];
                        }

                        if (itecnt == 0) {
                            // Compute internal energy from first iteration
                            intener1 = 0;
                            for (i = 0; i < NEQ; ++i) {
                                intener1 += dd[i] * (qtot[i] - fp[i]);
                            }
                        }


                        // Pass control to test function
                        errchk = test (d_temp, dd, f_temp, fp, qtot, f_ip, &intener1,
                                       &convchk, &toldisp, &tolforc, &tolener);

                        // Terminate program if errors encountered
                        if (errchk == 1) {
                            goto EXIT2;
                        }

                        // Update variables from previous iteration
                        // General
                        for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                            c1_ip[i] = c1_i[i];
                            c2_ip[i] = c2_i[i];
                            c3_ip[i] = c3_i[i];
                        }
                        // Truss
                        for (i = 0; i < NE_TR; ++i) {
                            defllen_ip[i] = defllen_i[i];
                        }
                        // Frame
                        for (i = 0; i < NE_FR; ++i) {
                            defllen_ip[NE_TR+i] = defllen_i[NE_TR+i];
                            for (j = 0; j < 14; ++j) {
                                efFE_ip[i*14+j] = efFE_i[i*14+j];
                            }
                        }
                        // Shell
                        for (i = 0; i < NE_SH; ++i) {
                            deffarea_ip[i] = deffarea_i[i];
                            for (j = 0; j < 3; ++j) {
                                defslen_ip[i*3+j] = defslen_i[i*3+j];
                            }
                        }
                        itecnt++; // Advance iteration counter
                    } while (convchk != 0 && frcchk_fr == 0 && frcchk_sh == 0 &&
                             itecnt <= itemax);

                    if (frcchk_fr == 2) {
                        dlpf = dlpfp; // Reset increment in load proportionality factor
                    } else if ((convchk != 0 || frcchk_fr != 0 || frcchk_sh != 0) &&
                               subcnt <= submax) {
                        if (lpf == lpfmax) {
                            fprintf(OFP[0], "\n***ERROR*** Maximum allowable load");
                            fprintf(OFP[0], " proportionality factor attempted without");
                            fprintf(OFP[0], " convergence\n");

                            goto EXIT2;
                        } else if (dlpfp == dlpfmin) {
                            fprintf(OFP[0], "\n***ERROR*** Minimum allowable increment of");
                            fprintf(OFP[0], " load proportionality factor reached\n");

                            goto EXIT2;
                        }

                        if (frcchk_fr != 1) {
                            // Decrease increment of load proportionality factor
                            dlpf = dlpfp / 2;
                        }

                        if (dlpf < dlpfmin) {
                            dlpf = dlpfmin;
                        }

                        // Step back load proportionality factor
                        lpf = lpf - dlpfp + dlpf;

                        subcnt++; // Advance subdivision counter
                        solcnt = 0; // Re-initialize converged solution counter
                    } else if (subcnt > submax) {
                        fprintf(OFP[0], "\n***ERROR*** Maximum allowable number of");
                        fprintf(OFP[0], " subdivisions exceeded\n");

                        goto EXIT2;
                    } else {
                        inccnt++; // Advance load increment counter

                        /* Update all permanent variables to values which represent structure
                         in its current configuration */
                        for (i = 0; i < NEQ; ++i) {
                            d[i] = d_temp[i];
                            f[i] = f_temp[i];
                        }
                        // General
                        for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                            ef[i] = ef_i[i];
                        }
                        // Frame
                        for (i = 0; i < NE_FR; ++i) {
                            for (j = 0; j < 14; ++j) {
                                efFE[i*14+j] = efFE_i[i*14+j];
                            }
                        }
                        // General
                        for (i = 0; i < NJ*3; ++i) {
                            x[i] = x_temp[i];
                        }
                        for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                            c1[i] = c1_i[i];
                            c2[i] = c2_i[i];
                            c3[i] = c3_i[i];
                        }
                        // Truss
                        for (i = 0; i < NE_TR; ++i) {
                            defllen[i] = defllen_i[i];
                        }
                        // Frame
                        for (i = 0; i < NE_FR; ++i) {
                            defllen[NE_TR+i] = defllen_i[NE_TR+i];
                            for (j = 0; j < 6; ++j) {
                                xfr[i*6+j] = xfr_temp[i*6+j];
                            }
                            if (yldflag[i*2] == 2) {
                                yldflag[i*2] = 0;
                            }
                            if (yldflag[i*2+1] == 2) {
                                yldflag[i*2+1] = 0;
                            }
                        }
                        // Shell
                        if (ANAFLAG == 2) {
                            for (i = 0; i < NE_SH; ++i) {
                                deffarea[i] = deffarea_i[i];
                                for (j = 0; j < 3; ++j) {
                                    defslen[i*3+j] = defslen_i[i*3+j];
                                }
                            }
                        } else {
                            for (i = 0; i < NE_SH; ++i) {
                                deffarea[i] = deffarea_i[i];
                                for (j = 0; j < 3; ++j) {
                                    defslen[i*3+j] = defslen_i[i*3+j];
                                    chi[i*3+j] = chi_temp[i*3+j];
                                }
                                for (j = 0; j < 9; ++j) {
                                    efN[i*9+j] = efN_temp[i*9+j];
                                    efM[i*9+j] = efM_temp[i*9+j];
                                }
                            }
                        }

                        // Pass control to output function
                        output (&lpf, &itecnt, d, ef, 1);

                        solcnt++;
                        subcnt = 0; // Re-initialize subdivision counter
                        /* If current load increment resulted in solmin converged solutions
                         in a row, increase increment in load proportionality factor */
                        if (solcnt >= solmin) {
                            dlpf *= 2; // Increase increment of load proportionality factor
                            if (dlpf > dlpfmax) {
                                dlpf = dlpfmax;
                            }
                            solcnt = 0; // Re-initialize solution counter
                        }
                        lpf += dlpf; // Increment load proportionality factor
                    }
                } while (lpf <= lpfmax);

                if (lpf >= lpfmax && convchk == 0) {
                    fprintf(OFP[0], "\nSolution successful\n");
                }
            }
        } else if (ALGFLAG == 3) {
            /*
             Define secondary variables which DO depend upon NEQ, specific to the MSAL
             algorithm
             */

            // Generalized total nodal displacement vector from previous load increment
            double *dp = alloc_dbl (NEQ);
            if (dp == NULL) {
                goto EXIT2;
            }
            p2p2d[nd] = dp;
            nd++;
            /* Generalized total nodal displacement vector from previously previous load
             increment */
            double *dpp = alloc_dbl (NEQ);
            if (dpp == NULL) {
                goto EXIT2;
            }
            p2p2d[nd] = dpp;
            nd++;
            // Miscellaneous incremental generalized nodal displacement vectors
            double *ddq = alloc_dbl (NEQ);
            if (ddq == NULL) {
                goto EXIT2;
            }
            p2p2d[nd] = ddq;
            nd++;
            double *ddr = alloc_dbl (NEQ);
            if (ddr == NULL) {
                goto EXIT2;
            }
            p2p2d[nd] = ddr;
            nd++;
            // Initial and current diagonals of tangent stiffness matrix
            double *ssd_o = alloc_dbl (NEQ);
            if (ssd_o == NULL) {
                goto EXIT2;
            }
            p2p2d[nd] = ssd_o;
            nd++;
            double *ssd = alloc_dbl (NEQ);
            if (ssd == NULL) {
                goto EXIT2;
            }
            p2p2d[nd] = ssd;
            nd++;

            /*
             Define secondary non-array variables, specific to MSAL algorithm
             */
            // Displacements at DOF k
            double dk; // Prescribed initial
            double dkc; // Current (absolute value)
            double dkimax; // Maximum allowable
            long dkdof; // DOF affected by prescribed initial displacement
            // MSAL parameters
            double arc, beta; // Arc length and arc length adjustment factor
            double a, b, c; // Coefficients of quadratic equation
            double alpha, temp, psi, psi_thresh; // Miscellaneous factors
            int iteopt; // Optimum number of eqilibrium iterations
            double dnorm, dnormallow; // Generic incremental displacement norms
            double dotprod; // Generic dot product
            // Load proportionality parameters
            double lpfc; // Current (absolute value)
            double lpf_temp;
            double lpfp, lpfpp; // Previous and previously previous
            // Subdivisions due to imaginary roots in arc length criterion
            int imagcnt, imagmax;
            // Subdivisions due to two neg. roots in arc length criterion
            int negcnt, negmax;
            int errchk2; // Error check on quad function

            // Pass control to msal function
            errchk = msal (&dk, &dkdof, jnt, jcode);

            // Terminate program if errors encountered
            if (errchk == 1) {
                goto EXIT2;
            }

            // Read in MSAL parameters from input file
            fscanf(IFP[0], "%lf\n", &alpha);
            fscanf(IFP[0], "%lf\n", &psi_thresh);
            fscanf(IFP[0], "%d\n", &iteopt);
            fscanf(IFP[0], "%lf,%lf\n", &lpfmax, &dkimax);
            fscanf(IFP[0], "%d,%d,%d,%d\n", &itemax, &submax, &imagmax, &negmax);
            fscanf(IFP[0], "%lf,%lf,%lf\n", &toldisp, &tolforc, &tolener);
            if (OPTFLAG == 2) {
                fprintf(IFP[1], "%le\n", alpha);
                fprintf(IFP[1], "%le\n", psi_thresh);
                fprintf(IFP[1], "%d\n", iteopt);
                fprintf(IFP[1], "%le,%le\n", lpfmax, dkimax);
                fprintf(IFP[1], "%d,%d,%d,%d\n", itemax, submax, imagmax, negmax);
                fprintf(IFP[1], "%le,%le,%le\n", toldisp, tolforc, tolener);
            }

            // Initialize generalized total nodal displacement and internal force vectors
            for (i = 0; i < NEQ; ++i) {
                d[i] = dp[i] = 0;
                f[i] = fp[i] = 0;
            }
            lpfp = 0;

            // Initialize element force vectors
            for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                ef_i[i] = ef_ip[i] = ef[i] = 0;
            }
            // Initialize direction cosines
            for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                c1_i[i] = c1_ip[i] = c1[i];
                c2_i[i] = c2_ip[i] = c2[i];
                c3_i[i] = c3_ip[i] = c3[i];
            }
            // Initialize truss deformed length variables
            for (i = 0; i < NE_TR; ++i) {
                defllen_i[i] = defllen_ip[i] = defllen[i] = llength[i];
            }
            // Initialize frame element variables
            for (i = 0; i < NE_FR; ++i) {
                yldflag[i*2] = yldflag[i*2+1] = 0;
                defllen_i[NE_TR+i] = defllen_ip[NE_TR+i] = defllen[NE_TR+i] =
                llength[NE_TR+i];
                for (j = 0; j < 14; ++j) {
                    efFE_i[i*14+j] = efFE_ip[i*14+j] = efFE[i*14+j] = 0;
                }
            }
            // Initialize shell element variables
            for (i = 0; i < NE_SH; ++i) {
                deffarea_i[i] = deffarea_ip[i] = deffarea[i] = farea[i];
                for (j = 0; j < 3; ++j) {
                    defslen_i[i*3+j] = defslen_ip[i*3+j] = defslen[i*3+j] = slength[i*3+j];
                }
            }

            /* Begin load incrementation; follows Bathe and Dvorkin's algorithm for the
             initial load increment described in Bathe and Dvorkin (1981) */

            // Initialize the iteration counter to first iteration of first load increment
            itecnt = 0;

            // Initialize tangent stiffness matrix to zero
            for (i = 0; i < lss; ++i) {
                ss[i] = 0;
            }

            if (NE_TR > 0) {
                // Pass control to stiff_tr function
                stiff_tr (ss, emod, carea, llength, defllen, yield, c1, c2, c3, ef, maxa,
                          mcode);
            }

            if (NE_FR > 0) {
                // Pass control to stiff_fr function
                stiff_fr (ss, emod, gmod, carea, offset, osflag, llength, defllen, istrong,
                          iweak, ipolar, iwarp, yldflag, yield, zstrong, zweak, c1, c2, c3, ef,
                          efFE, mendrel, maxa, mcode);
            }

            if (NE_SH > 0) {
                // Pass control to stiff_sh function
                stiff_sh (ss, emod, nu, x, xlocal, thick, farea, deffarea, slength, defslen,
                          yield, c1, c2, c3, ef, d, chi, efN, efM, maxa, minc, mcode);
            }

            // Solve the system for incremental displacements
            if (lss == 1) {
                /* Carry out computation of incremental displacement directly for
                 lss = 1 */
                ddq[0] = q[0] / ss[0];
                ssd[0] = ss[0];
            } else {
                // Pass control to solve function
                errchk = solve (jcode, ss, ss_fsi, sm, sm_fsi, sd_fsi, q, ddq, maxa, ssd, &det, um, vm, am, uc, vc, ac, pinpt, tinpt,
                                Keff, Reff, Meff, alpha, delta, ipiv, 0, 1);

                // Terminate program if errors encountered
                if (errchk == 1) {
                    goto EXIT2;
                }
            }

            // Compute load proportionality factor for first iteration
            lpf = dk / ddq[dkdof];

            intener1 = 0;
            for (i = 0; i < NEQ; ++i) {
                // Compute generalized incremental and total nodal displacement vectors
                d[i] = dd[i] = lpf * ddq[i];
                // Compute internal energy from the first iteration
                intener1 += dd[i] * (lpf * q[i]);
                // Store initial diagonals of tangent stiffness matrix
                ssd_o[i] = ssd[i];
            }

            // Pass control to updatc function
            updatc (x, x_ip, xfr, dd, defllen_i, deffarea_i, defslen_i, offset, osflag,
                    auxpt, c1_i, c2_i, c3_i, minc, jcode);

            if (NE_TR > 0) {
                // Pass control to forces_tr function
                forces_tr (f, ef_i, d, emod, carea, llength, defllen_i, yield, c1_i, c2_i,
                           c3_i, mcode);
            }

            if (NE_FR > 0) {
                // Pass control to forces_fr function
                frcchk_fr = forces_fr (f, ef_ip, ef_i, efFE_ref, efFE_ip, efFE_i, yldflag,
                                       dd, emod, gmod, carea, offset, osflag, llength, defllen_ip, istrong,
                                       iweak, ipolar, iwarp, yield, zstrong, zweak, c1_ip, c2_ip, c3_ip, c1_i,
                                       c2_i, c3_i, mendrel, mcode, &lpf, &itecnt);
            }

            if (NE_SH > 0) {
                // Pass control to forces_sh function
                frcchk_sh = forces_sh (f, ef_ip, ef_i, efN, efM, dd, d, chi, x, x_ip, emod,
                                       nu, xlocal, thick, farea, deffarea_ip, slength, defslen_ip, yield, c1_ip,
                                       c2_ip, c3_ip, c1_i, c2_i, c3_i, minc, mcode, jcode);
            }



            itecnt = 1; // Advance iteration counter

            // Update variables from previous iteration
            // General
            for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                ef_ip[i] = ef_i[i];
            }
            for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                c1_ip[i] = c1_i[i];
                c2_ip[i] = c2_i[i];
                c3_ip[i] = c3_i[i];
            }
            // Truss
            for (i = 0; i < NE_TR; ++i) {
                defllen_ip[i] = defllen_i[i];
            }
            // Frame
            for (i = 0; i < NE_FR; ++i) {
                defllen_ip[NE_TR+i] = defllen_i[NE_TR+i];
                for (j = 0; j < 14; ++j) {
                    efFE_ip[i*14+j] = efFE_i[i*14+j];
                }
            }
            // Shell
            for (i = 0; i < NE_SH; ++i) {
                deffarea_ip[i] = deffarea_i[i];
                for (j = 0; j < 3; ++j) {
                    defslen_ip[i*3+j] = defslen_i[i*3+j];
                }
            }

            frcchk_fr = frcchk_sh = 0;
            do {
                for (i = 0; i < NEQ; ++i) {
                    /* Compute generalized total external load vector, accounting for
                     generalized fixed-end load vector */
                    qtot[i] = q[i] * lpf;
                    // Compute residual force vector
                    r[i] = qtot[i] - f[i];
                }

                // Solve the system for incremental displacements
                if (lss == 1) {
                    /* Carry out computation of incremental displacement directly for
                     lss = 1 */
                    ddr[0] = r[0] / ss[0];
                    ssd[0] = ss[0];
                } else {
                    // Pass control to solve function
                    errchk = solve (jcode, ss, ss_fsi, sm, sm_fsi, sd_fsi, r, ddr, maxa, ssd, &det, um, vm, am, uc, vc, ac, pinpt, tinpt,
                                    Keff, Reff, Meff, alpha, delta, ipiv, 1, 1);

                    // Terminate program if errors encountered
                    if (errchk == 1) {
                        goto EXIT2;
                    }
                }

                // Compute increment in load proportionality factor, add to total
                dlpf = -ddr[dkdof] / ddq[dkdof];
                lpf += dlpf;

                /* Compute generalized incremental nodal displacement vector, update
                 generalized total nodal displacement vector, store generalized internal
                 force vector from previous iteration, and re-initialize generalized
                 internal force vector */
                for (i = 0; i < NEQ; ++i) {
                    dd[i] = ddr[i] + dlpf * ddq[i];
                    d[i] += dd[i];
                    f_ip[i] = f[i];
                    f[i] = 0;
                }

                // Pass control to updatc function
                updatc (x, x_ip, xfr, dd, defllen_i, deffarea_i, defslen_i, offset, osflag,
                        auxpt, c1_i, c2_i, c3_i, minc, jcode);

                if (NE_TR > 0) {
                    // Pass control to forces_tr function
                    forces_tr (f, ef_i, d, emod, carea, llength, defllen_i, yield, c1_i,
                               c2_i, c3_i, mcode);
                }

                if (NE_FR > 0) {
                    // Pass control to forces_fr function
                    frcchk_fr = forces_fr (f, ef_ip, ef_i, efFE_ref, efFE_ip, efFE_i,
                                           yldflag, dd, emod, gmod, carea, offset, osflag, llength, defllen_ip,
                                           istrong, iweak, ipolar, iwarp, yield, zstrong, zweak, c1_ip, c2_ip,
                                           c3_ip, c1_i, c2_i, c3_i, mendrel, mcode, &dlpf, &itecnt);
                }

                if (NE_SH > 0) {
                    // Pass control to forces_sh function
                    frcchk_sh = forces_sh (f, ef_ip, ef_i, efN, efM, dd, d, chi, x, x_ip,
                                           emod, nu, xlocal, thick, farea, deffarea_ip, slength, defslen_ip,
                                           yield, c1_ip, c2_ip, c3_ip, c1_i, c2_i, c3_i, minc, mcode, jcode);
                }

                // Pass control to test function
                errchk = test (d, dd, f, fp, qtot, f_ip, &intener1, &convchk, &toldisp,
                               &tolforc, &tolener);

                // Terminate program if errors encountered
                if (errchk == 1) {
                    goto EXIT2;
                }

                itecnt++; // Advance iteration counter

                // Update variables from previous iteration
                // General
                for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                    ef_ip[i] = ef_i[i];
                }
                for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                    c1_ip[i] = c1_i[i];
                    c2_ip[i] = c2_i[i];
                    c3_ip[i] = c3_i[i];
                }
                // Truss
                for (i = 0; i < NE_TR; ++i) {
                    defllen_ip[i] = defllen_i[i];
                }
                // Frame
                for (i = 0; i < NE_FR; ++i) {
                    defllen_ip[NE_TR+i] = defllen_i[NE_TR+i];
                    for (j = 0; j < 14; ++j) {
                        efFE_ip[i*14+j] = efFE_i[i*14+j];
                    }
                }
                // Shell
                for (i = 0; i < NE_SH; ++i) {
                    deffarea_ip[i] = deffarea_i[i];
                    for (j = 0; j < 3; ++j) {
                        defslen_ip[i*3+j] = defslen_i[i*3+j];
                    }
                }
            } while (convchk != 0 && frcchk_fr == 0 && frcchk_sh == 0 && itecnt <= itemax);

            if (convchk != 0 || frcchk_fr != 0 || frcchk_sh != 0) {
                fprintf(OFP[0], "\n***ERROR*** Initially presecribed displacement too");
                fprintf(OFP[0], " large\n");

                goto EXIT2;
            } else {
                /* Update all permanent variables to values which represent structure in its
                 current configuration */
                // General
                for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                    ef[i] = ef_i[i];
                }
                for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                    c1[i] = c1_i[i];
                    c2[i] = c2_i[i];
                    c3[i] = c3_i[i];
                }
                // Truss
                for (i = 0; i < NE_TR; ++i) {
                    defllen[i] = defllen_i[i];
                }
                // Frame
                for (i = 0; i < NE_FR; ++i) {
                    defllen[NE_TR+i] = defllen_i[NE_TR+i];
                    for (j = 0; j < 14; ++j) {
                        efFE[i*14+j] = efFE_i[i*14+j];
                    }
                }
                // Shell
                for (i = 0; i < NE_SH; ++i) {
                    deffarea[i] = deffarea_i[i];
                    for (j = 0; j < 3; ++j) {
                        defslen[i*3+j] = defslen_i[i*3+j];
                    }
                }

                // Pass control to output function
                output (&lpf, &itecnt, d, ef, 0);

                // Pass control to output function
                output (&lpf, &itecnt, d, ef, 1);
            }

            /* Store displacement at DOF "k" and current load proportionality factor for
             comparison with maximum allowable values */
            dkc = fabs(d[dkdof]);
            lpfc = fabs(lpf);

            /* Continue load incrementation employing Bathe and Dvorkin's arc length
             algorithm described in Bathe and Dvorkin (1981) with the addition of the "psi"
             term from Crisfield and Shi (1991) which eliminates load control in the arc
             length criterion; termed here as modified spherical arc length (MSAL) solution
             algorithm */

            // Store all variables from previously converged load increments
            for (i = 0; i < NEQ; ++i) {
                dpp[i] = dp[i];
                dp[i] = d[i];
                fp[i] = f[i];
            }
            lpfpp = lpfp;
            lpfp = lpf;

            /* Set all temporary variables, and variables which refer to the structure in its
             current configuration, to values obtained at last successful load increment;
             this step is required so as not to overwrite structure properties prematurely
             if load increment is unsuccessful / invalid */
            // General
            for (i = 0; i < NJ*3; ++i) {
                x_temp[i] = x[i];
            }
            for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                ef_i[i] = ef_ip[i] = ef[i];
            }
            for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                c1_i[i] = c1_ip[i] = c1[i];
                c2_i[i] = c2_ip[i] = c2[i];
                c3_i[i] = c3_ip[i] = c3[i];
            }
            // Truss
            for (i = 0; i < NE_TR; ++i) {
                defllen_i[i] = defllen_ip[i] = defllen[i];
            }
            // Frame
            for (i = 0; i < NE_FR; ++i) {
                defllen_i[NE_TR+i] = defllen_ip[NE_TR+i] = defllen[NE_TR+i];
                for (j = 0; j < 6; ++j) {
                    xfr_temp[i*6+j] = xfr[i*6+j];
                }
                for (j = 0; j < 14; ++j) {
                    efFE_i[i*14+j] = efFE_ip[i*14+j] = efFE[i*14+j];
                }
            }
            // Shell
            for (i = 0; i < NE_SH; ++i) {
                deffarea_i[i] = deffarea_ip[i] = deffarea[i];
                for (j = 0; j < 3; ++j) {
                    defslen_i[i*3+j] = defslen_ip[i*3+j] = defslen[i*3+j];
                    chi_temp[i*3+j] = chi[i*3+j];
                }
                for (j = 0; j < 9; ++j) {
                    efN_temp[i*9+j] = efN[i*9+j];
                    efM_temp[i*9+j] = efM[i*9+j];
                }
            }

            /* Compute arc length adjustment factor, beta, from Euclidean norm of current
             total displacement vector and allowable displacement parameter, i.e.
             alpha * (Euclidean norm of displacement vector from first load step) */
            dnorm = sqrt(dot (dp, dp, NEQ));
            dnormallow = alpha * sqrt(dot (dp, dp, NEQ));
            beta = sqrt(((double) iteopt) / (double) itecnt) * (dnormallow / dnorm);

            /* Compute factor on load control, psi, to eliminate load control in arc length
             criterion when in the vicinity of a critical point; this corresponds with the
             the ratio of initial to current diagonal members of the tangent stiffness
             matrix */
            psi = 1;
            for (i = 0; i < NEQ; ++i) {
                temp = fabs(ssd[i] / ssd_o[i]);
                if (temp < psi) {
                    psi = temp;
                }
            }

            /* Begin load incrementation; initialize errchk2 and all counter variables to
             zero */
            errchk2 = subcnt = imagcnt = negcnt = 0;
            while (lpfc <= lpfmax && dkc <= dkimax) {
                /* Store generalized total nodal displacement vector and load proportionality
                 factor from current configuration */
                for (i = 0; i < NEQ; ++i) {
                    d_temp[i] = d[i];
                }
                lpf_temp = lpf;

                if (errchk2 == 0) {
                    // Compute arc length
                    dotprod = 0;
                    for (i = 0; i < NEQ; ++i) {
                        dotprod += (dp[i] - dpp[i]) * (dp[i] - dpp[i]);
                    }
                    if (psi >= psi_thresh) {
                        arc = beta * sqrt(dotprod + (lpfp - lpfpp) * (lpfp - lpfpp));
                    } else {
                        arc = beta * sqrt(dotprod);
                    }
                } else {
                    /* Reset all temporary variables and variables which refer to the
                     structure in its current configuration to values obtained at last
                     successful load increment */
                    // General
                    for (i = 0; i < NJ*3; ++i) {
                        x_temp[i] = x[i];
                    }
                    for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                        ef_i[i] = ef_ip[i] = ef[i];
                    }
                    for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                        c1_i[i] = c1_ip[i] = c1[i];
                        c2_i[i] = c2_ip[i] = c2[i];
                        c3_i[i] = c3_ip[i] = c3[i];
                    }
                    // Truss
                    for (i = 0; i < NE_TR; ++i) {
                        defllen_i[i] = defllen_ip[i] = defllen[i];
                    }
                    // Frame
                    for (i = 0; i < NE_FR; ++i) {
                        defllen_i[NE_TR+i] = defllen_ip[NE_TR+i] = defllen[NE_TR+i];
                        for (j = 0; j < 6; ++j) {
                            xfr_temp[i*6+j] = xfr[i*6+j];
                        }
                        for (j = 0; j < 14; ++j) {
                            efFE_i[i*14+j] = efFE_ip[i*14+j] = efFE[i*14+j];
                        }
                    }
                    // Shell
                    for (i = 0; i < NE_SH; ++i) {
                        deffarea_i[i] = deffarea_ip[i] = deffarea[i];
                        for (j = 0; j < 3; ++j) {
                            defslen_i[i*3+j] = defslen_ip[i*3+j] = defslen[i*3+j];
                            chi_temp[i*3+j] = chi[i*3+j];
                        }
                        for (j = 0; j < 9; ++j) {
                            efN_temp[i*9+j] = efN[i*9+j];
                            efM_temp[i*9+j] = efM[i*9+j];
                        }
                    }

                    errchk2 = 0; // Re-initialize at start of each increment
                }

                // Initialize tangent stiffness matrix to zero
                for (i = 0; i < lss; ++i) {
                    ss[i] = 0;
                }

                if (NE_TR > 0) {
                    // Pass control to stiff_tr function
                    stiff_tr (ss, emod, carea, llength, defllen, yield, c1, c2, c3, ef, maxa,
                              mcode);
                }

                if (NE_FR > 0) {
                    // Pass control to stiff_fr function
                    stiff_fr (ss, emod, gmod, carea, offset, osflag, llength, defllen,
                              istrong, iweak, ipolar, iwarp, yldflag, yield, zstrong, zweak, c1,
                              c2, c3, ef, efFE, mendrel, maxa, mcode);
                }

                if (NE_SH > 0) {
                    // Pass control to stiff_sh function
                    stiff_sh (ss, emod, nu, x, xlocal, thick, farea, deffarea, slength,
                              defslen, yield, c1, c2, c3, ef, d, chi, efN, efM, maxa, minc, mcode);
                }

                // Solve the system for incremental displacements
                if (lss == 1) {
                    // Carry out computation of incremental displacement directly for lss = 1
                    ddq[0] = q[0] / ss[0];
                    ssd[0] = ss[0];
                    if (ss[0] > 0) {
                        det = 0;
                    } else {
                        det = 1;
                    }
                } else {
                    // Pass control to solve function
                    errchk = solve (jcode, ss, ss_fsi, sm, sm_fsi, sd_fsi, q, ddq, maxa, ssd, &det, um, vm, am, uc, vc, ac, pinpt, tinpt,
                                    Keff, Reff, Meff, alpha, delta, ipiv, 0, 1);

                    // Terminate program if errors encountered
                    if (errchk == 1) {
                        goto EXIT2;
                    }
                }

                // Compute a-coefficient of quadratic equation for solution of dlpf
                if (psi >= psi_thresh) {
                    a = dot(q, q, NEQ) + dot (ddq, ddq, NEQ);
                } else {
                    a = dot (ddq, ddq, NEQ);
                }

                // Compute load proportionality factor for first iteration
                if (det == 0) {
                    dlpf = arc * sqrt(1 / a);
                } else {
                    dlpf = -arc * sqrt(1 / a);
                }
                lpf_temp += dlpf;

                intener1 = 0;
                for (i = 0; i < NEQ; ++i) {
                    /* Compute generalized incremental nodal displacement vector and update
                     generalized total nodal displacement vector */
                    dd[i] = dlpf * ddq[i];
                    d_temp[i] += dd[i];
                    // Compute internal energy from the first iteration
                    intener1 += dd[i] * (dlpf * q[i]);
                    // Re-initialize generalized internal force vector
                    f_temp[i] = 0;
                }

                // Pass control to updatc function
                updatc (x_temp, x_ip, xfr_temp, dd, defllen_i, deffarea_i, defslen_i, offset,
                        osflag, auxpt, c1_i, c2_i, c3_i, minc, jcode);

                if (NE_TR > 0) {
                    // Pass control to forces_tr function
                    forces_tr (f_temp, ef_i, d, emod, carea, llength, defllen_i, yield, c1_i,
                               c2_i, c3_i, mcode);
                }

                if (NE_FR > 0) {
                    // Pass control to forces_fr function
                    forces_fr (f_temp, ef_ip, ef_i, efFE_ref, efFE_ip, efFE_i, yldflag, dd,
                               emod, gmod, carea, offset, osflag, llength, defllen_ip, istrong,
                               iweak, ipolar, iwarp, yield, zstrong, zweak, c1_ip, c2_ip, c3_ip,
                               c1_i, c2_i, c3_i, mendrel, mcode, &dlpf, &itecnt);
                }

                if (NE_SH > 0) {
                    // Pass control to forces_sh function
                    forces_sh (f_temp, ef_ip, ef_i, efN_temp, efM_temp, dd, d_temp, chi_temp,
                               x_temp, x_ip, emod, nu, xlocal, thick, farea, deffarea_ip, slength,
                               defslen_ip, yield, c1_ip, c2_ip, c3_ip, c1_i, c2_i, c3_i, minc,
                               mcode, jcode);
                }

                // Update variables from previous iteration
                // General
                for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                    ef_ip[i] = ef_i[i];
                }
                for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                    c1_ip[i] = c1_i[i];
                    c2_ip[i] = c2_i[i];
                    c3_ip[i] = c3_i[i];
                }
                // Truss
                for (i = 0; i < NE_TR; ++i) {
                    defllen_ip[i] = defllen_i[i];
                }
                // Frame
                for (i = 0; i < NE_FR; ++i) {
                    defllen_ip[NE_TR+i] = defllen_i[NE_TR+i];
                    for (j = 0; j < 14; ++j) {
                        efFE_ip[i*14+j] = efFE_i[i*14+j];
                    }
                }
                // Shell
                for (i = 0; i < NE_SH; ++i) {
                    deffarea_ip[i] = deffarea_i[i];
                    for (j = 0; j < 3; ++j) {
                        defslen_ip[i*3+j] = defslen_i[i*3+j];
                    }
                }

                /* Initialize the iteration counter to second iteration of current load
                 increment */
                itecnt = 1;

                frcchk_fr = frcchk_sh = 0;
                do {
                    for (i = 0; i < NEQ; ++i) {
                        /* Compute generalized total external load vector, accounting for
                         generalized fixed-end load vector */
                        qtot[i] = q[i] * lpf_temp;
                        // Compute residual force vector
                        r[i] = qtot[i] - f_temp[i];
                        /* Store generalized internal force vector from previous iteration
                         and re-initialize generalized internal force vector */
                        f_ip[i] = f_temp[i];
                        f_temp[i] = 0;
                    }

                    // Solve the system for incremental displacements
                    if (lss == 1) {
                        /* Carry out computation of incremental displacement directly for
                         lss = 1 */
                        ddr[0] = r[0] / ss[0];
                    } else {
                        // Pass control to solve function
                        errchk = solve (jcode, ss, ss_fsi, sm, sm_fsi, sd_fsi, r, ddr, maxa, ssd, &det, um, vm, am, uc, vc, ac, pinpt, tinpt,
                                        Keff, Reff, Meff, alpha, delta, ipiv, 1, 1);

                        // Terminate program if errors encountered
                        if (errchk == 1) {
                            goto EXIT2;
                        }
                    }

                    /* Compute b- and c-coefficients of quadratic equation for solution of
                     dlpf */
                    if (psi >= psi_thresh) {
                        b = 2 * (dot (d_temp, ddq, NEQ) - dot (dp, ddq, NEQ) +
                                 dot (ddr, ddq, NEQ) + (lpf_temp - lpfp) * dot(q, q, NEQ));
                        c = 2 * (dot (d_temp, ddr, NEQ) - dot (dp, ddr, NEQ) -
                                 dot (d_temp, dp, NEQ)) + dot (d_temp, d_temp, NEQ) +
                        dot (dp, dp, NEQ) + dot (ddr, ddr, NEQ) +
                        (lpf_temp - lpfp) * (lpf_temp - lpfp) * dot(q, q, NEQ) -
                        arc * arc;
                    } else {
                        b = 2 * (dot (d_temp, ddq, NEQ) - dot (dp, ddq, NEQ) +
                                 dot (ddr, ddq, NEQ));
                        c = 2 * (dot (d_temp, ddr, NEQ) - dot (dp, ddr, NEQ) -
                                 dot (d_temp, dp, NEQ)) + dot (d_temp, d_temp, NEQ) +
                        dot (dp, dp, NEQ) + dot (ddr, ddr, NEQ) - arc * arc;
                    }

                    // Pass control to quad function
                    errchk2 = quad (&a, &b, &c, d_temp, dp, ddr, ddq, dd, &dlpf, &lpf_temp);

                    if (errchk2 == 0) {
                        // Pass control to updatc function
                        updatc (x_temp, x_ip, xfr_temp, dd, defllen_i, deffarea_i, defslen_i,
                                offset, osflag, auxpt, c1_i, c2_i, c3_i, minc, jcode);

                        if (NE_TR > 0) {
                            // Pass control to forces_tr function
                            forces_tr (f_temp, ef_i, d, emod, carea, llength, defllen_i,
                                       yield, c1_i, c2_i, c3_i, mcode);
                        }

                        if (NE_FR > 0) {
                            // Pass control to forces_fr function
                            frcchk_fr = forces_fr (f_temp, ef_ip, ef_i, efFE_ref, efFE_ip,
                                                   efFE_i, yldflag, dd, emod, gmod, carea, offset, osflag,
                                                   llength, defllen_ip, istrong, iweak, ipolar, iwarp, yield,
                                                   zstrong, zweak, c1_ip, c2_ip, c3_ip, c1_i, c2_i, c3_i,
                                                   mendrel, mcode, &dlpf, &itecnt);
                        }

                        if (NE_SH > 0) {
                            // Pass control to forces_sh function
                            frcchk_sh = forces_sh (f_temp, ef_ip, ef_i, efN_temp, efM_temp,
                                                   dd, d_temp, chi_temp, x_temp, x_ip, emod, nu, xlocal, thick,
                                                   farea, deffarea_ip, slength, defslen_ip, yield, c1_ip, c2_ip,
                                                   c3_ip, c1_i, c2_i, c3_i, minc, mcode, jcode);
                        }

                        // Pass control to test function
                        errchk = test (d_temp, dd, f_temp, fp, qtot, f_ip, &intener1,
                                       &convchk, &toldisp, &tolforc, &tolener);

                        // Terminate program if errors encountered
                        if (errchk == 1) {
                            goto EXIT2;
                        }

                        if (convchk == 0) {
                            /* Compute norm of displacement increment for comparison against
                             allowable */
                            dnorm = 0;
                            for (i = 0; i < NEQ; ++i) {
                                dnorm += (d_temp[i] - dp[i]) * (d_temp[i] - dp[i]);
                            }
                            dnorm = sqrt(dnorm);

                            if (dnorm > 100 * dnormallow) {
                                /* If incremental displacement norm exceeds allowable, reduce
                                 arc length and re-attempt solution */
                                arc /= beta;
                                beta = dnormallow / dnorm;
                                arc *= beta;

                                // Re-initialize generalized internal force vector
                                for (i = 0; i < NEQ; ++i) {
                                    f_temp[i] = f[i];
                                }

                                errchk2 = 1;
                            }
                        } else if (frcchk_fr == 2) {
                            // Re-initialize generalized internal force vector
                            for (i = 0; i < NEQ; ++i) {
                                f_temp[i] = f[i];
                            }

                            errchk2 = 1;
                        } else if ((frcchk_fr != 0 || frcchk_sh != 0) && subcnt <= submax) {
                            arc *= 0.5; // Reduce arc length and re-attempt load increment

                            // Re-initialize generalized internal force vector
                            for (i = 0; i < NEQ; ++i) {
                                f_temp[i] = f[i];
                            }

                            subcnt++;
                            errchk2 = 1;
                        } else {
                            itecnt++; // Advance solution counter

                            if (itecnt > itemax) {
                                // Reduce arc length and re-attempt load increment
                                arc *= 0.5;

                                // Re-initialize generalized internal force vector
                                for (i = 0; i < NEQ; ++i) {
                                    f_temp[i] = f[i];
                                }

                                subcnt++;
                                errchk2 = 1;
                            }

                            // Update variables from previous iteration
                            // General
                            for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                                ef_ip[i] = ef_i[i];
                            }
                            for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                                c1_ip[i] = c1_i[i];
                                c2_ip[i] = c2_i[i];
                                c3_ip[i] = c3_i[i];
                            }
                            // Truss
                            for (i = 0; i < NE_TR; ++i) {
                                defllen_ip[i] = defllen_i[i];
                            }
                            // Frame
                            for (i = 0; i < NE_FR; ++i) {
                                defllen_ip[NE_TR+i] = defllen_i[NE_TR+i];
                                for (j = 0; j < 14; ++j) {
                                    efFE_ip[i*14+j] = efFE_i[i*14+j];
                                }
                            }
                            // Shell
                            for (i = 0; i < NE_SH; ++i) {
                                deffarea_ip[i] = deffarea_i[i];
                                for (j = 0; j < 3; ++j) {
                                    defslen_ip[i*3+j] = defslen_i[i*3+j];
                                }
                            }
                        }
                    } else if (errchk2 == 2) {
                        arc *= 0.5; // Reduce arc length and re-attempt load increment

                        // Re-initialize generalized internal force vector
                        for (i = 0; i < NEQ; ++i) {
                            f_temp[i] = f[i];
                        }

                        imagcnt++;
                    } else if (errchk2 == 3) {
                        arc *= 0.5; // Reduce arc length and re-attempt load increment

                        // Re-initialize generalized internal force vector
                        for (i = 0; i < NEQ; ++i) {
                            f_temp[i] = f[i];
                        }

                        negcnt++;
                    } else if (errchk2 == 4) {
                        goto EXIT2;
                    }
                } while (convchk != 0 && errchk2 == 0 && subcnt <= submax &&
                         imagcnt <= imagmax && negcnt <= negmax);

                if (subcnt > submax || imagcnt > imagmax || negcnt > negmax) {
                    fprintf(OFP[0], "\n***ERROR*** Maximum allowable number of");
                    fprintf(OFP[0], " subdivisions exceeded\n");

                    goto EXIT2;
                } else if (errchk2 == 0) {
                    /* If incremental displacement norm less than allowable, increase arc
                     length */
                    beta = sqrt(((double) iteopt) / ((double) itecnt)) *
                    (dnormallow / dnorm);

                    /* Update all permanent variables to values which represent structure in
                     its current configuration */
                    for (i = 0; i < NEQ; ++i) {
                        dpp[i] = dp[i];
                        d[i] = dp[i] = d_temp[i];
                        f[i] = fp[i] = f_temp[i];
                    }
                    lpfpp = lpfp;
                    lpf = lpfp = lpf_temp;

                    // General
                    for (i = 0; i < NJ*3; ++i) {
                        x[i] = x_temp[i];
                    }
                    for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                        ef[i] = ef_i[i];
                    }
                    for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                        c1[i] = c1_i[i];
                        c2[i] = c2_i[i];
                        c3[i] = c3_i[i];
                    }
                    // Truss
                    for (i = 0; i < NE_TR; ++i) {
                        defllen[i] = defllen_i[i];
                    }
                    // Frame
                    for (i = 0; i < NE_FR; ++i) {
                        defllen[NE_TR+i] = defllen_i[NE_TR+i];
                        for (j = 0; j < 6; ++j) {
                            xfr[i*6+j] = xfr_temp[i*6+j];
                        }
                        for (j = 0; j < 14; ++j) {
                            efFE[i*14+j] = efFE_i[i*14+j];
                        }
                        if (yldflag[i*2] == 2) {
                            yldflag[i*2] = 0;
                        }
                        if (yldflag[i*2+1] == 2) {
                            yldflag[i*2+1] = 0;
                        }
                    }
                    // Shell
                    if (ANAFLAG == 2) {
                        for (i = 0; i < NE_SH; ++i) {
                            deffarea[i] = deffarea_i[i];
                            for (j = 0; j < 3; ++j) {
                                defslen[i*3+j] = defslen_i[i*3+j];
                            }
                        }
                    } else if (ANAFLAG == 3) {
                        for (i = 0; i < NE_SH; ++i) {
                            deffarea[i] = deffarea_i[i];
                            for (j = 0; j < 3; ++j) {
                                defslen[i*3+j] = defslen_i[i*3+j];
                                chi[i*3+j] = chi_temp[i*3+j];
                            }
                            for (j = 0; j < 9; ++j) {
                                efN[i*9+j] = efN_temp[i*9+j];
                                efM[i*9+j] = efM_temp[i*9+j];
                            }
                        }
                    }

                    /* Store displacement at DOF "k" and current load proportionality factor
                     for comparison with maximum allowable values */
                    dkc = fabs(d[dkdof]);
                    lpfc = fabs(lpf);

                    // Re-compute factor on load control
                    psi = 1;
                    for (i = 0; i < NEQ; ++i) {
                        temp = fabs(ssd[i] / ssd_o[i]);
                        if (temp < psi) {
                            psi = temp;
                        }
                    }

                    subcnt = imagcnt = negcnt = 0; // Re-initialize all counters to zero

                    // Pass control to output function
                    output (&lpf, &itecnt, d, ef, 1);
                }
            }

            /* If maximum load proportionality factor or displacement at DOF "k" is exceeded,
             i.e. if solution is successful, report statistics from solution algorithm */
            if (lpfc >= lpfmax || dkc >= dkimax) {
                fprintf(OFP[0], "\nSolution successful\n");
            }
        }

        else if (ALGFLAG == 4){ // Dynamic analysis: linear Newmark Intergration Method

            // Pass control to output function
            output (&lpfmax, &itecnt, d, ef, 0);

            // Newmark integration constants
            double alpham, alphaf, numopt, spectrds;

            // Initialize generalized total nodal displacement and internal force vectors
            for (i = 0; i < NEQ; ++i) {
                d[i] = 0;
                f[i] = 0;
            }

            // Initialize element force vectors
            for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                ef[i] = 0;
            }
            // Initialize truss deformed length variables
            for (i = 0; i < NE_TR; ++i) {
                defllen[i] = llength[i];
            }
            // Initialize frame element variables
            for (i = 0; i < NE_FR; ++i) {
                yldflag[i*2] = yldflag[i*2+1] = 0;
                defllen[NE_TR+i] = llength[NE_TR+i];
                for (j = 0; j < 14; ++j) {
                    efFE[i*14+j] = 0;
                }
            }
            // Initialize shell element variables
            for (i = 0; i < NE_SH; ++i) {
                deffarea[i] = farea[i];
                for (j = 0; j < 3; ++j) {
                    defslen[i*3+j] = slength[i*3+j];
                    chi[i*3+j] = 0;
                }
                for (j = 0; j < 9; ++j) {
                    efN[i*9+j] = 0;
                    efM[i*9+j] = 0;
                }
            }

            // Time integration parameters
            fscanf(IFP[0], "%lf,%lf\n", &numopt, &spectrds);
            
            if (numopt == 0 && spectrds != 1) {
                fprintf(OFP[0], "\n***ERROR*** Invalid spectral radius value for Newmark");
                fprintf(OFP[0], "  analysis without numerical dissipation.\n");
                goto EXIT1;
            }

            // Read in solver parameters from input file
            fscanf(IFP[0], "%lf\n", &lpfmax);
            if (OPTFLAG == 2) {
                fprintf(IFP[1], "%le\n", lpfmax);
            }

            /* Determine Newmark integration variables given numerical dissipation options
             and spectral radius */
            if(numopt == 0){
                alpham = 0;
                alphaf = 0;
            } else if (numopt == 1){
                alpham = (2*spectrds-1)/(spectrds+1);
                alphaf = spectrds/(spectrds+1);
            } else if (numopt == 2){
                alpham = 0;
                alphaf = (1-spectrds)/(1+spectrds);
            } else if (numopt == 3){
                alpham = (spectrds-1)/(spectrds+1);
                alphaf = 0;
            }
            
            /* Compute generalized total external load vector, accounting for
             generalized fixed-end load vector */
            for (i = 0; i < NEQ; ++i) {
                qtot[i] = q[i] * lpfmax;
            }

            // Initialize tangent stiffness matrix to zero
            for (i = 0; i < lss; ++i) {
                ss[i] = 0;
                sm[i] = 0;
            }

            if (NE_TR > 0) {
                // Pass control to stiff_tr function
                stiff_tr (ss, emod, carea, llength, defllen, yield, c1, c2, c3, ef, maxa,
                          mcode);
                mass_tr (sm, carea, llength, dens, x, minc, mcode, jac);
            }

            if (NE_FR > 0) {
                // Pass control to stiff_fr function
                stiff_fr (ss, emod, gmod, carea, offset, osflag, llength, defllen,
                          istrong, iweak, ipolar, iwarp, yldflag, yield, zstrong, zweak, c1,
                          c2, c3, ef, efFE, mendrel, maxa, mcode);
                mass_fr (sm, carea, llength, istrong, iweak, ipolar, iwarp, dens, osflag,
                         offset, x, xfr, minc, mcode, jac);
            }


            if (NE_SH > 0) {
                // Pass control to stiff_sh function
                stiff_sh (ss, emod, nu, x, xlocal, thick, farea, deffarea, slength,
                          defslen, yield, c1, c2, c3, ef, d, chi, efN, efM, maxa, minc, mcode);
                mass_sh (sm, carea, dens, thick, farea, slength, x, minc, mcode, jac);
            }

            if (NE_BR > 0) {
                // Pass control to stiff and mass functions
                stiff_br (ss, x, emod, nu, minc, mcode, jcode, Jinv, jac);
                mass_br (sm, dens, x, minc, mcode, jac);
            }

            /* Evaluate expression for actual dt. If actual dt < input dt, then linearlly
             interpolate between the input loads to get load, pressure and fluid acceleration
             values at each dt */
            double dtmax;
            dtmax = 1*dt;

            if (dtmax < dt) {
                dt = dtmax;
            }

            NTSTPS = ttot/dt + 1;
            double ssd;

            // Pass control to solve function
            errchk = solve (jcode, ss, ss, sm, sm, sd_fsi, r, dd, maxa, &ssd, &det, um, vm, am, uc, vc, ac, pinpt, tinpt,
                            Keff, Reff, Meff, alpham, alphaf, ipiv, 0, 1);

            // Terminate program if errors encountered
            if (errchk == 1) {
                goto EXIT2;
            }

            if (NE_TR > 0) {
                // Pass control to forces_tr function
                forces_tr (f, ef, d, emod, carea, llength, defllen, yield, c1, c2, c3,
                           mcode);
            }

            if (NE_FR > 0) {
                // Pass control to forces_fr function
                forces_fr (f, ef, ef, efFE_ref, efFE, efFE, yldflag, d, emod, gmod,
                           carea, offset, osflag, llength, defllen, istrong, iweak, ipolar,
                           iwarp, yield, zstrong, zweak, c1, c2, c3, c1, c2, c3, mendrel, mcode,
                           &dlpf, &itecnt);
            }

            if (NE_SH > 0) {
                // Pass control to forces_sh function
                forces_sh (f, ef, ef, efN, efM, d, d, chi, x, x, emod, nu, xlocal, thick,
                           farea, deffarea, slength, defslen, yield, c1, c2, c3, c1, c2, c3,
                           minc, mcode, jcode);
            }

            // Pass control to output function
            output (&lpfmax, &itecnt, d, ef, 1);

            fprintf(OFP[0], "\nSolution successful!!\n");


        }else if (ALGFLAG == 5){ // Dynamic analysis: Nonlinear Newmark Integration Method

            //Define secondary non-array variables, specific to NR algorithm
            double lpfi, lpfp, dlpfi; // Variables for NR iteration
            double ssd, sum; // Dummy variables for solve function
            double time, ddt, dt_temp, sub_dt, tsflag; // Variables for time stepping scheme
            double  a0, a1, a2, a3, a4, a5, a6, a7; // Variables for Newmark constants
            double numopt, spectrds, alpham, alphaf; //numerical dissipation options and spectral radius

            int inccnt; //Load increment counter
            int solcnt, solmin; // Minimum number of solutions
            int i, k;

            //Read in Newmark integration constants
            fscanf(IFP[0], "\n%lf,%lf\n", &numopt, &spectrds);
            
            // Read in solver parameters from input file
            fscanf(IFP[0], "%lf,%lf,%lf,%lf,%lf\n", &lpfmax, &lpf, &dlpf, &dlpfmax,
                   &dlpfmin);
            fscanf(IFP[0], "%d,%d,%d\n", &itemax, &submax, &solmin);
            fscanf(IFP[0], "%lf,%lf,%lf\n", &toldisp, &tolforc, &tolener);

            if (OPTFLAG == 2) {
                fprintf(IFP[1], "%d,%d,%d\n", itemax, submax, solmin);
            }
            
            if (numopt == 0 && spectrds != 1) {
                fprintf(OFP[0], "\n***ERROR*** Invalid spectral radius value for Newmark");
                fprintf(OFP[0], "  analysis without numerical dissipation.\n");
                goto EXIT1;
            }
            
            //Determine Newmark integration variables given numerical dissipation options
            //and spectral radius
            if(numopt == 0){
                alpham = 0;
                alphaf = 0;
            } else if (numopt == 1){
                alpham = (2*spectrds-1)/(spectrds+1);
                alphaf = spectrds/(spectrds+1);
            } else if (numopt == 2){
                alpham = 0;
                alphaf = (1-spectrds)/(1+spectrds);
            } else if (numopt == 3){
                alpham = (spectrds-1)/(spectrds+1);
                alphaf = 0;
            }
            
            alpha = pow(1-alpham+alphaf, 2)/4;
            delta = 0.5-alpham+alphaf;

            /* Evaluate expression for actual dt. If actual dt < input dt, then linearlly
             interpolate between the input loads to get load, pressure and fluid acceleration
             values at each dt */
            double dtmax;
            dtmax = 1*dt;
            if (dtmax < dt) {
                dt = dtmax;
            }

            NTSTPS = ttot/dt + 1;

            // Initialize time stepping variables
            ddt = sub_dt = 1;
            dt_temp = dt;
            tsflag = 0;
            
            lpfp = lpfi = lpf;
            dlpfi = dlpf;

            // Pass control to output function
            output (&time, &itecnt, d, ef, 0);
            k = 0;

            // Initialize generalized total nodal displacement and internal force vectors
            for (i = 0; i < NEQ; ++i) {
                d[i] = 0;
                f[i] = 0;
            }

            // Initialize element force vectors
            for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                ef[i] = 0;
            }

            // Initialize truss deformed length variables
            for (i = 0; i < NE_TR; ++i) {
                defllen[i] = llength[i];
            }

            // Initialize frame element variables
            for (i = 0; i < NE_FR; ++i) {
                yldflag[i*2] = yldflag[i*2+1] = 0;
                defllen[NE_TR+i] = llength[NE_TR+i];
                llength_temp[NE_TR+i] = llength[NE_TR+i];
                for (j = 0; j < 14; ++j) {
                    efFE[i*14+j] = 0;
                }
            }
            // Initialize shell element variables
            for (i = 0; i < NE_SH; ++i) {
                deffarea[i] = farea[i];
                for (j = 0; j < 3; ++j) {
                    defslen[i*3+j] = slength[i*3+j];
                    chi[i*3+j] = 0;
                }
                for (j = 0; j < 9; ++j) {
                    efN[i*9+j] = 0;
                    efM[i*9+j] = 0;
                }
            }

            // Initialize new displacement, velocity, and accelearation arrays
            for (i = 0; i < NEQ; ++i){
                uc[i] = vc[i] = ac [i] = 0;
            }

            //Loop through each time step to update stiffness and mass matrix
            do {
                
                //Initialize temporary time increment interval variables
                dt_temp = ddt * dt;
                sub_dt = ddt;
                tsflag = 0;

                do {
                    a0 = 1/(alpha*pow(dt_temp,2));
                    a1 = delta/(alpha*dt_temp);
                    a2 = 1/(alpha*dt_temp);
                    a3 = 1/(2*alpha) - 1;
                    a4 = delta/alpha - 1;
                    a5 = (dt_temp)/2*(delta/alpha - 2);
                    a6 = (dt_temp)*(1-delta);
                    a7 = delta*(dt_temp);

                    // Initialize load step, converged solution, and subdivision counters
                    inccnt = solcnt = subcnt = 0;
                    /* Begin load incrementation; load will be incremented until load
                     proportionality factor is equal to user specified maximum */
                    lpf = lpfp = lpfi; //reset lpf and dlpf
                    dlpf = dlpfi;

                    do {
                        // If load proportionality factor exceeds maximum, set equal to maximum
                        if (lpf > lpfmax) {
                            lpf = lpfmax;
                        }
                        /* Set all temporary variables, and variables which refer to the structure
                         in its current configuration, to values obtained at last successful
                         load increment; this step is required so as not to overwrite structure
                         properties prematurely if load increment is unsuccessful / invalid */
                        for (i = 0; i < NEQ; ++i) {
                            /* Store generalized internal force vector from previous
                             configuration */
                            fp[i] = f[i];
                            d_temp[i] = d[i];
                            f_temp[i] = f[i];
                        }
                        dlpfp = dlpf;

                        // General
                        for (i = 0; i < NJ*3; ++i) {
                            x_temp[i] = x[i];
                        }
                        for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                            ef_i[i] = ef_ip[i] = ef[i];
                        }
                        for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                            c1_i[i] = c1_ip[i] = c1[i];
                            c2_i[i] = c2_ip[i] = c2[i];
                            c3_i[i] = c3_ip[i] = c3[i];
                        }

                        for (i = 0; i < NEQ; ++i){
                            um[i] = uc_i[i] = uc[i];
                            vm[i] = vc_i[i] = vc[i];
                            am[i] = ac_i[i] = ac[i];
                        }

                        // Truss
                        for (i = 0; i < NE_TR; ++i) {
                            defllen_i[i] = defllen_ip[i] = defllen[i];
                        }
                        // Frame
                        for (i = 0; i < NE_FR; ++i) {
                            defllen_i[NE_TR+i] = defllen_ip[NE_TR+i] = defllen[NE_TR+i];
                            llength[NE_TR+i] = llength_temp[NE_TR+i];
                            for (j = 0; j < 6; ++j) {
                                xfr_temp[i*6+j] = xfr[i*6+j];
                            }
                            for (j = 0; j < 14; ++j) {
                                efFE_i[i*14+j] = efFE_ip[i*14+j] = efFE[i*14+j];
                            }
                        }
                        // Shell
                        for (i = 0; i < NE_SH; ++i) {
                            deffarea_i[i] = deffarea_ip[i] = deffarea[i];
                            for (j = 0; j < 3; ++j) {
                                defslen_i[i*3+j] = defslen_ip[i*3+j] = defslen[i*3+j];
                                chi_temp[i*3+j] = chi[i*3+j];
                            }
                            for (j = 0; j < 9; ++j) {
                                efN_temp[i*9+j] = efN[i*9+j];
                                efM_temp[i*9+j] = efM[i*9+j];
                            }
                        }

                        // Re-initialize iteration counter at the start of each increment
                        itecnt = 0;

                        /* Start of each equilibrium iteration within load increment; iterations
                         will continue until convergence is reached or iteration count exceeds
                         user specified maximum */
                        frcchk_fr = frcchk_sh = 0;

                        do {
                            /* The load history is under linear interpolation assumption in the case when damping scheme is applied. Limit the size of delta T that may be used to maintain fidelity with the loading history. */
                            if (itecnt == 0) {// Predictor step
                                if (k == 0){
                                    for (i = 0; i < NEQ; ++i){
                                        /* Compute generalized total external load vector, accounting for
                                         generalized fixed-end load vector */
                                        qtot[i] = pinpt[i*NTSTPS+k]*sub_dt*lpf;
                                        // Compute residual force vector
                                        r[i] = qtot[i] - f_temp[i];
                                    }
                                }else {
                                    for (i = 0; i < NEQ; ++i) {
                                        /* Compute generalized total external load vector, accounting for
                                         generalized fixed-end load vector */
                                        qtot[i] = (pinpt[i*NTSTPS+k-1]+(pinpt[i*NTSTPS+k]-pinpt[i*NTSTPS+k-1])*(sub_dt))*lpf;
                                        // Compute residual force vector
                                        r[i] = (qtot[i]-f_temp[i]) + alphaf/(1-alphaf)*(pinpt[i*NTSTPS+k-1]-f_temp[i]);
                                    }
                                }
                            } else { // Corrector steps
                                // Compute residual force vector
                                if (SLVFLAG == 0) {
                                    for (i = 0; i< NEQ; ++i){
                                        r[i] = (f_temp[i]-qtot[i]);
                                        r[i] = r[i]+sm[i]*ac_i[i]-sm[i]*((a2*vc_i[i]+a3*ac_i[i])*(1-alpham)-ac_i[i]*alpham)/(1-alphaf);
                                    }
                                } else if (SLVFLAG == 1){
                                    for (i = 0; i < NEQ; ++i) {
                                        sum = 0;
                                        for (j = 0; j < NEQ; ++j) {
                                            sum += sm[i*NEQ+j];
                                        }
                                        r[i] = (f_temp[i]-qtot[i]);
                                        r[i] = r[i]+sum*ac_i[i]-sum*((a2*vc_i[i]+a3*ac_i[i])*(1-alpham)-ac_i[i]*alpham)/(1-alphaf);
                                    }
                                }
                            }
                            
                            for (i = 0; i < lss; ++i) {
                                ss[i] = 0;
                                sm[i] = 0;
                            }

                            if (NE_TR > 0) {
                                // Pass control to stiff_tr and mass_tr function
                                stiff_tr (ss, emod, carea, llength_temp, defllen_ip, yield, c1_ip,
                                          c2_ip, c3_ip, ef_ip, maxa, mcode);
                                mass_tr (sm, carea, llength_temp, dens, x, minc, mcode, jac);
                            }

                            if (NE_FR > 0) {
                                // Pass control to stiff_fr and mass_fr function
                                stiff_fr (ss, emod, gmod, carea, offset, osflag, llength_temp,
                                          defllen_ip, istrong, iweak, ipolar, iwarp, yldflag,
                                          yield, zstrong, zweak, c1_ip, c2_ip, c3_ip, ef_ip,
                                          efFE_ip, mendrel, maxa, mcode);
                                mass_fr (sm, carea, llength_temp, istrong, iweak, ipolar, iwarp, dens, osflag,
                                         offset, x, xfr, minc, mcode, jac);
                            }

                            if (NE_SH > 0) {
                                // Pass control to stiff_sh and mass_sh function
                                stiff_sh (ss, emod, nu, x_temp, xlocal, thick, farea,
                                          deffarea_ip, slength, defslen_ip, yield, c1_ip, c2_ip,
                                          c3_ip, ef_ip, d_temp, chi_temp, efN_temp, efM_temp, maxa,
                                          minc, mcode);
                                mass_sh (sm, carea, dens, thick, farea, slength, x, minc, mcode, jac);
                            }

                            if (lss == 1) {
                                /* Carry out computation of incremental displacement directly for
                                 lss = 1 */
                                dd[0] = r[0] / ss[0];
                            } else {
                                // Pass control to solve function
                                errchk = solve (jcode, ss, ss, sm, sm, sd_fsi, r, dd, maxa, &ssd, &det, uc_i, vc_i, ac_i, um, vm, am, qtot, tinpt,
                                                Keff, Reff, Meff, alpham, alphaf, ipiv, 0, ddt);

                                // Terminate program if errors encountered
                                if (errchk == 1) {
                                    goto EXIT2;
                                }
                            }

                            /* Update generalized total nodal displacement vector, store
                             generalized internal force vector from previous iteration, and
                             re-initialize generalized internal force vector */
                            for (i = 0; i < NEQ; ++i) {
                                if (itecnt > 0 ){
                                    dd[i] = dd[i] * (-1);
                                }
                                d_temp[i] += dd[i];
                                f_ip[i] = f_temp[i];
                                f_temp[i] = 0;
                            }

                            // Calculate displacements, velocities and accelerations
                            for (i = 0; i < NEQ; ++i) {
                                uc_i[i] = d_temp[i];
                                ac_i[i] = (uc_i[i] - um[i])*a0 - a2*vm[i] - a3*am[i];
                                vc_i[i] = vm[i] + a6*am[i] + a7*ac_i[i];
                            }


                            // Pass control to updatc function
                            updatc (x_temp, x_ip, xfr_temp, dd, defllen_i, deffarea_i, defslen_i,
                                    offset, osflag, auxpt, c1_i, c2_i, c3_i, minc, jcode);

                            if (NE_TR > 0) {
                                // Pass control to forces_tr function
                                forces_tr (f_temp, ef_i, d, emod, carea, llength_temp, defllen_i,
                                           yield, c1_i, c2_i, c3_i, mcode);
                            }

                            if (NE_FR > 0) {
                                // Pass control to forces_fr function
                                frcchk_fr = forces_fr (f_temp, ef_ip, ef_i, efFE_ref, efFE_ip,
                                                       efFE_i, yldflag, dd, emod, gmod, carea, offset, osflag,
                                                       llength_temp, defllen_ip, istrong, iweak, ipolar, iwarp, yield,
                                                       zstrong, zweak, c1_ip, c2_ip, c3_ip, c1_i, c2_i, c3_i,
                                                       mendrel, mcode, &dlpf, &itecnt);
                            }

                            if (NE_SH > 0) {
                                // Pass control to forces_sh function
                                frcchk_sh = forces_sh (f_temp, ef_ip, ef_i, efN_temp, efM_temp,
                                                       dd, d_temp, chi_temp, x_temp, x_ip, emod, nu, xlocal, thick,
                                                       farea, deffarea_ip, slength, defslen_ip, yield, c1_ip, c2_ip,
                                                       c3_ip, c1_i, c2_i, c3_i, minc, mcode, jcode);
                            }

                            //Compute out-of-balance dynamic forces
                            if (SLVFLAG == 0) {
                                for (i = 0; i < NEQ; ++i){
                                    dyn[i] = qtot[i] - sm[i]*ac_i[i];
                                }
                            }
                            else if (SLVFLAG == 1) { // using CLAPACK solver
                                for (i = 0; i < NEQ; ++i) {
                                    sum = 0;
                                    for (j = 0; j < NEQ; ++j) {
                                        sum += sm[i*NEQ+j];
                                    }
                                    dyn[i] = qtot[i] - sum*ac_i[i];
                                }
                            }

                            if (itecnt == 0) {
                                // Compute internal energy from first iteration
                                intener1 = 0;
                                for (i = 0; i < NEQ; ++i) {
                                    intener1 += dd[i] * (dyn[i] - fp[i]);
                                }
                            }

                            errchk = test (d_temp, dd, f_temp, fp, dyn, f_ip, &intener1,
                                           &convchk, &toldisp, &tolforc, &tolener);

                            // Terminate program if errors encountered
                            if (errchk == 1) {
                                goto EXIT2;
                            }

                            // Update element internal forces from previous iteration
                            for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                                ef_ip[i] = ef_i[i];
                            }

                            // Update variables from previous iteration
                            // General
                            for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                                c1_ip[i] = c1_i[i];
                                c2_ip[i] = c2_i[i];
                                c3_ip[i] = c3_i[i];
                            }

                            // Truss
                            for (i = 0; i < NE_TR; ++i) {
                                defllen_ip[i] = defllen_i[i];
                            }
                            // Frame
                            for (i = 0; i < NE_FR; ++i) {
                                defllen_ip[NE_TR+i] = defllen_i[NE_TR+i];
                                for (j = 0; j < 14; ++j) {
                                    efFE_ip[i*14+j] = efFE_i[i*14+j];
                                }
                            }
                            // Shell
                            for (i = 0; i < NE_SH; ++i) {
                                deffarea_ip[i] = deffarea_i[i];
                                for (j = 0; j < 3; ++j) {
                                    defslen_ip[i*3+j] = defslen_i[i*3+j];
                                }
                            }

                            itecnt++; // Advance iteration counter
                            
                        } while (convchk != 0 && frcchk_fr == 0 && frcchk_sh == 0 && itecnt <= itemax);

                        if (frcchk_fr == 2) {
                            dlpf = dlpfp; // Reset increment in load proportionality factor
                        } else if ((convchk != 0 || frcchk_fr != 0 || frcchk_sh != 0) &&
                                   subcnt <= submax) {
                            if (lpf == lpfmax) {
                                fprintf(OFP[0], "\n***ERROR*** Maximum allowable load");
                                fprintf(OFP[0], " proportionality factor attempted without");
                                fprintf(OFP[0], " convergence at current time step\n");

                                break;
                            } else if (dlpfp == dlpfmin) {
                                fprintf(OFP[0], "\n***ERROR*** Minimum allowable increment of");
                                fprintf(OFP[0], " load proportionality factor reached at current time step\n");

                                break;
                            }
                            if (frcchk_fr != 1) {
                                // Decrease increment of load proportionality factor
                                dlpf = dlpfp / 2;
                            }
                            if (dlpf < dlpfmin) {
                                dlpf = dlpfmin;
                            }

                            // Step back load proportionality factor
                            lpf = lpf - dlpfp + dlpf;

                            subcnt++; // Advance subdivision counter
                            solcnt = 0; // Re-initialize converged solution counter
                        } else if (subcnt > submax) {
                            fprintf(OFP[0], "\n***ERROR*** Maximum allowable number of");
                            fprintf(OFP[0], " subdivisions exceeded at current time step\n");

                            break;
                        } else {
                            inccnt++; // Advance load increment counter

                            solcnt++;
                            subcnt = 0; // Re-initialize subdivision counter
                            /* If current load increment resulted in solmin converged solutions
                             in a row, increase increment in load proportionality factor */
                            if (solcnt >= solmin) {
                                dlpf *= 2; // Increase increment of load proportionality factor
                                solcnt = 0; // Re-initialize solution counter
                            }

                            lpf += dlpf; // Increment load proportionality factor
                        }
                    } while (lpf <= lpfmax);

                    if (convchk != 0) {
                        if (frcchk_fr != 0 || frcchk_sh != 0 || convchk != 0){
                            fprintf(OFP[0], "\nTime increment interval reduced\n");
                            
                            // Reduce time increment interval by half if solution does not converge
                            ddt = ddt/2;
                            // Compute temporary time increment interval
                            dt_temp = ddt * dt;

                            if (tsflag == 0){
                                // Initialize sub-time increment interval to time increment multiplier if solution exceeded yield surface
                                sub_dt = ddt;
                            } else if (tsflag == 1){
                                // Set sub-time increment interval to its previous value if solution exceeded yield surface after some time increment iterations
                                sub_dt = sub_dt - ddt;
                            }

                        } else {
                            fprintf(OFP[0], "\n***ERROR*** Solutions failed to converge\n");
                            goto EXIT2;
                        }

                    } else {
                        // Update all permanent variables to values which represent structure
                        //in its current configuration
                        for (i = 0; i < NEQ; ++i) {
                            d[i] = d_temp[i];
                            f[i] = f_temp[i];
                        }

                        for (i = 0; i < NEQ; ++i){
                            uc[i] = uc_i[i];
                            vc[i] = vc_i[i];
                            ac[i] = ac_i[i];
                        }
                        
                        // General
                        for (i = 0; i < NE_TR*2+NE_FR*14+NE_SH*18; ++i) {
                            ef[i] = ef_i[i];
                        }
                        for (i = 0; i < NJ*3; ++i) {
                            x[i] = x_temp[i];
                        }
                        for (i = 0; i < NE_TR+NE_FR*3+NE_SH*3; ++i) {
                            c1[i] = c1_i[i];
                            c2[i] = c2_i[i];
                            c3[i] = c3_i[i];
                        }

                        // Truss
                        for (i = 0; i < NE_TR; ++i) {
                            defllen[i] = defllen_i[i];
                        }

                        // Frame
                        for (i = 0; i < NE_FR; ++i) {
                            defllen[NE_TR+i] = defllen_i[NE_TR+i];
                            llength[NE_TR+i] = llength_temp[NE_TR+i];
                            for (j = 0; j < 14; ++j) {
                                efFE[i*14+j] = efFE_i[i*14+j];
                            }
                            for (j = 0; j < 6; ++j) {
                                xfr[i*6+j] = xfr_temp[i*6+j];
                            }
                            if (yldflag[i*2] == 2) {
                                yldflag[i*2] = 0;
                            }
                            if (yldflag[i*2+1] == 2) {
                                yldflag[i*2+1] = 0;
                            }
                        }

                        // Shell
                        if (ANAFLAG == 2) {
                            for (i = 0; i < NE_SH; ++i) {
                                deffarea[i] = deffarea_i[i];
                                for (j = 0; j < 3; ++j) {
                                    defslen[i*3+j] = defslen_i[i*3+j];
                                }
                            }
                        } else {
                            for (i = 0; i < NE_SH; ++i) {
                                deffarea[i] = deffarea_i[i];
                                for (j = 0; j < 3; ++j) {
                                    defslen[i*3+j] = defslen_i[i*3+j];
                                    chi[i*3+j] = chi_temp[i*3+j];
                                }
                                for (j = 0; j < 9; ++j) {
                                    efN[i*9+j] = efN_temp[i*9+j];
                                    efM[i*9+j] = efM_temp[i*9+j];
                                }
                            }
                        }
                    }
                    
                    // Update sub-time increment interval and time stepping flag if solution passes force check
                    if (frcchk_fr == 0 && frcchk_sh == 0){
                        if (ddt < 1) {
                            if (sub_dt <= 1) {
                                sub_dt = sub_dt + ddt;
                                tsflag = 1;
                            }
                            if (sub_dt > 1){
                                ddt = ddt * 2;
                                tsflag = 2;
                            }
                        } else if (ddt == 1){
                            tsflag = 2;
                        }
                    }
                } while (ddt >= 0.0001 && tsflag != 2);

                if (convchk != 0 || frcchk_fr != 0 || frcchk_sh != 0) {

                    fprintf(OFP[0], "\n***ERROR*** Minimum time increment interval reached");
                    fprintf(OFP[0], " without convergence\n");

                    goto EXIT2;
                }

                time = k * dt;
                
                //Pass control to output function
                output (&time, &itecnt, d, ef, 1);

                ++k;

            } while (k < NTSTPS);

            if (convchk == 0) {
                fprintf(OFP[0], "\nSolution successful!!\n");
            }

        }
    }
    // Pass control to free_all function
    return free_all (p2p2i, ni, p2p2l, nl, p2p2d, nd, 0);

EXIT1:
    // Pass control to closeio function
    return closeio(1);
EXIT2:
    // Pass control to free_all function
    return free_all (p2p2i, ni, p2p2l, nl, p2p2d, nd, 1);
}


