//********************************************************************************
//**																			**
//**  Pertains to CU-BEN ver 4.0												**
//**																			**
//**  Copyright (c) 2018 C. J. Earls                                            **
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

#define phitol 1e-4 // Allowable +/- deviation from 1.0 of phi

extern long NJ, NE_TR, NE_FR, NEQ;
extern int ANAFLAG, ALGFLAG, OPTFLAG, SLVFLAG;
extern FILE *IFP[4], *OFP[8];

void prop_fr (double *px, double *pxfr, double *pemod, double *pgmod, double *pdens, double *poffset,
    int *posflag, double *pauxpt, double *pcarea, double *pllength, double *pistrong,
    double *piweak, double *pipolar, double *piwarp, double *pyield, double *pzstrong,
    double *pzweak, double *pc1, double *pc2, double *pc3, int *pmendrel, long *pminc)
{
    // Initialize function variables
    int xmendrel[4];
    long i, j, k, l, ptr;
    double el[3], temp[3], localx[3], localy[3], localz[3];
    double xoffset[6];

    fprintf(OFP[0], "\nFrame Element Properties:\n\tElement\t\tElastic Modulus\t\t");
    fprintf(OFP[0], "Shear Modulus\t\tArea\t\t___________________________Moments of");
    fprintf(OFP[0], " Inertia__________________________\n\t\t\t\t\t\t\t\t\t\t\t");
    fprintf(OFP[0], "Strong-Axis Bending\tWeak-Axis Bending\tPolar\t\tWarping\n");

    ptr = NE_TR * 2;
    for (i = 0; i < NE_FR; ++i) {
        // Read in element properties from input file
        fscanf(IFP[0], "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", pemod+NE_TR+i, pgmod+i, pdens+i,
            pcarea+NE_TR+i, pistrong+i, piweak+i, pipolar+i, piwarp+i);
        // Write out element properties to output file
        fprintf(OFP[0], "\t%ld\t\t%lf\t\t%lf\t\t%lf\t%lf\t\t%lf\t\t%lf\t%lf\n", i + 1,
            *(pemod+NE_TR+i), *(pgmod+i), *(pcarea+NE_TR+i), *(pistrong+i), *(piweak+i),
            *(pipolar+i), *(piwarp+i));

        /* Initialize member end offset flags; zero indicates no member end offsets
           present */
        *(posflag+i) = 0;
    }
    if (OPTFLAG == 2) {
		for (i = 0; i < NE_FR; ++i) {
			fprintf(IFP[1], "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
				*(pemod+NE_TR+i), *(pgmod+i), *(pcarea+NE_TR+i), *(pistrong+i), *(piweak+i),
				*(pipolar+i), *(piwarp+i));
		}
    }

    // Read in member end offsets from input file
    fscanf(IFP[0], "%ld,%lf,%lf,%lf,%lf,%lf,%lf\n", &j, &xoffset[0], &xoffset[1],
        &xoffset[2], &xoffset[3], &xoffset[4], &xoffset[5]);
    if (j != 0) {
        fprintf(OFP[0], "\nFrame Member End Offsets:\n\t\t___________________End-1");
        fprintf(OFP[0], "___________________\t___________________End-2");
        fprintf(OFP[0], "___________________\n\tMember\tDirection-1\tDirection-2\t");
        fprintf(OFP[0], "Direction-3\tDirection-1\tDirection-2\tDirection-3\n");
        do {
        	if (OPTFLAG == 2) {
				fprintf(IFP[1], "%ld,%e,%e,%e,%e,%e,%e\n",
					j, xoffset[0], xoffset[1], xoffset[2], xoffset[3], xoffset[4],
					xoffset[5]);
        	}

            j--; // Decrement for convenience
            *(posflag+j) = 1;
            *(poffset+j*6) = xoffset[0];
            *(poffset+j*6+1) = xoffset[1];
            *(poffset+j*6+2) = xoffset[2];
            *(poffset+j*6+3) = xoffset[3];
            *(poffset+j*6+4) = xoffset[4];
            *(poffset+j*6+5) = xoffset[5];

            // Write out member end offsets to output file
            fprintf(OFP[0], "\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\n", j + 1,
                *(poffset+j*6), *(poffset+j*6+1), *(poffset+j*6+2), *(poffset+j*6+3),
                *(poffset+j*6+4), *(poffset+j*6+5));
            // Read in member end offsets from input file
            fscanf(IFP[0], "%ld,%lf,%lf,%lf,%lf,%lf,%lf\n", &j, &xoffset[0], &xoffset[1],
                &xoffset[2], &xoffset[3], &xoffset[4], &xoffset[5]);
        } while (j != 0);
    }
    if (OPTFLAG == 2) {
		fprintf(IFP[1], "0,0,0,0,0,0,0\n");
    }

    fprintf(OFP[0], "\nFrame Element Orientation:\n\tElement\t\tLength\t\t");
    fprintf(OFP[0], "______________Auxiliary Points_____________\n\t\t\t\t\t");
    fprintf(OFP[0], "Direction-1\tDirection-2\tDirection-3\n");

    for (i = 0; i < NE_FR; ++i) {
        // Compute element length
        j = *(pminc+ptr+i*2) - 1;
        k = *(pminc+ptr+i*2+1) - 1;
        if (*(posflag+i) == 0) {
            for (l = 0; l < 3; ++l) {
                *(pxfr+i*6+l) = *(px+j*3+l);
                *(pxfr+i*6+3+l) = *(px+k*3+l);
            }
        } else {
            for (l = 0; l < 3; ++l) {
                *(pxfr+i*6+l) = *(px+j*3+l) + *(poffset+i*6+l);
                *(pxfr+i*6+3+l) = *(px+k*3+l) + *(poffset+i*6+3+l);
            }
        }
        el[0] = *(pxfr+i*6+3) - *(pxfr+i*6);
        el[1] = *(pxfr+i*6+3+1) - *(pxfr+i*6+1);
        el[2] = *(pxfr+i*6+3+2) - *(pxfr+i*6+2);
        *(pllength+NE_TR+i) = sqrt(dot(el,el,3));

        // Read in element auxiliary point from input file
        fscanf(IFP[0], "%lf,%lf,%lf\n", pauxpt+i*3, pauxpt+i*3+1, pauxpt+i*3+2);

        // Compute direction cosines
        localx[0] = el[0] / *(pllength+NE_TR+i);
        localx[1] = el[1] / *(pllength+NE_TR+i);
        localx[2] = el[2] / *(pllength+NE_TR+i);
        temp[0] = *(pauxpt+i*3) - *(pxfr+i*6);
        temp[1] = *(pauxpt+i*3+1) - *(pxfr+i*6+1);
        temp[2] = *(pauxpt+i*3+2) - *(pxfr+i*6+2);
        cross(localx,temp,localz,1);
        cross(localz,localx,localy,1);
        for (j = 0; j < 3; ++j) {
            *(pc1+NE_TR+i*3+j) = localx[j];
            *(pc2+NE_TR+i*3+j) = localy[j];
            *(pc3+NE_TR+i*3+j) = localz[j];
        }

        // Write out element length and auxiliary point from input file
        fprintf(OFP[0], "\t%ld\t\t%lf\t%lf\t%lf\t%lf\n", i + 1, *(pllength+NE_TR+i),
            *(pauxpt+i*3), *(pauxpt+i*3+1), *(pauxpt+i*3+2));
    }
    if (OPTFLAG == 2) {
		for (i = 0; i < NE_FR; ++i) {
			fprintf(IFP[1], "%lf,%lf,%lf\n", *(pauxpt+i*3), *(pauxpt+i*3+1),
				*(pauxpt+i*3+2));
		}
    }

    fprintf(OFP[0], "\nFrame Element Yield Criteria:\n\tElement\t\tYield Stress\t\t");
    fprintf(OFP[0], "_____Plastic Section Moduli______\n\t\t\t\t\t\tStrong-Axis\t\t");
    fprintf(OFP[0], "Weak-Axis\n");

    for (i = 0; i < NE_FR; ++i) {
        // Read in element yield criteria from input file
        fscanf(IFP[0], "%lf,%lf,%lf\n", pyield+NE_TR+i, pzstrong+i, pzweak+i);
        // Write out element yield criteria to output file
        fprintf(OFP[0], "\t%ld\t\t%lf\t\t%lf\t\t%lf\n", i + 1, *(pyield+NE_TR+i),
            *(pzstrong+i), *(pzweak+i));

        /* Initialize member end bending releases; zero indicates no member end bending
           releases present */
        *(pmendrel+i*5) = 0;
    }
    if (OPTFLAG == 2) {
		for (i = 0; i < NE_FR; ++i) {
			fprintf(IFP[1], "%lf,%lf,%lf\n", *(pyield+NE_TR+i), *(pzstrong+i),
				*(pzweak+i));
		}
    }

    // Read in member end releases from input file
    fscanf(IFP[0], "%ld,%d,%d,%d,%d\n", &j, &xmendrel[0], &xmendrel[1], &xmendrel[2],
        &xmendrel[3]);
    if (j != 0) {
        fprintf(OFP[0], "\nFrame Member End Bending Releases:\n\t\t____________End-1");
        fprintf(OFP[0], "____________\t____________End-2____________\n\tMember\t");
        fprintf(OFP[0], "Strong-Axis\tWeak-Axis\tStrong-Axis\tWeak-Axis\n");
        do {
        	if (OPTFLAG == 2) {
				fprintf(IFP[1], "%ld,%d,%d,%d,%d\n", j, xmendrel[0], xmendrel[1],
					xmendrel[2], xmendrel[3]);
        	}

            j--; // Decrement for convenience
            *(pmendrel+j*5) = 1;
            *(pmendrel+j*5+1) = xmendrel[0];
            *(pmendrel+j*5+2) = xmendrel[1];
            *(pmendrel+j*5+3) = xmendrel[2];
            *(pmendrel+j*5+4) = xmendrel[3];

            // Write out member end releases to output file
            fprintf(OFP[0], "\t%ld\t%d\t\t%d\t\t%d\t\t%d\n", j + 1, *(pmendrel+j*5+1),
                *(pmendrel+j*5+2), *(pmendrel+j*5+3), *(pmendrel+j*5+4));
            // Read in member end releases from input file
            fscanf(IFP[0], "%ld,%d,%d,%d,%d\n", &j, &xmendrel[0], &xmendrel[1],
                &xmendrel[2], &xmendrel[3]);
        } while (j != 0);
    }
    if (OPTFLAG == 2) {
		fprintf(IFP[1], "0,0,0,0,0\n");
    }
}

void stiff_fr (double *pss, double *pemod, double *pgmod, double *pcarea,
    double *poffset, int *posflag, double *pllength, double *pdefllen_ip,
    double *pistrong, double *piweak, double *pipolar, double *piwarp, int *pyldflag,
    double *pyield, double *pzstrong, double *pzweak, double *pc1_ip, double *pc2_ip,
    double *pc3_ip, double *pef_ip, double *pefFE_ip, int *pmendrel, long *pmaxa,
    long *pmcode)
{
    // Initialize primary function variables
    long i, j, k, n, je, ie, ptr, ptr2;

    double k_fr[14][14]; // General element stiffness matrix in local coordinate system
    double T_ip[14][14]; // Coordinate transformation matrix
    /* General element stiffness matrix in global coordinate system, w.r.t. local joints
       i and j */
    double Kij_fr[14][14];
    double T_rl[14][14]; // Rigid link transformation matrix
    /* General element stiffness matrix in global coordinate system, transformed (?) to
       global joints 1 and 2 */
    double K12_fr[14][14];
    // Total element force vector in previous local coordinate system
    double eftot_ip[14];

    ptr = NE_TR * 2;
    ptr2 = NE_TR * 6;
    for (n = 0; n < NE_FR; ++n) {
        // Initialize all elements to zero and compute total element force vector
        for (i = 0; i < 14; ++i) {
            for (j = 0; j < 14; ++j) {
                k_fr[i][j] = 0;
                T_ip[i][j] = T_rl[i][j] = 0;
            }
            eftot_ip[i] = *(pef_ip+ptr+n*14+i) + *(pefFE_ip+n*14+i);
        }

        // Pass control to stiffe_fr function
        stiffe_fr (&k_fr[0][0], pemod, pgmod, pcarea, pllength, pistrong, piweak,
            pipolar, piwarp, n);

        /* Include nonlinear subroutines in element stiffness matrix assembly depending
           upon user-requested analysis */
        if (ANAFLAG == 2) {
            // Pass control to stiffg_fr function
            stiffg_fr (&k_fr[0][0], eftot_ip, pdefllen_ip, pcarea, pipolar, n);
        } else if (ANAFLAG == 3) {
            // Pass control to stiffg_fr function
            stiffg_fr (&k_fr[0][0], eftot_ip, pdefllen_ip, pcarea, pipolar, n);

            if (*(pyldflag+n*2) != 2 || *(pyldflag+n*2+1) != 2) {
	            // Pass control to stiffm_fr function
	            stiffm_fr (&k_fr[0][0], eftot_ip, pyldflag, pyield, pcarea, pzstrong,
	                pzweak, n);
			}
        }

        // Check for member end bending releases
        if (*(pmendrel+n*5) == 1) {
            release (&k_fr[0][0], pmendrel, n);
        }

        // Assign non-zero elements of coordinate transformation matrix
        T_ip[0][0] = T_ip[3][3] = T_ip[7][7] = T_ip[10][10] = *(pc1_ip+NE_TR+n*3);
        T_ip[0][1] = T_ip[3][4] = T_ip[7][8] = T_ip[10][11] = *(pc1_ip+NE_TR+n*3+1);
        T_ip[0][2] = T_ip[3][5] = T_ip[7][9] = T_ip[10][12] = *(pc1_ip+NE_TR+n*3+2);
        T_ip[1][0] = T_ip[4][3] = T_ip[8][7] = T_ip[11][10] = *(pc2_ip+NE_TR+n*3);
        T_ip[1][1] = T_ip[4][4] = T_ip[8][8] = T_ip[11][11] = *(pc2_ip+NE_TR+n*3+1);
        T_ip[1][2] = T_ip[4][5] = T_ip[8][9] = T_ip[11][12] = *(pc2_ip+NE_TR+n*3+2);
        T_ip[2][0] = T_ip[5][3] = T_ip[9][7] = T_ip[12][10] = *(pc3_ip+NE_TR+n*3);
        T_ip[2][1] = T_ip[5][4] = T_ip[9][8] = T_ip[12][11] = *(pc3_ip+NE_TR+n*3+1);
        T_ip[2][2] = T_ip[5][5] = T_ip[9][9] = T_ip[12][12] = *(pc3_ip+NE_TR+n*3+2);
        T_ip[6][6] = T_ip[13][13] = 1;

        if (*(posflag+n) == 0) {
            // Pass control to transform function
            transform (&k_fr[0][0], &T_ip[0][0], &K12_fr[0][0], 14);
        } else {
            // Pass control to transform function
            transform (&k_fr[0][0], &T_ip[0][0], &Kij_fr[0][0], 14);

            /* Assign non-zero elements of TRANSPOSE of rigid link transformation matrix;
               this allows for use of transform function without modification */
            for (i = 0; i < 14; ++i) {
                T_rl[i][i] = 1;
            }
            T_rl[1][3] = -(*(poffset+n*6+2));
            T_rl[2][3] = *(poffset+n*6+1);
            T_rl[0][4] = *(poffset+n*6+2);
            T_rl[2][4] = -(*(poffset+n*6));
            T_rl[0][5] = -(*(poffset+n*6+1));
            T_rl[1][5] = *(poffset+n*6);
            T_rl[8][10] = -(*(poffset+n*6+5));
            T_rl[9][10] = *(poffset+n*6+4);
            T_rl[7][11] = *(poffset+n*6+5);
            T_rl[9][11] = -(*(poffset+n*6+3));
            T_rl[7][12] = -(*(poffset+n*6+4));
            T_rl[8][12] = *(poffset+n*6+3);

            // Pass control to transform function
            transform (&Kij_fr[0][0], &T_rl[0][0], &K12_fr[0][0], 14);
        }
		
		if (SLVFLAG == 0){
			/*Initialize index and then assign element tangent stiffness coefficients of
			element n to the structure stiffness matrix by index, mcode, and maxa */
			for (je = 0; je < 14; ++je) {
				j = *(pmcode+ptr2+n*14+je);
				if (j != 0) {
					// Check mcode above current entry to find rank of "j"
					for (ie = 0; ie <= je; ++ie) {
						i = *(pmcode+ptr2+n*14+ie);
						if (i != 0) {
							if (i > j) { // Find element address as diagonal address + delta
								k = *(pmaxa+i-1) + (i - j);
							} else {
								k = *(pmaxa+j-1) + (j - i);
							}
							/* Add current element stiffness to previous elements'
							 contributions to the given DOFs */
							*(pss+k-1) += K12_fr[ie][je];
						}
					}
				}
			}
		}
		else {
			/* Build the full order (i.e. [NEQ][NEQ] system stiffness matrix using mcode */
			for (ie = 0; ie < 14; ++ie) {
				for (je = 0; je < 14; ++je) {
					i = *(pmcode+ptr2+n*14+ie);
					j = *(pmcode+ptr2+n*14+je);
					if ((i != 0) && (j != 0)) {
						*(pss+(i-1)*NEQ+j-1) += K12_fr[je][ie];
					}
				}
			}
		}
	}
}

void stiffe_fr (double *pk_fr, double *pemod, double *pgmod, double *pcarea,
    double *pllength, double *pistrong, double *piweak, double *pipolar, double *piwarp,
    long n)
{
    *(pk_fr+0*14+0) = *(pk_fr+7*14+7) =
        *(pemod+NE_TR+n) * (*(pcarea+NE_TR+n)) / *(pllength+NE_TR+n);
    *(pk_fr+7*14+0) = *(pk_fr+0*14+7) = -(*(pk_fr+0*14+0));
    *(pk_fr+1*14+1) = *(pk_fr+8*14+8) =
        12 * (*(pemod+NE_TR+n)) * (*(pistrong+n)) / pow(*(pllength+NE_TR+n),3);
    *(pk_fr+8*14+1) = *(pk_fr+1*14+8) = -(*(pk_fr+1*14+1));
    *(pk_fr+2*14+2) = *(pk_fr+9*14+9) =
        12 * (*(pemod+NE_TR+n)) * (*(piweak+n)) / pow(*(pllength+NE_TR+n),3);
    *(pk_fr+9*14+2) = *(pk_fr+2*14+9) = -(*(pk_fr+2*14+2));
    *(pk_fr+3*14+3) = *(pk_fr+10*14+10) =
        6 * (*(pgmod+n)) * (*(pipolar+n)) / (5 * (*(pllength+NE_TR+n))) +
        12 * (*(pemod+NE_TR+n)) * (*(piwarp+n)) / pow(*(pllength+NE_TR+n),3);
    *(pk_fr+10*14+3) = *(pk_fr+3*14+10) = -(*(pk_fr+3*14+3));
    *(pk_fr+5*14+5) = *(pk_fr+12*14+12) =
        4 * (*(pemod+NE_TR+n)) * (*(pistrong+n)) / *(pllength+NE_TR+n);
    *(pk_fr+4*14+4) = *(pk_fr+11*14+11) =
        4 * (*(pemod+NE_TR+n)) * (*(piweak+n)) / *(pllength+NE_TR+n);
    *(pk_fr+6*14+6) = *(pk_fr+13*14+13) =
        2 * (*(pgmod+n)) * (*(pipolar+n)) * (*(pllength+NE_TR+n)) / 15 +
        4 * (*(pemod+NE_TR+n)) * (*(piwarp+n)) / *(pllength+NE_TR+n);
    *(pk_fr+5*14+1) = *(pk_fr+1*14+5) = *(pk_fr+12*14+1) = *(pk_fr+1*14+12) =
        6 * (*(pemod+NE_TR+n)) * (*(pistrong+n)) / pow(*(pllength+NE_TR+n),2);
    *(pk_fr+8*14+5) = *(pk_fr+5*14+8) = *(pk_fr+12*14+8) = *(pk_fr+8*14+12) =
        -(*(pk_fr+5*14+1));
    *(pk_fr+9*14+4) = *(pk_fr+4*14+9) = *(pk_fr+11*14+9) = *(pk_fr+9*14+11) =
        6 * (*(pemod+NE_TR+n)) * (*(piweak+n)) / pow(*(pllength+NE_TR+n),2);
    *(pk_fr+4*14+2) = *(pk_fr+2*14+4) = *(pk_fr+11*14+2) = *(pk_fr+2*14+11) =
        -(*(pk_fr+9*14+4));
    *(pk_fr+6*14+3) = *(pk_fr+3*14+6) = *(pk_fr+13*14+3) = *(pk_fr+3*14+13) =
        *(pgmod+n) * (*(pipolar+n)) / 10 +
        6 * (*(pemod+NE_TR+n)) * (*(piwarp+n)) / pow(*(pllength+NE_TR+n),2);
    *(pk_fr+10*14+6) = *(pk_fr+6*14+10) = *(pk_fr+13*14+10) = *(pk_fr+10*14+13) =
        -(*(pk_fr+6*14+3));
    *(pk_fr+12*14+5) = *(pk_fr+5*14+12) =
        2 * (*(pemod+NE_TR+n)) * (*(pistrong+n)) / *(pllength+NE_TR+n);
    *(pk_fr+11*14+4) = *(pk_fr+4*14+11) =
        2 * (*(pemod+NE_TR+n)) * (*(piweak+n)) / *(pllength+NE_TR+n);
    *(pk_fr+13*14+6) = *(pk_fr+6*14+13) =
        -(*(pgmod+n) * (*(pipolar+n)) * (*(pllength+NE_TR+n)) / 30 -
        2 * (*(pemod+NE_TR+n)) * (*(piwarp+n)) / *(pllength+NE_TR+n));
}

void stiffg_fr (double *pk_fr, double *peftot_ip, double *pdefllen_ip, double *pcarea,
    double *pipolar, long n)
{
    *(pk_fr+0*14+0) += *(peftot_ip+7) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+7*14+7) += *(peftot_ip+7) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+7*14+0) -= *(peftot_ip+7) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+0*14+7) -= *(peftot_ip+7) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+1*14+1) += 6 * (*(peftot_ip+7)) / (5 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+8*14+8) += 6 * (*(peftot_ip+7)) / (5 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+2*14+2) += 6 * (*(peftot_ip+7)) / (5 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+9*14+9) += 6 * (*(peftot_ip+7)) / (5 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+8*14+1) -= 6 * (*(peftot_ip+7)) / (5 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+1*14+8) -= 6 * (*(peftot_ip+7)) / (5 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+9*14+2) -= 6 * (*(peftot_ip+7)) / (5 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+2*14+9) -= 6 * (*(peftot_ip+7)) / (5 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+3*14+3) += 6 * (*(peftot_ip+7)) * (*(pipolar+n)) /
        (5 * (*(pcarea+NE_TR+n)) * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+10*14+10) += 6 * (*(peftot_ip+7)) * (*(pipolar+n)) /
        (5 * (*(pcarea+NE_TR+n)) * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+10*14+3) -= 6 * (*(peftot_ip+7)) * (*(pipolar+n)) /
        (5 * (*(pcarea+NE_TR+n)) * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+3*14+10) -= 6 * (*(peftot_ip+7)) * (*(pipolar+n)) /
        (5 * (*(pcarea+NE_TR+n)) * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+4*14+4) += 2 * (*(peftot_ip+7)) * (*(pdefllen_ip+NE_TR+n)) / 15;
    *(pk_fr+11*14+11) += 2 * (*(peftot_ip+7)) * (*(pdefllen_ip+NE_TR+n)) / 15;
    *(pk_fr+5*14+5) += 2 * (*(peftot_ip+7)) * (*(pdefllen_ip+NE_TR+n)) / 15;
    *(pk_fr+12*14+12) += 2 * (*(peftot_ip+7)) * (*(pdefllen_ip+NE_TR+n)) / 15;
    *(pk_fr+6*14+6) += 2 * (*(peftot_ip+7)) * (*(pipolar+n)) /
        (15 * (*(pcarea+NE_TR+n)));
    *(pk_fr+13*14+13) += 2 * (*(peftot_ip+7)) * (*(pipolar+n)) /
        (15 * (*(pcarea+NE_TR+n)));
    *(pk_fr+3*14+1) += (11 * (*(peftot_ip+4)) - *(peftot_ip+11)) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+1*14+3) += (11 * (*(peftot_ip+4)) - *(peftot_ip+11)) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+8*14+3) -= (11 * (*(peftot_ip+4)) - *(peftot_ip+11)) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+3*14+8) -= (11 * (*(peftot_ip+4)) - *(peftot_ip+11)) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+4*14+1) += *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+1*14+4) += *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+5*14+2) += *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+2*14+5) += *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+11*14+8) += *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+8*14+11) += *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+12*14+9) += *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+9*14+12) += *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+11*14+1) -= *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+1*14+11) -= *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+12*14+2) -= *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+2*14+12) -= *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+8*14+4) -= *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+4*14+8) -= *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+9*14+5) -= *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+5*14+9) -= *(peftot_ip+10) / *(pdefllen_ip+NE_TR+n);
    *(pk_fr+5*14+1) += *(peftot_ip+7) / 10;
    *(pk_fr+1*14+5) += *(peftot_ip+7) / 10;
    *(pk_fr+12*14+1) += *(peftot_ip+7) / 10;
    *(pk_fr+1*14+12) += *(peftot_ip+7) / 10;
    *(pk_fr+9*14+4) += *(peftot_ip+7) / 10;
    *(pk_fr+4*14+9) += *(peftot_ip+7) / 10;
    *(pk_fr+11*14+9) += *(peftot_ip+7) / 10;
    *(pk_fr+9*14+11) += *(peftot_ip+7) / 10;
    *(pk_fr+4*14+2) -= *(peftot_ip+7) / 10;
    *(pk_fr+2*14+4) -= *(peftot_ip+7) / 10;
    *(pk_fr+11*14+2) -= *(peftot_ip+7) / 10;
    *(pk_fr+2*14+11) -= *(peftot_ip+7) / 10;
    *(pk_fr+8*14+5) -= *(peftot_ip+7) / 10;
    *(pk_fr+5*14+8) -= *(peftot_ip+7) / 10;
    *(pk_fr+12*14+8) -= *(peftot_ip+7) / 10;
    *(pk_fr+8*14+12) -= *(peftot_ip+7) / 10;
    *(pk_fr+6*14+1) += *(peftot_ip+4) / 10;
    *(pk_fr+1*14+6) += *(peftot_ip+4) / 10;
    *(pk_fr+8*14+6) -= *(peftot_ip+4) / 10;
    *(pk_fr+6*14+8) -= *(peftot_ip+4) / 10;
    *(pk_fr+10*14+8) += (*(peftot_ip+4) - 11 * (*(peftot_ip+11))) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+8*14+10) += (*(peftot_ip+4) - 11 * (*(peftot_ip+11))) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+10*14+1) -= (*(peftot_ip+4) - 11 * (*(peftot_ip+11))) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+1*14+10) -= (*(peftot_ip+4) - 11 * (*(peftot_ip+11))) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+13*14+8) += *(peftot_ip+11) / 10;
    *(pk_fr+8*14+13) += *(peftot_ip+11) / 10;
    *(pk_fr+13*14+1) -= *(peftot_ip+11) / 10;
    *(pk_fr+1*14+13) -= *(peftot_ip+11) / 10;
    *(pk_fr+3*14+2) += (11 * (*(peftot_ip+5)) - *(peftot_ip+12)) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+2*14+3) += (11 * (*(peftot_ip+5)) - *(peftot_ip+12)) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+9*14+3) -= (11 * (*(peftot_ip+5)) - *(peftot_ip+12)) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+3*14+9) -= (11 * (*(peftot_ip+5)) - *(peftot_ip+12)) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+6*14+2) += *(peftot_ip+5) / 10;
    *(pk_fr+2*14+6) += *(peftot_ip+5) / 10;
    *(pk_fr+9*14+6) -= *(peftot_ip+5) / 10;
    *(pk_fr+6*14+9) -= *(peftot_ip+5) / 10;
    *(pk_fr+10*14+9) += (*(peftot_ip+5) - 11 * (*(peftot_ip+12))) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+9*14+10) += (*(peftot_ip+5) - 11 * (*(peftot_ip+12))) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+10*14+2) -= (*(peftot_ip+5) - 11 * (*(peftot_ip+12))) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+2*14+10) -= (*(peftot_ip+5) - 11 * (*(peftot_ip+12))) /
        (10 * (*(pdefllen_ip+NE_TR+n)));
    *(pk_fr+13*14+9) += *(peftot_ip+12) / 10;
    *(pk_fr+9*14+13) += *(peftot_ip+12) / 10;
    *(pk_fr+13*14+2) -= *(peftot_ip+12) / 10;
    *(pk_fr+2*14+13) -= *(peftot_ip+12) / 10;
    *(pk_fr+4*14+3) -= (2 * (*(peftot_ip+5)) - *(peftot_ip+12)) / 5;
    *(pk_fr+3*14+4) -= (2 * (*(peftot_ip+5)) - *(peftot_ip+12)) / 5;
    *(pk_fr+5*14+3) += (2 * (*(peftot_ip+4)) - *(peftot_ip+11)) / 5;
    *(pk_fr+3*14+5) += (2 * (*(peftot_ip+4)) - *(peftot_ip+11)) / 5;
    *(pk_fr+6*14+3) += *(peftot_ip+7) * (*(pipolar+n)) / (10 * (*(pcarea+NE_TR+n)));
    *(pk_fr+3*14+6) += *(peftot_ip+7) * (*(pipolar+n)) / (10 * (*(pcarea+NE_TR+n)));
    *(pk_fr+13*14+3) += *(peftot_ip+7) * (*(pipolar+n)) / (10 * (*(pcarea+NE_TR+n)));
    *(pk_fr+3*14+13) += *(peftot_ip+7) * (*(pipolar+n)) / (10 * (*(pcarea+NE_TR+n)));
    *(pk_fr+10*14+6) -= *(peftot_ip+7) * (*(pipolar+n)) / (10 * (*(pcarea+NE_TR+n)));
    *(pk_fr+6*14+10) -= *(peftot_ip+7) * (*(pipolar+n)) / (10 * (*(pcarea+NE_TR+n)));
    *(pk_fr+13*14+10) -= *(peftot_ip+7) * (*(pipolar+n)) / (10 * (*(pcarea+NE_TR+n)));
    *(pk_fr+10*14+13) -= *(peftot_ip+7) * (*(pipolar+n)) / (10 * (*(pcarea+NE_TR+n)));
    *(pk_fr+11*14+3) -= (2 * (*(peftot_ip+5)) + *(peftot_ip+12)) / 10;
    *(pk_fr+3*14+11) -= (2 * (*(peftot_ip+5)) + *(peftot_ip+12)) / 10;
    *(pk_fr+12*14+3) += (2 * (*(peftot_ip+4)) + *(peftot_ip+11)) / 10;
    *(pk_fr+3*14+12) += (2 * (*(peftot_ip+4)) + *(peftot_ip+11)) / 10;
    *(pk_fr+6*14+4) -= (3 * (*(peftot_ip+5)) - *(peftot_ip+12)) *
        (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+4*14+6) -= (3 * (*(peftot_ip+5)) - *(peftot_ip+12)) *
        (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+10*14+4) -= (*(peftot_ip+5) + 2 * (*(peftot_ip+12))) / 10;
    *(pk_fr+4*14+10) -= (*(peftot_ip+5) + 2 * (*(peftot_ip+12))) / 10;
    *(pk_fr+11*14+4) -= *(peftot_ip+7) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+4*14+11) -= *(peftot_ip+7) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+12*14+5) -= *(peftot_ip+7) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+5*14+12) -= *(peftot_ip+7) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+12*14+4) += *(peftot_ip+10) / 2;
    *(pk_fr+4*14+12) += *(peftot_ip+10) / 2;
    *(pk_fr+11*14+5) -= *(peftot_ip+10) / 2;
    *(pk_fr+5*14+11) -= *(peftot_ip+10) / 2;
    *(pk_fr+13*14+4) += *(peftot_ip+5) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+4*14+13) += *(peftot_ip+5) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+6*14+5) += (3 * (*(peftot_ip+4)) - *(peftot_ip+11)) *
        (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+5*14+6) += (3 * (*(peftot_ip+4)) - *(peftot_ip+11)) *
        (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+10*14+5) += (*(peftot_ip+4) + 2 * (*(peftot_ip+11))) / 10;
    *(pk_fr+5*14+10) += (*(peftot_ip+4) + 2 * (*(peftot_ip+11))) / 10;
    *(pk_fr+13*14+5) -= *(peftot_ip+4) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+5*14+13) -= *(peftot_ip+4) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+11*14+6) -= *(peftot_ip+12) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+6*14+11) -= *(peftot_ip+12) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+12*14+6) += *(peftot_ip+11) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+6*14+12) += *(peftot_ip+11) * (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+13*14+6) -= *(peftot_ip+7) * (*(pipolar+n)) / (30 * (*(pcarea+NE_TR+n)));
    *(pk_fr+6*14+13) -= *(peftot_ip+7) * (*(pipolar+n)) / (30 * (*(pcarea+NE_TR+n)));
    *(pk_fr+11*14+10) += (*(peftot_ip+5) - 2 * (*(peftot_ip+12))) / 5;
    *(pk_fr+10*14+11) += (*(peftot_ip+5) - 2 * (*(peftot_ip+12))) / 5;
    *(pk_fr+12*14+10) -= (*(peftot_ip+4) - 2 * (*(peftot_ip+11))) / 5;
    *(pk_fr+10*14+12) -= (*(peftot_ip+4) - 2 * (*(peftot_ip+11))) / 5;
    *(pk_fr+13*14+11) -= (*(peftot_ip+5) - 3 * (*(peftot_ip+12))) *
        (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+11*14+13) -= (*(peftot_ip+5) - 3 * (*(peftot_ip+12))) *
        (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+13*14+12) += (*(peftot_ip+4) - 3 * (*(peftot_ip+11))) *
        (*(pdefllen_ip+NE_TR+n)) / 30;
    *(pk_fr+12*14+13) += (*(peftot_ip+4) - 3 * (*(peftot_ip+11))) *
        (*(pdefllen_ip+NE_TR+n)) / 30;
}

void stiffm_fr (double *pk_fr, double *peftot_ip, int *pyldflag, double *pyield,
    double *pcarea, double *pzstrong, double *pzweak, long n)
{
    // Initialize function variables
    int i, j, k; // Counter variables

    // Yield surface parameters
    double Py = *(pcarea+NE_TR+n) * (*(pyield+NE_TR+n)); // Squash load
    double Mpy = *(pzweak+n) * (*(pyield+NE_TR+n)); // Weak-axis plastic moment
    double Mpz = *(pzstrong+n) * (*(pyield+NE_TR+n)); // Strong-axis plastic moment
    double p[2], my[2], mz[2]; // Ratios used in yield function
    double phi[2]; // Values of yield function at Member Ends 1 and 2

    double G[14][2]; // Matrix of yield surface gradients
    // Results of matrix multiplications for computation of plastic reduction matrix
    double k_G[14][2], GT_k[2][14], GT_k_G[2][2], kG_GTkG[14][2];
    double det, sum, temp; // Variables for matrix operations

    // Initialize all elements to zero
    for (i = 0; i < 14; ++i) {
        G[i][0] = G[i][1] = 0;
    }

    // Ratio of internal axial force to squash load at Member Ends 1 and 2
    p[0] = *(peftot_ip) / Py;
    p[1] = *(peftot_ip+7) / Py;
    /* Ratio of weak-axis internal moment to weak-axis plastic moment at Member Ends 1
       and 2 */
    my[0] = *(peftot_ip+4) / Mpy;
    my[1] = *(peftot_ip+11) / Mpy;
    /* Ratio of strong-axis internal moment to strong-axis plastic moment at Member Ends
       1 and 2 */
    mz[0] = *(peftot_ip+5) / Mpz;
    mz[1] = *(peftot_ip+12) / Mpz;

    // Compute values of yield function at Member Ends 1 and 2
    for (i = 0; i < 2; ++i) {
        phi[i] = pow(p[i],2) + pow(mz[i],2) + pow(my[i],4) +
            3.5 * pow(p[i],2) * pow(mz[i],2) + 3 * pow(p[i],6) * pow(my[i],2) +
            4.5 * pow(mz[i],4) * pow(my[i],2);
    }

    if (phi[0] >= 1 - phitol && phi[1] >= 1 - phitol) {
        // Assign non-zero elements of matrix of yield surface gradients
        if (*(pyldflag+n*2) != 2) {
            G[0][0] =
                2 * p[0] / Py + 7 * p[0] * pow(mz[0],2) / Py +
                18 * pow(p[0],5) * pow(my[0],2) / Py;
            G[4][0] =
                4 * pow(my[0],3) / Mpy + 6 * pow(p[0],6) * my[0] / Mpy +
                9 * pow(mz[0],4) * my[0] / Mpy;
            G[5][0] =
                2 * mz[0] / Mpz + 7 * pow(p[0],2) * mz[0] / Mpz +
                18 * pow(mz[0],3) * pow(my[0],2) / Mpz;
        }
        if (*(pyldflag+n*2+1) != 2) {
            G[7][1] =
                2 * p[1] / Py + 7 * p[1] * pow(mz[1],2) / Py +
                18 * pow(p[1],5) * pow(my[1],2) / Py;
            G[11][1] =
                4 * pow(my[1],3) / Mpy + 6 * pow(p[1],6) * my[1] / Mpy +
                9 * pow(mz[1],4) * my[1] / Mpy;
            G[12][1] =
                2 * mz[1] / Mpz + 7 * pow(p[1],2) * mz[1] / Mpz +
                18 * pow(mz[1],3) * pow(my[1],2) / Mpz;
        }

        /* (I) Multiply element tangent stiffness matrix by matrix of yield surface
           gradients and (II) multiply transpose of matrix of yield surface gradients by
           element tangent stiffness matrix, i.e. take transpose of (I) */
        for (i = 0; i < 14; ++i) {
            for (j = 0; j < 2; ++j) {
                sum = 0;
                for (k = 0; k < 14; ++k) {
                    sum += *(pk_fr+i*14+k) * G[k][j];
                }
                k_G[i][j] = sum;
                GT_k[j][i] = k_G[i][j];
            }
        }

        /* (III) Multiply above result (II) by matrix of yield surface gradients, i.e.
           carry out matrix multiplications within inverse */
        for (i = 0; i < 2; ++i) {
            for (j = 0; j < 2; ++j) {
                sum = 0;
                for (k = 0; k < 14; ++k) {
                    sum += GT_k[i][k] * G[k][j];
                }
                GT_k_G[i][j] = sum;
            }
        }

        // (IV) Compute the inverse of the above result (III)
        det = GT_k_G[0][0] * GT_k_G[1][1] - GT_k_G[0][1] * GT_k_G[1][0];
        temp = GT_k_G[0][0];
        GT_k_G[0][0] = GT_k_G[1][1] / det;
        GT_k_G[1][1] = temp / det;
        GT_k_G[0][1] *= -1 / det;
        GT_k_G[1][0] *= -1 / det;

        // (V) Multiply above result (I) by other above result (IV)
        for (i = 0; i < 14; ++i) {
            for (j = 0; j < 2; ++j) {
                sum = 0;
                for (k = 0; k < 2; ++k) {
                    sum += k_G[i][k] * GT_k_G[k][j];
                }
                kG_GTkG[i][j] = sum;
            }
        }

        // Multiply above result (V) by other above result (II)
        for (i = 0; i < 14; ++i) {
            for (j = 0; j < 14; ++j) {
                sum = 0;
                for (k = 0; k < 2; ++k) {
                    sum += kG_GTkG[i][k] * GT_k[k][j];
                }
                *(pk_fr+i*14+j) -= sum;
            }
        }
    } else if (phi[0] >= 1 - phitol && *(pyldflag+n*2) != 2) {
        // Assign non-zero elements of matrix of yield surface gradients
        G[0][0] =
            2 * p[0] / Py + 7 * p[0] * pow(mz[0],2) / Py +
            18 * pow(p[0],5) * pow(my[0],2) / Py;
        G[4][0] =
            4 * pow(my[0],3) / Mpy + 6 * pow(p[0],6) * my[0] / Mpy +
            9 * pow(mz[0],4) * my[0] / Mpy;
        G[5][0] =
            2 * mz[0] / Mpz + 7 * pow(p[0],2) * mz[0] / Mpz +
            18 * pow(mz[0],3) * pow(my[0],2) / Mpz;

        /* (I) Multiply element tangent stiffness matrix by matrix of yield surface
           gradients and (II) multiply transpose of matrix of yield surface gradients by
           element tangent stiffness matrix, i.e. take transpose of (I) */
        for (i = 0; i < 14; ++i) {
            sum = 0;
            for (k = 0; k < 14; ++k) {
                sum += *(pk_fr+i*14+k) * G[k][0];
            }
            k_G[i][0] = sum;
            GT_k[0][i] = k_G[i][0];
        }

        /* (III) Multiply above result (II) by matrix of yield surface gradients, i.e.
           carry out matrix multiplications within inverse */
        sum = 0;
        for (k = 0; k < 14; ++k) {
            sum += GT_k[0][k] * G[k][0];
        }
        GT_k_G[0][0] = sum;

        // (IV) Compute the inverse of the above result (III)
        GT_k_G[0][0] = 1 / GT_k_G[0][0];

        // (V) Multiply above result (I) by other above result (IV)
        for (i = 0; i < 14; ++i) {
            kG_GTkG[i][0] = k_G[i][0] * GT_k_G[0][0];
        }

        // Multiply above result (V) by other above result (II)
        for (i = 0; i < 14; ++i) {
            for (j = 0; j < 14; ++j) {
                *(pk_fr+i*14+j) -= kG_GTkG[i][0] * GT_k[0][j];
            }
        }
    } else if (phi[1] >= 1 - phitol && *(pyldflag+n*2+1) != 2) {
        // Assign non-zero elements of matrix of yield surface gradients
        G[7][0] =
            2 * p[1] / Py + 7 * p[1] * pow(mz[1],2) / Py +
            18 * pow(p[1],5) * pow(my[1],2) / Py;
        G[11][0] =
            4 * pow(my[1],3) / Mpy + 6 * pow(p[1],6) * my[1] / Mpy +
            9 * pow(mz[1],4) * my[1] / Mpy;
        G[12][0] =
            2 * mz[1] / Mpz + 7 * pow(p[1],2) * mz[1] / Mpz +
            18 * pow(mz[1],3) * pow(my[1],2) / Mpz;

        /* (I) Multiply element tangent stiffness matrix by matrix of yield surface
           gradients and (II) multiply transpose of matrix of yield surface gradients by
           element tangent stiffness matrix, i.e. take transpose of (I) */
        for (i = 0; i < 14; ++i) {
            sum = 0;
            for (k = 0; k < 14; ++k) {
                sum += *(pk_fr+i*14+k) * G[k][0];
            }
            k_G[i][0] = sum;
            GT_k[0][i] = k_G[i][0];
        }

        /* (III) Multiply above result (II) by matrix of yield surface gradients, i.e.
           carry out matrix multiplications within inverse */
        sum = 0;
        for (k = 0; k < 14; ++k) {
            sum += GT_k[0][k] * G[k][0];
        }
        GT_k_G[0][0] = sum;

        // (IV) Compute the inverse of the above result (III)
        GT_k_G[0][0] = 1 / GT_k_G[0][0];

        // (V) Multiply above result (I) by other above result (IV)
        for (i = 0; i < 14; ++i) {
            kG_GTkG[i][0] = k_G[i][0] * GT_k_G[0][0];
        }

        // Multiply above result (V) by other above result (II)
        for (i = 0; i < 14; ++i) {
            for (j = 0; j < 14; ++j) {
                *(pk_fr+i*14+j) -= kG_GTkG[i][0] * GT_k[0][j];
            }
        }
    }
}

void release (double *pk_fr, int *pmendrel, int n)
{
    // Initialize function variables
    int i, j, k, rel = 0;
    double det, sum, temp;

    // Determine number of member ends to be released
    for (i = 0; i < 4; ++i) {
        rel += *(pmendrel+n*5+1+i);
    }

    double G[14][rel]; // Member end bending release matrix
    /* Results of matrix multiplications for computation of member end bending release
       reduction matrix */
    double k_G[14][rel], GT_k[rel][14], GT_k_G[rel][rel], kG_GTkG[14][rel];

    // Initialize member end bending release matrix
    for (i = 0; i < 14; ++i) {
        for (j = 0; j < rel; ++j) {
            G[i][j] = 0;
        }
    }

    // Build member end bending release matrix
    j = 0;
    if (*(pmendrel+n*5+1) == 1) {
        G[5][j] = 1;
        j++;
    }
    if (*(pmendrel+n*5+2) == 1) {
        G[4][j] = 1;
        j++;
    }
    if (*(pmendrel+n*5+3) == 1) {
        G[12][j] = 1;
        j++;
    }
    if (*(pmendrel+n*5+4) == 1) {
        G[11][j] = 1;
    }

    /* (I) Multiply element tangent stiffness matrix by matrix of member end bending
       release matrix and (II) multiply transpose of member end bending release matrix by
       element tangent stiffness matrix, i.e. take transpose of (I) */
    for (i = 0; i < 14; ++i) {
        for (j = 0; j < rel; ++j) {
            sum = 0;
            for (k = 0; k < 14; ++k) {
                sum += *(pk_fr+i*14+k) * G[k][j];
            }
            k_G[i][j] = sum;
            GT_k[j][i] = k_G[i][j];
        }
    }

    /* (III) Multiply above result (II) by member end bending release matrix, i.e. carry
       out matrix multiplications within inverse */
    for (i = 0; i < rel; ++i) {
        for (j = 0; j < rel; ++j) {
            sum = 0;
            for (k = 0; k < 14; ++k) {
                sum += GT_k[i][k] * G[k][j];
            }
            GT_k_G[i][j] = sum;
        }
    }

    // (IV) Compute the inverse of the above result (III)
    if (rel == 1) {
        GT_k_G[0][0] = 1 / GT_k_G[0][0];
    } else if (rel == 2) {
        det = GT_k_G[0][0] * GT_k_G[1][1] - GT_k_G[0][1] * GT_k_G[1][0];
        temp = GT_k_G[0][0];
        GT_k_G[0][0] = GT_k_G[1][1] / det;
        GT_k_G[1][1] = temp / det;
        GT_k_G[0][1] *= -1 / det;
        GT_k_G[1][0] *= -1 / det;
    } else {
        inverse (&GT_k_G[0][0], rel);
    }

    // (V) Multiply above result (I) by other above result (IV)
    for (i = 0; i < 14; ++i) {
        for (j = 0; j < rel; ++j) {
            sum = 0;
            for (k = 0; k < rel; ++k) {
                sum += k_G[i][k] * GT_k_G[k][j];
            }
            kG_GTkG[i][j] = sum;
        }
    }

    // Multiply above result (V) by other above result (II)
    for (i = 0; i < 14; ++i) {
        for (j = 0; j < 14; ++j) {
            sum = 0;
            for (k = 0; k < rel; ++k) {
                sum += kG_GTkG[i][k] * GT_k[k][j];
            }
            *(pk_fr+i*14+j) -= sum;
        }
    }
}

int forces_fr (double *pf_temp, double *pef_ip, double *pef_i, double *pefFE_ref,
    double *pefFE_ip, double *pefFE_i, int *pyldflag, double *pdd, double *pemod,
    double *pgmod, double *pcarea, double *poffset, int *posflag, double *pllength,
    double *pdefllen_ip, double *pistrong, double *piweak, double *pipolar,
    double *piwarp, double *pyield, double *pzstrong, double *pzweak, double *pc1_ip,
    double *pc2_ip, double *pc3_ip, double *pc1_i, double *pc2_i, double *pc3_i,
    int *pmendrel, long *pmcode, double *pdlpf, int *pitecnt)
{
    // Initialize function variables
    int unlchk;
    long i, j, k, n, ptr, ptr2;
    double sum;

    /* General element stiffness matrix in local coordinate system; may serve as ke,
       ke + kg, ke - km, or ke + kg - km */
    double k_fr[14][14];

    /* Incremental element nodal displacements in global coordinate system, w.r.t. global
       joints 1 and 2 */
    double DD12[14];
    /* Incremental element nodal displacements in global coordinate system, w.r.t. local
       joints i and j */
    double DDij[14];
    // Incremental element nodal displacements in local coordinate system
    double dd[14];

    // Yield surface parameters
    double Py, Mpy, Mpz;
    double p[2], my[2], mz[2]; // Ratios used in yield surface equation
    double dp[2], dmy[2], dmz[2]; // Ratio increments used in yield surface equation
    double phi[2]; // Values of yield function at Member Ends 1 and 2

    double def[14]; // Incremental element force vector in local coordinate system
    // Total element force vector in previous local coordinate system
    double eftot_ip[14];
    double eftot_i[14]; // Total element force vector in current local coordinate system
    // Element force vector in global coordinate system, w.r.t. local joints i and j
    double EFij[14];

    // Transformation matrices
    double T_ip[14][14]; // From global to previous configuration
    double T_i[14][14]; // From global to current configuration
    double Ti_Tip[14][14]; // From previous configuration to current configuration
    double T_rl[14][14]; // Rigid link transformation matrix

    ptr = NE_TR * 2;
    ptr2 = NE_TR * 6;
    for (n = 0; n < NE_FR; ++n) {
        // Initialize all elements to zero and compute total element force vector
        for (i = 0; i < 14; ++i) {
            for (j = 0; j < 14; ++j) {
                k_fr[i][j] = 0;
                T_i[i][j] = T_ip[i][j] = T_rl[i][j] = 0;
            }
            eftot_ip[i] = *(pef_ip+ptr+n*14+i) + *(pefFE_ip+n*14+i);
        }

        // Pass control to stiffe_fr function
        stiffe_fr (&k_fr[0][0], pemod, pgmod, pcarea, pllength, pistrong, piweak,
            pipolar, piwarp, n);

        /* Include nonlinear subroutines in element stiffness matrix assembly depending
           upon user-requested analysis */
        if (ANAFLAG == 2) {
            // Pass control to stiffg_fr function
            stiffg_fr (&k_fr[0][0], eftot_ip, pdefllen_ip, pcarea, pipolar, n);
        } else if (ANAFLAG == 3) {
            // Pass control to stiffg_fr function
            stiffg_fr (&k_fr[0][0], eftot_ip, pdefllen_ip, pcarea, pipolar, n);

            if (*(pyldflag+n*2) != 2 || *(pyldflag+n*2+1) != 2) {
	            // Pass control to stiffm_fr function
	            stiffm_fr (&k_fr[0][0], eftot_ip, pyldflag, pyield, pcarea, pzstrong,
	                pzweak, n);
			}
        }

        // Check for member end bending releases
        if (*(pmendrel+n*5) == 1) {
            release (&k_fr[0][0], pmendrel, n);
        }

        /* Retrieve element incremental nodal displacements from generalized nodal
           displacement vector */
        for (i = 0; i < 14; ++i) {
            DD12[i] = 0;
            j = *(pmcode+ptr2+n*14+i);
            if (j != 0) {
                DD12[i] = *(pdd+j-1);
            }
        }

        // Assign non-zero elements of coordinate transformation matrix
        T_ip[0][0] = T_ip[3][3] = T_ip[7][7] = T_ip[10][10] = *(pc1_ip+NE_TR+n*3);
        T_ip[0][1] = T_ip[3][4] = T_ip[7][8] = T_ip[10][11] = *(pc1_ip+NE_TR+n*3+1);
        T_ip[0][2] = T_ip[3][5] = T_ip[7][9] = T_ip[10][12] = *(pc1_ip+NE_TR+n*3+2);
        T_ip[1][0] = T_ip[4][3] = T_ip[8][7] = T_ip[11][10] = *(pc2_ip+NE_TR+n*3);
        T_ip[1][1] = T_ip[4][4] = T_ip[8][8] = T_ip[11][11] = *(pc2_ip+NE_TR+n*3+1);
        T_ip[1][2] = T_ip[4][5] = T_ip[8][9] = T_ip[11][12] = *(pc2_ip+NE_TR+n*3+2);
        T_ip[2][0] = T_ip[5][3] = T_ip[9][7] = T_ip[12][10] = *(pc3_ip+NE_TR+n*3);
        T_ip[2][1] = T_ip[5][4] = T_ip[9][8] = T_ip[12][11] = *(pc3_ip+NE_TR+n*3+1);
        T_ip[2][2] = T_ip[5][5] = T_ip[9][9] = T_ip[12][12] = *(pc3_ip+NE_TR+n*3+2);
        T_ip[6][6] = T_ip[13][13] = 1;

        if (*(posflag+n) == 0) {
            /* Transform element incremental nodal displacement vector from global into
               local coordinate system */
            for (i = 0; i < 14; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += T_ip[i][j] * DD12[j];
                }
                dd[i] = sum;
            }
        } else {
            // Assign non-zero elements of rigid link transformation matrix
            for (i = 0; i < 14; ++i) {
                T_rl[i][i] = 1;
            }
            T_rl[3][1] = -(*(poffset+n*6+2));
            T_rl[3][2] = *(poffset+n*6+1);
            T_rl[4][0] = *(poffset+n*6+2);
            T_rl[4][2] = -(*(poffset+n*6));
            T_rl[5][0] = -(*(poffset+n*6+1));
            T_rl[5][1] = *(poffset+n*6);
            T_rl[10][8] = -(*(poffset+n*6+5));
            T_rl[10][9] = *(poffset+n*6+4);
            T_rl[11][7] = *(poffset+n*6+5);
            T_rl[11][9] = -(*(poffset+n*6+3));
            T_rl[12][7] = -(*(poffset+n*6+4));
            T_rl[12][8] = *(poffset+n*6+3);

            /* Transform element incremental nodal displacement vector from global joints
               1 and 2, to local joints i, j */
            for (i = 0; i < 14; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += T_rl[j][i] * DD12[j];
                }
                DDij[i] = sum;
            }

            /* Transform element incremental nodal displacement vector from global into
               local coordinate system */
            for (i = 0; i < 14; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += T_ip[i][j] * DDij[j];
                }
                dd[i] = sum;
            }
        }

        // Compute incremental element force vector
        for (i = 0; i < 14; ++i) {
            sum = 0;
            for (j = 0; j < 14; ++j) {
                sum += k_fr[i][j] * dd[j];
            }
            def[i] = sum;
        }

        // Assign non-zero elements of coordinate transformation matrix
        T_i[0][0] = T_i[3][3] = T_i[7][7] = T_i[10][10] = *(pc1_i+NE_TR+n*3);
        T_i[0][1] = T_i[3][4] = T_i[7][8] = T_i[10][11] = *(pc1_i+NE_TR+n*3+1);
        T_i[0][2] = T_i[3][5] = T_i[7][9] = T_i[10][12] = *(pc1_i+NE_TR+n*3+2);
        T_i[1][0] = T_i[4][3] = T_i[8][7] = T_i[11][10] = *(pc2_i+NE_TR+n*3);
        T_i[1][1] = T_i[4][4] = T_i[8][8] = T_i[11][11] = *(pc2_i+NE_TR+n*3+1);
        T_i[1][2] = T_i[4][5] = T_i[8][9] = T_i[11][12] = *(pc2_i+NE_TR+n*3+2);
        T_i[2][0] = T_i[5][3] = T_i[9][7] = T_i[12][10] = *(pc3_i+NE_TR+n*3);
        T_i[2][1] = T_i[5][4] = T_i[9][8] = T_i[12][11] = *(pc3_i+NE_TR+n*3+1);
        T_i[2][2] = T_i[5][5] = T_i[9][9] = T_i[12][12] = *(pc3_i+NE_TR+n*3+2);
        T_i[6][6] = T_i[13][13] = 1;

        /* Construct coordinate transformation matrix which transforms previous
           configuration to current configuration, i.e. T_i * T_ip^T */
        for (i = 0; i < 14; ++i) {
            for (j = 0; j < 14; ++j) {
                sum = 0;
                for (k = 0; k < 14; ++k) {
                    sum += T_i[i][k] * T_ip[j][k];
                }
                Ti_Tip[i][j] = sum;
            }
        }

        /* Add contribution of incremental element force vector to previous element force
           vector and update reference configuration */
        for (i = 0; i < 14; ++i) {
            sum = 0;
            for (j = 0; j < 14; ++j) {
                sum += Ti_Tip[i][j] * (*(pef_ip+ptr+n*14+j) + def[j]);
            }
            *(pef_i+ptr+n*14+i) = eftot_i[i] = sum;
        }

        /* Add contribution of increment in cumulative element fixed-end force vector to
           previous element fixed-end force vector and update reference configuration */
        if (*pitecnt == 0 && *(pyldflag+n*2) == 0 && *(pyldflag+n*2+1) == 0) {
            for (i = 0; i < 14; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += Ti_Tip[i][j] * (*(pefFE_ip+n*14+j) +
                        *pdlpf * (*(pefFE_ref+n*14+j)));
                }
                *(pefFE_i+n*14+i) = sum;
                eftot_i[i] = *(pef_i+ptr+n*14+i) + sum;
            }
        } else if (*pitecnt == 0 && *(pyldflag+n*2) == 0) {
            for (i = 0; i < 7; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += Ti_Tip[i][j] * (*(pefFE_ip+n*14+j) +
                        *pdlpf * (*(pefFE_ref+n*14+j)));
                }
                *(pefFE_i+n*14+i) = sum;
                eftot_i[i] = *(pef_i+ptr+n*14+i) + sum;
            }
            for (i = 7; i < 14; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += Ti_Tip[i][j] * (*(pefFE_ip+n*14+j));
                }
                *(pefFE_i+n*14+i) = sum;
                eftot_i[i] = *(pef_i+ptr+n*14+i) + sum;
            }
        } else if (*pitecnt == 0 && *(pyldflag+n*2+1) == 0) {
            for (i = 0; i < 7; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += Ti_Tip[i][j] * (*(pefFE_ip+n*14+j));
                }
                *(pefFE_i+n*14+i) = sum;
                eftot_i[i] = *(pef_i+ptr+n*14+i) + sum;
            }
            for (i = 7; i < 14; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += Ti_Tip[i][j] * (*(pefFE_ip+n*14+j) +
                        *pdlpf * (*(pefFE_ref+n*14+j)));
                }
                *(pefFE_i+n*14+i) = sum;
                eftot_i[i] = *(pef_i+ptr+n*14+i) + sum;
            }
        } else {
            for (i = 0; i < 14; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += Ti_Tip[i][j] * (*(pefFE_ip+n*14+j));
                }
                *(pefFE_i+n*14+i) = sum;
                eftot_i[i] = *(pef_i+ptr+n*14+i) + sum;
            }
        }

        if (ANAFLAG == 3) {
            Py = *(pcarea+NE_TR+n) * (*(pyield+NE_TR+n)); // Squash load
            Mpy = *(pzweak+n) * (*(pyield+NE_TR+n)); // Weak-axis plastic moment
            Mpz = *(pzstrong+n) * (*(pyield+NE_TR+n)); // Strong-axis plastic moment

            // Ratio of internal axial force to squash load at Member Ends 1 and 2
            p[0] = eftot_i[0] / Py;
            p[1] = eftot_i[7] / Py;
            /* Ratio of weak-axis internal moment to weak-axis plastic moment at Member
               Ends 1 and 2 */
            my[0] = eftot_i[4] / Mpy;
            my[1] = eftot_i[11] / Mpy;
            /* Ratio of strong-axis internal moment to strong-axis plastic moment at
               Member Ends 1 and 2 */
            mz[0] = eftot_i[5] / Mpz;
            mz[1] = eftot_i[12] / Mpz;

            // Compute values of yield function at Member Ends 1 and 2
            for (i = 0; i < 2; ++i) {
                phi[i] =
                    pow(p[i],2) + pow(mz[i],2) + pow(my[i],4) +
                    3.5 * pow(p[i],2) * pow(mz[i],2) +
                    3 * pow(p[i],6) * pow(my[i],2) +
                    4.5 * pow(mz[i],4) * pow(my[i],2);
            }

            // Check whether phi is on or exceeding the yield surface
            if (phi[0] > phi[1] && phi[0] > 1 + phitol && *(pyldflag+n*2) != 2) {
                // Ratio of internal axial force to squash load at Member End-1
                p[0] = eftot_ip[0] / Py;
                dp[0] = (eftot_i[0] - eftot_ip[0]) / Py;
                /* Ratio of weak-axis internal moment to weak-axis plastic moment at
                   Member End-1 */
                my[0] = eftot_ip[4] / Mpy;
                dmy[0] = (eftot_i[4] - eftot_ip[4]) / Mpy;
                /* Ratio of strong-axis internal moment to strong-axis plastic moment at
                   Member End-1 */
                mz[0] = eftot_ip[5] / Mpz;
                dmz[0] = (eftot_i[5] - eftot_ip[5]) / Mpz;

                /* Compute updated dlpf such that corrected load increment results in
                   "phi >= 1 - phitol" and "phi <= 1 + phitol" */
                *pdlpf *= regula_falsi(&p[0], &dp[0], &my[0], &dmy[0], &mz[0], &dmz[0]);
                *(pyldflag+n*2) = 1;
                return 1;
            } else if (phi[1] > phi[0] && phi[1] > 1 + phitol &&
                *(pyldflag+n*2+1) != 2) {
                // Ratio of internal axial force to squash load at Member End-2
                p[1] = eftot_ip[7] / Py;
                dp[1] = (eftot_i[7] - eftot_ip[7]) / Py;
                /* Ratio of weak-axis internal moment to weak-axis plastic moment at
                   Member End-2 */
                my[1] = eftot_ip[11] / Mpy;
                dmy[1] = (eftot_i[11] - eftot_ip[11]) / Mpy;
                /* Ratio of strong-axis internal moment to strong-axis plastic moment at
                   Member End-2 */
                mz[1] = eftot_ip[12] / Mpz;
                dmz[1] = (eftot_i[12] - eftot_ip[12]) / Mpz;

                /* Compute updated dlpf such that corrected load increment results in
                   "phi >= 1 - phitol" and "phi <= 1 + phitol" */
                *pdlpf *= regula_falsi(&p[1], &dp[1], &my[1], &dmy[1], &mz[1], &dmz[1]);
                *(pyldflag+n*2+1) = 1;
                return 1;
            } else if ((phi[0] >= 1 - phitol && *(pyldflag+n*2) != 2) &&
                (phi[1] >= 1 - phitol && *(pyldflag+n*2+1) != 2)) {
                *(pyldflag+n*2) = *(pyldflag+n*2+1) = 1;
            } else if (phi[0] >= 1 - phitol && *(pyldflag+n*2) != 2) {
                *(pyldflag+n*2) = 1;
            } else if (phi[1] >= 1 - phitol && *(pyldflag+n*2+1) != 2) {
                *(pyldflag+n*2+1) = 1;
            }

            // Check for elastic unloading of yielded member ends
            if (*(pyldflag+n*2) == 1 || *(pyldflag+n*2+1) == 1) {
		        // Initialize all elements to zero
		        for (i = 0; i < 14; ++i) {
		            for (j = 0; j < 14; ++j) {
		                k_fr[i][j] = 0;
		            }
		        }

		        // Pass control to stiffe_fr function
		        stiffe_fr (&k_fr[0][0], pemod, pgmod, pcarea, pllength, pistrong, piweak,
		            pipolar, piwarp, n);

	            // Pass control to stiffg_fr function
	            stiffg_fr (&k_fr[0][0], eftot_ip, pdefllen_ip, pcarea, pipolar, n);

	            // Ratio of internal axial force to squash load at Member Ends 1 and 2
	            p[0] = eftot_i[0] / Py;
	            p[1] = eftot_i[7] / Py;
	            /* Ratio of weak-axis internal moment to weak-axis plastic moment at
	               Member Ends 1 and 2 */
	            my[0] = eftot_i[4] / Mpy;
	            my[1] = eftot_i[11] / Mpy;
	            /* Ratio of strong-axis internal moment to strong-axis plastic moment at
	               Member Ends 1 and 2 */
	            mz[0] = eftot_i[5] / Mpz;
	            mz[1] = eftot_i[12] / Mpz;

                // Pass control to unload function
                unlchk = unload(phi, p, my, mz, &Py, &Mpy, &Mpz, &k_fr[0][0], dd, n);
                if (unlchk == 1) {
                    *(pyldflag+n*2) = *(pyldflag+n*2+1) = 2;
                    return 2;
                } else if (unlchk == 2) {
                    *(pyldflag+n*2) = 2;
                    return 2;
                } else if (unlchk == 3) {
                    *(pyldflag+n*2+1) = 2;
                    return 2;
                }
            }
        }

        if (*(posflag+n) == 0) {
            /* Transform element force vector from local into global coordinate system
               and add element contribution to generalized internal force vector */
            for (i = 0; i < 14; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += T_i[j][i] * (*(pef_i+ptr+n*14+j));
                }
                j = *(pmcode+ptr2+n*14+i);
                if (j != 0) {
                    *(pf_temp+j-1) += sum;
                }
            }
        } else {
            /* Transform element force vector from local into global coordinate system
               w.r.t. local joints i and j */
            for (i = 0; i < 14; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += T_i[j][i] * (*(pef_i+ptr+n*14+j));
                }
                EFij[i] = sum;
            }

            /* Transform element force vector to global joints 1 and 2 and add element
               contribution to generalized internal force vector */
            for (i = 0; i < 14; ++i) {
                sum = 0;
                for (j = 0; j < 14; ++j) {
                    sum += T_rl[i][j] * EFij[j];
                }
                j = *(pmcode+ptr2+n*14+i);
                if (j != 0) {
                    *(pf_temp+j-1) += sum;
                }
            }
        }
    }
    return 0;
}

void mass_fr (double *psm, double *pcarea, double *pllength, double *pistrong, double *piweak, 
              double *pipolar, double *piwarp, double *pdens, int *posflag, double *poffset, 
              double *px,  double *pxfr, long *pminc, long *pmcode, double *pjac)
{
	long i, j, k, l,ie, je, ptr;
	double el[3];
    
    ptr = NE_TR * 2;
	double m_fr[14][14]; // General element mass matrix 
	
	for (i = 0; i < NE_FR; ++i) {
		
        // Compute element length
        j = *(pminc+ptr+i*2) - 1;
        k = *(pminc+ptr+i*2+1) - 1;
        if (*(posflag+i) == 0) {
            for (l = 0; l < 3; ++l) {
                *(pxfr+i*6+l) = *(px+j*3+l);
                *(pxfr+i*6+3+l) = *(px+k*3+l);
            }
        } else {
            for (l = 0; l < 3; ++l) {
                *(pxfr+i*6+l) = *(px+j*3+l) + *(poffset+i*6+l);
                *(pxfr+i*6+3+l) = *(px+k*3+l) + *(poffset+i*6+3+l);
            }
        }
        el[0] = *(pxfr+i*6+3) - *(pxfr+i*6);
        el[1] = *(pxfr+i*6+3+1) - *(pxfr+i*6+1);
        el[2] = *(pxfr+i*6+3+2) - *(pxfr+i*6+2);
        *(pllength+NE_TR+i) = sqrt(dot(el,el,3));
        
		// Initialize element mass array to zero
		for (j = 0; j < 14; ++j) {
			for (k = 0; k < 14; k++) {
				m_fr[j][k] = 0;
			}
		}
        
        /* Assemble lumped element mass matrix (Bathe Table 9.5) */
        
        // Displacement degrees of freedom
        m_fr[0][0] = m_fr[1][1] = m_fr[2][2] = m_fr[7][7] = m_fr[8][8] = m_fr[9][9] = 
        (*(pdens+i) * (*(pcarea+i)) * (*(pllength+i)))/24 * 12;
        
        // Rotational degrees of freedom
        m_fr[3][3] = m_fr[4][4] = m_fr[5][5] = m_fr[10][10] = m_fr[11][11] = m_fr[12][12] = 
        (*(pdens+i) * (*(pcarea+i)) * (*(pllength+i)))/24 * (pow(*(pllength+i),2));

		
		// Assemble system mass array - size [lss] - for use in skyline solver or system mass matrix - full order [NEQ][NEQ]
        if (SLVFLAG == 0) {
            /* Initialize index and then assign element mass components to structure mass array by index and mcode */
            for (ie = 0; ie < 14; ++ie) {
                for (je = 0; je < 14; ++je) {
                    j = *(pmcode+i*14+ie);
                    k = *(pmcode+i*14+je);
                    /* Add current element mass component to previous elements' components to the given DOFs */
                    if (j != 0) {
                        if (j == k) {
                            *(psm+j-1) += m_fr[ie][je];
                        }
                    }
                }
            }
        }
        else {
            /*Build the full order (i.e. [NEQ][NEQ] mass mastrix using mcode*/
            for (ie = 0; ie < 14; ++ie) {
                for (je = 0; je < 14; ++je) {
                    j = *(pmcode+i*14+ie);
                    k = *(pmcode+i*14+je);
                    
                    if ((j != 0) && (k != 0)) {
                        *(psm+(j-1)*NEQ+k-1) += m_fr[je][ie];
                    }
                }
            }
        }
    }
}



double regula_falsi (double *pp, double *pdp, double *pmy, double *pdmy, double *pmz,
    double *pdmz)
{
    // Initialize function variables
    double tau_u, tau_l, tau_r;
    double ptot, mytot, mztot;
    double phi_u, phi_l, phi_r;

    /* Guess at upper and lower bounds on root and compute associated upper and lower
       bounds on phi */
    // Upper bound of root conservatively taken to be 1, for initial guess
    tau_u = 1;
    ptot = *pp + tau_u * (*pdp);
    mytot = *pmy + tau_u * (*pdmy);
    mztot = *pmz + tau_u * (*pdmz);
    phi_u =
        pow(ptot,2) + pow(mztot,2) + pow(mytot,4) + 3.5 * pow(ptot,2) * pow(mztot,2) +
        3 * pow(ptot,6) * pow(mytot,2) + 4.5 * pow(mztot,4) * pow(mytot,2);
    // Lower bound of root conservatively taken to be 0, for initial guess
    tau_l = 0;
    ptot = *pp + tau_l * (*pdp);
    mytot = *pmy + tau_l * (*pdmy);
    mztot = *pmz + tau_l * (*pdmz);
    phi_l =
        pow(ptot,2) + pow(mztot,2) + pow(mytot,4) + 3.5 * pow(ptot,2) * pow(mztot,2) +
        3 * pow(ptot,6) * pow(mytot,2) + 4.5 * pow(mztot,4) * pow(mytot,2);

    // Compute estimate of the scalar multiplier
    tau_r = tau_u - (phi_u - 1) * (tau_l - tau_u) / (phi_l - phi_u);
    ptot = *pp + tau_r * (*pdp);
    mytot = *pmy + tau_r * (*pdmy);
    mztot = *pmz + tau_r * (*pdmz);
    phi_r =
        pow(ptot,2) + pow(mztot,2) + pow(mytot,4) + 3.5 * pow(ptot,2) * pow(mztot,2) +
        3 * pow(ptot,6) * pow(mytot,2) + 4.5 * pow(mztot,4) * pow(mytot,2);

    // Iterate on the scalar multiplier such that phi_r is within the specified tolerance
    do {
        if ((phi_l - 1 > 0 && phi_r - 1 > 0) || (phi_l - 1 < 0 && phi_r - 1 < 0)) {
            tau_l = tau_r;
            phi_l = phi_r;
        } else {
            tau_u = tau_r;
            phi_u = phi_r;
        }
        tau_r = tau_u - (phi_u - 1) * (tau_l - tau_u) / (phi_l - phi_u);
        ptot = *pp + tau_r * (*pdp);
        mytot = *pmy + tau_r * (*pdmy);
        mztot = *pmz + tau_r * (*pdmz);
        phi_r =
            pow(ptot,2) + pow(mztot,2) + pow(mytot,4) +
            3.5 * pow(ptot,2) * pow(mztot,2) +
            3 * pow(ptot,6) * pow(mytot,2) +
            4.5 * pow(mztot,4) * pow(mytot,2);

    } while (phi_r <= 1 - phitol && phi_r >= 1 + phitol);

    return tau_r;
}

int unload (double *pphi, double *pp, double *pmy, double *pmz, double *pPy,
    double *pMpy, double *pMpz, double *pk_fr, double *pdd, long n)
{
    // Initialize function variables
    int i, j, k; // Counter variables

    double G[14][2]; // Matrix of yield surface gradients
    /* Results of matrix multiplications for computation of magnitudes of plastic
       deformations */
    double GT_k[2][14], GT_k_G[2][2], GTkG_GTk[2][14];
    double lambda[2]; // Magnitude of plastic deformation
    double det, sum, temp; // Variables for matrix operations

    // Initialize all elements to zero
    for (i = 0; i < 14; ++i) {
        G[i][0] = G[i][1] = 0;
    }

    if (*pphi >= 1 - phitol && *(pphi+1) >= 1 - phitol) {
        // Assign non-zero elements of matrix of yield surface gradients
        G[0][0] =
            2 * (*pp) / (*pPy) + 7 * (*pp) * pow((*pmz),2) / (*pPy) +
            18 * pow((*pp),5) * pow((*pmy),2) / (*pPy);
        G[4][0] =
            4 * pow((*pmy),3) / (*pMpy) + 6 * pow((*pp),6) * (*pmy) / (*pMpy) +
            9 * pow((*pmz),4) * (*pmy) / (*pMpy);
        G[5][0] =
            2 * (*pmz) / (*pMpz) + 7 * pow((*pp),2) * (*pmz) / (*pMpz) +
            18 * pow((*pmz),3) * pow((*pmy),2) / (*pMpz);
        G[7][1] =
            2 * (*(pp+1)) / (*pPy) + 7 * (*(pp+1)) * pow((*(pmz+1)),2) / (*pPy) +
            18 * pow((*(pp+1)),5) * pow((*(pmy+1)),2) / (*pPy);
        G[11][1] =
            4 * pow((*(pmy+1)),3) / (*pMpy) +
            6 * pow((*(pp+1)),6) * (*(pmy+1)) / (*pMpy) +
            9 * pow((*(pmz+1)),4) * (*(pmy+1)) / (*pMpy);
        G[12][1] =
            2 * (*(pmz+1)) / (*pMpz) + 7 * pow((*(pp+1)),2) * (*(pmz+1)) / (*pMpz) +
            18 * pow((*(pmz+1)),3) * pow((*(pmy+1)),2) / (*pMpz);

        /* (I) Multiply element tangent stiffness matrix by matrix of yield surface
           gradients and (II) multiply transpose of matrix of yield surface gradients by
           element tangent stiffness matrix, i.e. take transpose of (I) */
        for (i = 0; i < 2; ++i) {
            for (j = 0; j < 14; ++j) {
                sum = 0;
                for (k = 0; k < 14; ++k) {
                    sum += G[k][i] * (*(pk_fr+k*14+j));
                }
                GT_k[i][j] = sum;
            }
        }

        /* (III) Multiply above result (II) by matrix of yield surface gradients, i.e.
           carry out matrix multiplications within inverse */
        for (i = 0; i < 2; ++i) {
            for (j = 0; j < 2; ++j) {
                sum = 0;
                for (k = 0; k < 14; ++k) {
                    sum += GT_k[i][k] * G[k][j];
                }
                GT_k_G[i][j] = sum;
            }
        }

        // (IV) Compute the inverse of the above result (III)
        det = GT_k_G[0][0] * GT_k_G[1][1] - GT_k_G[0][1] * GT_k_G[1][0];
        temp = GT_k_G[0][0];
        GT_k_G[0][0] = GT_k_G[1][1] / det;
        GT_k_G[1][1] = temp / det;
        GT_k_G[0][1] *= -1 / det;
        GT_k_G[1][0] *= -1 / det;

        // (V) Multiply above result (I) by other above result (IV)
        for (i = 0; i < 2; ++i) {
            for (j = 0; j < 14; ++j) {
                sum = 0;
                for (k = 0; k < 2; ++k) {
                    sum += GT_k_G[i][k] * GT_k[k][j];
                }
                GTkG_GTk[i][j] = sum;
            }
        }

        /* Multiply above result (V) by incremental displacements to obtain magnitude of
           plastic deformation */
        for (i = 0; i < 2; ++i) {
            sum = 0;
            for (j = 0; j < 14; ++j) {
                sum += GTkG_GTk[i][j] * (*(pdd+j));
            }
            lambda[i] = sum;
        }

        if (lambda[0] < -1e-8 && lambda[1] < -1e-8) {
            return 1;
        } else if (lambda[0] < -1e-8) {
            return 2;
        } else if (lambda[1] < -1e-8) {
            return 3;
        }
    } else if (*pphi >= 1 - phitol) {
        // Assign non-zero elements of matrix of yield surface gradients
        G[0][0] =
            2 * (*pp) / (*pPy) + 7 * (*pp) * pow((*pmz),2) / (*pPy) +
            18 * pow((*pp),5) * pow((*pmy),2) / (*pPy);
        G[4][0] =
            4 * pow((*pmy),3) / (*pMpy) + 6 * pow((*pp),6) * (*pmy) / (*pMpy) +
            9 * pow((*pmz),4) * (*pmy) / (*pMpy);
        G[5][0] =
            2 * (*pmz) / (*pMpz) + 7 * pow((*pp),2) * (*pmz) / (*pMpz) +
            18 * pow((*pmz),3) * pow((*pmy),2) / (*pMpz);

        /* (I) Multiply element tangent stiffness matrix by matrix of yield surface
           gradients and (II) multiply transpose of matrix of yield surface gradients by
           element tangent stiffness matrix, i.e. take transpose of (I) */
        for (j = 0; j < 14; ++j) {
            sum = 0;
            for (k = 0; k < 14; ++k) {
                sum += G[k][0] * (*(pk_fr+k*14+j));
            }
            GT_k[0][j] = sum;
        }

        /* (III) Multiply above result (II) by matrix of yield surface gradients, i.e.
           carry out matrix multiplications within inverse */
        sum = 0;
        for (k = 0; k < 14; ++k) {
            sum += GT_k[0][k] * G[k][0];
        }
        GT_k_G[0][0] = sum;

        // (IV) Compute the inverse of the above result (III)
        GT_k_G[0][0] = 1 / GT_k_G[0][0];

        // (V) Multiply above result (I) by other above result (IV)
        for (j = 0; j < 14; ++j) {
            GTkG_GTk[0][j] = GT_k_G[0][0] * GT_k[0][j];
        }

        /* Multiply above result (V) by incremental displacements to obtain magnitude of
           plastic deformation */
        sum = 0;
        for (j = 0; j < 14; ++j) {
            sum += GTkG_GTk[0][j] * (*(pdd+j));
        }
        if (sum < -1e-8) {
            return 2;
        }
    } else if (*(pphi+1) >= 1 - phitol) {
        // Assign non-zero elements of matrix of yield surface gradients
        G[7][0] =
            2 * (*(pp+1)) / (*pPy) + 7 * (*(pp+1)) * pow((*(pmz+1)),2) / (*pPy) +
            18 * pow((*(pp+1)),5) * pow((*(pmy+1)),2) / (*pPy);
        G[11][0] =
            4 * pow((*(pmy+1)),3) / (*pMpy) +
            6 * pow((*(pp+1)),6) * (*(pmy+1)) / (*pMpy) +
            9 * pow((*(pmz+1)),4) * (*(pmy+1)) / (*pMpy);
        G[12][0] =
            2 * (*(pmz+1)) / (*pMpz) + 7 * pow((*(pp+1)),2) * (*(pmz+1)) / (*pMpz) +
            18 * pow((*(pmz+1)),3) * pow((*(pmy+1)),2) / (*pMpz);

        /* (I) Multiply element tangent stiffness matrix by matrix of yield surface
           gradients and (II) multiply transpose of matrix of yield surface gradients by
           element tangent stiffness matrix, i.e. take transpose of (I) */
        for (j = 0; j < 14; ++j) {
            sum = 0;
            for (k = 0; k < 14; ++k) {
                sum += G[k][0] * (*(pk_fr+k*14+j));
            }
            GT_k[0][j] = sum;
        }

        /* (III) Multiply above result (II) by matrix of yield surface gradients, i.e.
           carry out matrix multiplications within inverse */
        sum = 0;
        for (k = 0; k < 14; ++k) {
            sum += GT_k[0][k] * G[k][0];
        }
        GT_k_G[0][0] = sum;

        // (IV) Compute the inverse of the above result (III)
        GT_k_G[0][0] = 1 / GT_k_G[0][0];

        // (V) Multiply above result (I) by other above result (IV)
        for (j = 0; j < 14; ++j) {
            GTkG_GTk[0][j] = GT_k_G[0][0] * GT_k[0][j];
        }

        /* Multiply above result (V) by incremental displacements to obtain magnitude of
           plastic deformation */
        sum = 0;
        for (j = 0; j < 14; ++j) {
            sum += GTkG_GTk[0][j] * (*(pdd+j));
        }
        if (sum < -1e-8) {
            return 3;
        }
    }

    return 0;
}
