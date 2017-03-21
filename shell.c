//********************************************************************************
//**																			**
//**  Pertains to CU-BEN ver 3.141												**
//**																			**
//**  Copyright (c) 2017 C. J. Earls                                            **
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

extern long NJ, NE_TR, NE_FR, NE_SH, NEQ;
extern int ANAFLAG, ALGFLAG, OPTFLAG, SLVFLAG;
extern FILE *IFP[2], *OFP[5];

void prop_sh (double *px, double *pemod, double *pnu, double *pxlocal, double *pthick, double *pdens,
    double *pfarea, double *pslength, double *pyield, double *pc1, double *pc2,
    double *pc3, long *pminc)
{
    // Initialize function variables
    long i, j, k, l, m, ptr, ptr2, ptr3;
    double el12[3], el23[3], el31[3], normal[3], localx[3], localy[3], localz[3];

    fprintf(OFP[0], "\nShell Element Properties:\n\tElement\t\tElastic Modulus\t\t");
    fprintf(OFP[0], "Poisson's Ratio\t\tThickness\tDensity\t\tArea\t\tLength-12\tLength-23\t");
    fprintf(OFP[0], "Length-31\tYield Stress\n");

    ptr = NE_TR + NE_FR;
    ptr2 = NE_TR * 2 + NE_FR * 2;
    ptr3 = NE_TR + NE_FR * 3;
    for (i = 0; i < NE_SH; ++i) {
        
        // Read in element properties from input file
        fscanf(IFP[0], "%lf,%lf,%lf,%lf,%lf\n", pemod+ptr+i, pnu+i, pthick+i, pdens+i, pyield+ptr+i);

        // Compute element side-lengths
        j = *(pminc+ptr2+i*3) - 1;
        k = *(pminc+ptr2+i*3+1) - 1;
        l = *(pminc+ptr2+i*3+2) - 1;
        /* Note that the order is reversed to be from Vertex-1 to Vertex-3 so that localz
           is properly oriented */
        for (m = 0; m < 3; ++m) {
            el23[m] = *(px+l*3+m) - *(px+k*3+m);
            el31[m] = *(px+l*3+m) - *(px+j*3+m);
            el12[m] = *(px+k*3+m) - *(px+j*3+m);
        }
        *(pslength+i*3+1) = sqrt(dot(el23,el23,3));
        *(pslength+i*3+2) = sqrt(dot(el31,el31,3));
        *(pslength+i*3) = sqrt(dot(el12,el12,3));

        // Compute element face area
        cross(el12,el31,normal,0); // Compute vector normal to shell surface
        *(pfarea+i) = 0.5 * sqrt(dot(normal,normal,3));

        // Compute direction cosines
        for (m = 0; m < 3; ++m) {
            localx[m] = el12[m] / *(pslength+i*3);
            localz[m] = normal[m] / (2 * (*(pfarea+i)));
        }
        cross(localz,localx,localy,1);
        for (m = 0; m < 3; ++m) {
            *(pc1+ptr3+i*3+m) = localx[m];
            *(pc2+ptr3+i*3+m) = localy[m];
            *(pc3+ptr3+i*3+m) = localz[m];
        }

        // Pass control to the mem_coord function
        mem_coord (pxlocal, i, i, j, k, l, px, pc1, pc2, pc3, ptr3);

        // Write out element properties to output file
        fprintf(OFP[0], "\t%ld\t\t%lf\t\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", i + 1,
            *(pemod+ptr+i), *(pnu+i), *(pthick+i), *(pdens+i), *(pfarea+i), *(pslength+i*3),
            *(pslength+i*3+1), *(pslength+i*3+2), *(pyield+ptr+i));
    }
    if (OPTFLAG == 2) {
    	for (i = 0; i < NE_SH; ++i) {
    		fprintf(IFP[1], "%e,%e,%e,%e\n", *(pemod+ptr+i), *(pnu+i),
				*(pthick+i), *(pyield+ptr+i));
    	}
    }
}

void stiff_sh (double *pss, double *pemod, double *pnu, double *px_temp, double *pxlocal,
    double *pthick, double *pfarea, double *pdeffarea_ip, double *pslength,
    double *pdefslen_ip, double *pyield, double *pc1_ip, double *pc2_ip, double *pc3_ip,
    double *pef_ip, double *pd_temp, double *pchi_temp, double *pefN_temp,
    double *pefM_temp, long *pmaxa, long *pminc, long *pmcode)
{
    // Initialize function variables
    long i, j, k, n, je, ie, ptr, ptr2, ptr3, ptr4;

    double k_sh[18][18]; // General element stiffness matrix in global coordinate system
    double T_ip[18][18]; // Coordinate transformation matrix
    double K_sh[18][18]; // Total element stiffness matrix in global coordinate system

    // Yield surface parameters
    double dm[6]; // Element membrane nodal displacements in local coordinate system
    // Ivanov's yield criteria parameters
    double No; // Uniaxial yield force per unit width
    double alpha[3]; // Partial plastification factor
    double Me[3]; // Modified uniaxial yield moment per unit width
    double Nbar[3], Mbar[3], MNbar[3]; // Quadratic stress intensities
    // Factors for computation of Ivanov's yield criteria
    double q_fact[3], r_fact[3], s_fact[3];
    int h_fact[3];
    double phi[3]; // Values of yield function at Vertices 1, 2, and 3
    int yv; // Yielded vertex; "0" indicates no yielded vertex

    ptr = NE_TR + NE_FR;
    ptr2 = NE_TR * 2 + NE_FR * 2;
    ptr3 = NE_TR + NE_FR * 3;
    ptr4 = NE_TR * 6 + NE_FR * 14;
    for (n = 0; n < NE_SH; ++n) {
        // Initialize all elements to zero
        for (i = 0; i < 18; ++i) {
            for (j = 0; j < 18; ++j) {
                k_sh[i][j] = 0;
                T_ip[i][j] = 0;
            }
        }

        // Compute element stiffness matrix depending upon user-requested analysis
        if (ANAFLAG == 1 || ANAFLAG == 4) {
            // Pass control to stiffe_sh function
            stiffe_sh (&k_sh[0][0], pemod, pnu, pxlocal, pthick, pfarea, pslength, ptr,
                n);
        } else if (ANAFLAG == 2) {
            // Pass control to stiffe_sh function
            stiffe_sh (&k_sh[0][0], pemod, pnu, pxlocal, pthick, pfarea, pslength, ptr,
                n);

            // Pass control to mem_coord function
            mem_coord (dm, n, 0, *(pminc+ptr2+n*3) - 1, *(pminc+ptr2+n*3+1) - 1,
                *(pminc+ptr2+n*3+2) - 1, px_temp, pc1_ip, pc2_ip, pc3_ip, ptr3);

            // Assign element membrane nodal displacements
            dm[5] = dm[2] - *(pxlocal+n*3+2);
            dm[4] = dm[1] - *(pxlocal+n*3+1);
            dm[2] = dm[0] - *(pxlocal+n*3);
            dm[0] = dm[1] = dm[3] = 0;

            // Pass control to stiffg_sh function
            stiffg_sh (&k_sh[0][0], pemod, pnu, pxlocal, pthick, pdeffarea_ip, dm, ptr,
                n);
        } else {
            // Compute uniaxial yield force per unit length
            No = *(pyield+ptr+n) * (*(pthick+n));
            yv = 0; // Reset yielded vertex flag

            // Compute Ivanov's yield criteria at Vertices 1, 2, and 3
            for (i = 0; i < 3; ++i) {
                // Compute modified uniaxial yield moment per unit width
                alpha[i] = 1.0 - 0.4 * exp(-2.6 * sqrt(*(pchi_temp+n*3+i)));
                Me[i] = alpha[i] * 0.25 * (*(pyield+ptr+n)) * pow(*(pthick+n),2);

                // Compuate quadratic stress intensities
                Nbar[i] = pow(*(pefN_temp+n*9+i*3),2) + pow(*(pefN_temp+n*9+i*3+1),2) -
                    *(pefN_temp+n*9+i*3) * (*(pefN_temp+n*9+i*3+1)) +
                    3 * pow(*(pefN_temp+n*9+i*3+2),2);
                Mbar[i] = pow(*(pefM_temp+n*9+i*3),2) + pow(*(pefM_temp+n*9+i*3+1),2) -
                    *(pefM_temp+n*9+i*3) * (*(pefM_temp+n*9+i*3+1)) +
                    3 * pow(*(pefM_temp+n*9+i*3+2),2);
                MNbar[i] = *(pefM_temp+n*9+i*3) * (*(pefN_temp+n*9+i*3)) +
                    *(pefM_temp+n*9+i*3+1) * (*(pefN_temp+n*9+i*3+1)) -
                    0.5 * (*(pefM_temp+n*9+i*3)) * (*(pefN_temp+n*9+i*3+1)) -
                    0.5 * (*(pefM_temp+n*9+i*3+1)) * (*(pefN_temp+n*9+i*3)) +
                    3 * (*(pefM_temp+n*9+i*3+2)) * (*(pefN_temp+n*9+i*3+2));

                // Compute factors for computation of Ivanov's yield criteria
                q_fact[i] = Nbar[i] * pow(Me[i],2) + 0.48 * Mbar[i] * pow(No,2);
                if (q_fact[i] >= 1e-4) {
                    r_fact[i] = sqrt(pow(No,2) * pow(Mbar[i],2) +
                        4 * pow(Me[i],2) * pow(MNbar[i],2));
                    if (r_fact[i] / (2 * pow(Me[i],2) * No) >= 1e-4) {
                        h_fact[i] = 1;
                    } else {
                        h_fact[i] = 0;
                    }
                    s_fact[i] = Nbar[i] * Mbar[i] - pow(MNbar[i],2);

                    /* Compute Ivanov's yield criteria and if value is greater than one,
                       flag vertex as yielded */
                    if (h_fact[i] == 1) {
                        phi[i] = Nbar[i] / pow(No,2) + 0.5 * Mbar[i] / pow(Me[i],2) -
                            0.25 * s_fact[i] / q_fact[i] +
                            r_fact[i] / (2 * pow(Me[i],2) * No);
                        if (phi[i] >= 1 - phitol) {
                            if (yv == 0) {
                                yv = i + 1;
                            } else {
                                if (phi[i] < phi[yv - 1]) {
                                    yv = i + 1;
                                }
                            }
                        }
                    } else {
                        phi[i] = Nbar[i] / pow(No,2) + 0.5 * Mbar[i] / pow(Me[i],2) -
                            0.25 * s_fact[i] / q_fact[i];
                        if (phi[i] >= 1 - phitol) {
                            if (yv == 0) {
                                yv = i + 1;
                            } else {
                                if (phi[i] < phi[yv - 1]) {
                                    yv = i + 1;
                                }
                            }
                        }
                    }
                } else {
                    phi[i] = 0;
                }
            }

            /* Check if yielding has occured at any vertex and compute appropriate
               element stiffness matrix */
            if (yv == 0) {
                // Pass control to stiffe_sh function
                stiffe_sh (&k_sh[0][0], pemod, pnu, pxlocal, pthick, pfarea, pslength,
                    ptr, n);

                // Pass control to mem_coord function
                mem_coord (dm, n, 0, *(pminc+ptr2+n*3) - 1, *(pminc+ptr2+n*3+1) - 1,
                    *(pminc+ptr2+n*3+2) - 1, px_temp, pc1_ip, pc2_ip, pc3_ip, ptr3);

                // Assign element membrane nodal displacements
                dm[5] = dm[2] - *(pxlocal+n*3+2);
                dm[4] = dm[1] - *(pxlocal+n*3+1);
                dm[2] = dm[0] - *(pxlocal+n*3);
                dm[0] = dm[1] = dm[3] = 0;

                // Pass control to stiffg_sh function
                stiffg_sh (&k_sh[0][0], pemod, pnu, pxlocal, pthick, pdeffarea_ip, dm,
                    ptr, n);
            } else {
                // Pass control to stiffm_sh function
                stiffm_sh (&k_sh[0][0], pemod, pnu, pxlocal, pthick, pdeffarea_ip,
                    pdefslen_ip, pyield, pefN_temp, pefM_temp, pchi_temp, yv-1, &No,
                    &alpha[yv-1], &Me[yv-1], &Nbar[yv-1], &Mbar[yv-1], &MNbar[yv-1],
                    &q_fact[yv-1], &r_fact[yv-1], &s_fact[yv-1], &h_fact[yv-1], ptr, n);

                // Pass control to mem_coord function
                mem_coord (dm, n, 0, *(pminc+ptr2+n*3) - 1, *(pminc+ptr2+n*3+1) - 1,
                    *(pminc+ptr2+n*3+2) - 1, px_temp, pc1_ip, pc2_ip, pc3_ip, ptr3);

                // Assign element membrane nodal displacements
                dm[5] = dm[2] - *(pxlocal+n*3+2);
                dm[4] = dm[1] - *(pxlocal+n*3+1);
                dm[2] = dm[0] - *(pxlocal+n*3);
                dm[0] = dm[1] = dm[3] = 0;

                // Pass control to stiffg_sh function
                stiffg_sh (&k_sh[0][0], pemod, pnu, pxlocal, pthick, pdeffarea_ip, dm,
                    ptr, n);
            }
        }

        // Assign non-zero elements of coordinate transformation matrix
        T_ip[0][0] = T_ip[3][3] = T_ip[6][6] = T_ip[9][9] = T_ip[12][12] =
            T_ip[15][15] = *(pc1_ip+ptr3+n*3);
        T_ip[0][1] = T_ip[3][4] = T_ip[6][7] = T_ip[9][10] = T_ip[12][13] =
            T_ip[15][16] = *(pc1_ip+ptr3+n*3+1);
        T_ip[0][2] = T_ip[3][5] = T_ip[6][8] = T_ip[9][11] = T_ip[12][14] =
            T_ip[15][17] = *(pc1_ip+ptr3+n*3+2);
        T_ip[1][0] = T_ip[4][3] = T_ip[7][6] = T_ip[10][9] = T_ip[13][12] =
            T_ip[16][15] = *(pc2_ip+ptr3+n*3);
        T_ip[1][1] = T_ip[4][4] = T_ip[7][7] = T_ip[10][10] = T_ip[13][13] =
            T_ip[16][16] = *(pc2_ip+ptr3+n*3+1);
        T_ip[1][2] = T_ip[4][5] = T_ip[7][8] = T_ip[10][11] = T_ip[13][14] =
            T_ip[16][17] = *(pc2_ip+ptr3+n*3+2);
        T_ip[2][0] = T_ip[5][3] = T_ip[8][6] = T_ip[11][9] = T_ip[14][12] =
            T_ip[17][15] = *(pc3_ip+ptr3+n*3);
        T_ip[2][1] = T_ip[5][4] = T_ip[8][7] = T_ip[11][10] = T_ip[14][13] =
            T_ip[17][16] = *(pc3_ip+ptr3+n*3+1);
        T_ip[2][2] = T_ip[5][5] = T_ip[8][8] = T_ip[11][11] = T_ip[14][14] =
            T_ip[17][17] = *(pc3_ip+ptr3+n*3+2);

        // Pass control to transform function
        transform (&k_sh[0][0], &T_ip[0][0], &K_sh[0][0], 18);
		
		if (SLVFLAG == 0) {
			/* Initialize index and then assign element tangent stiffness coefficients of
			 element n to the structure stiffness matrix by index, mcode, and maxa */
            for (je = 0; je < 18; ++je) {
				j = *(pmcode+ptr4+n*18+je);
				if (j != 0) {
					// Check mcode above current entry to find rank of "j"
					for (ie = 0; ie <= je; ++ie) {
						i = *(pmcode+ptr4+n*18+ie);
						if (i != 0) {
							if (i > j) { // Find element address as diagonal address + delta
								k = *(pmaxa+i-1) + (i - j);
							} else {
								k = *(pmaxa+j-1) + (j - i);
							}
							/* Add current element stiffness to previous elements'
							 contributions to the given DOFs */
							*(pss+k-1) += K_sh[ie][je];
						}
					}
				}
			}
		}
		else {
			/* Build the full order (i.e. [NEQ][NEQ] system stiffness matrix using mcode */
			for (ie = 0; ie < 18; ++ie) {
				for (je = 0; je < 18; ++je) {
                    i = *(pmcode+ptr4+n*18+ie);
					j = *(pmcode+ptr4+n*18+je);

					if ((i != 0) && (j != 0)) {
						*(pss+(i-1)*NEQ+j-1) += K_sh[je][ie];
					}
				}
			}
		}       
    }
}

void stiffe_sh (double *pk_sh, double *pemod, double *pnu, double *pxlocal,
    double *pthick, double *pfarea, double *pslength, long ptr, long n)
{
    // Initialize function variables
    double ke_m_sh[6][6], ke_b_sh[9][9]; // Membrane and plate bending stiffness matrices

    // Pass control to stiff_m_sh function
    stiffe_m_sh (&ke_m_sh[0][0], pemod, pnu, pxlocal, pthick, pfarea, ptr, n);

    /* Add contribution of linear-elastic membrane stiffness matrix to element stiffness
       matrix */
    *(pk_sh+0*18+0)   = ke_m_sh[0][0];
    *(pk_sh+0*18+1)   = ke_m_sh[0][1];
    *(pk_sh+1*18+0)   = ke_m_sh[1][0];
    *(pk_sh+1*18+1)   = ke_m_sh[1][1];
    *(pk_sh+0*18+6)   = ke_m_sh[0][2];
    *(pk_sh+0*18+7)   = ke_m_sh[0][3];
    *(pk_sh+1*18+6)   = ke_m_sh[1][2];
    *(pk_sh+1*18+7)   = ke_m_sh[1][3];
    *(pk_sh+0*18+12)  = ke_m_sh[0][4];
    *(pk_sh+0*18+13)  = ke_m_sh[0][5];
    *(pk_sh+1*18+12)  = ke_m_sh[1][4];
    *(pk_sh+1*18+13)  = ke_m_sh[1][5];
    *(pk_sh+6*18+0)   = ke_m_sh[2][0];
    *(pk_sh+6*18+1)   = ke_m_sh[2][1];
    *(pk_sh+7*18+0)   = ke_m_sh[3][0];
    *(pk_sh+7*18+1)   = ke_m_sh[3][1];
    *(pk_sh+6*18+6)   = ke_m_sh[2][2];
    *(pk_sh+6*18+7)   = ke_m_sh[2][3];
    *(pk_sh+7*18+6)   = ke_m_sh[3][2];
    *(pk_sh+7*18+7)   = ke_m_sh[3][3];
    *(pk_sh+6*18+12)  = ke_m_sh[2][4];
    *(pk_sh+6*18+13)  = ke_m_sh[2][5];
    *(pk_sh+7*18+12)  = ke_m_sh[3][4];
    *(pk_sh+7*18+13)  = ke_m_sh[3][5];
    *(pk_sh+12*18+0)  = ke_m_sh[4][0];
    *(pk_sh+12*18+1)  = ke_m_sh[4][1];
    *(pk_sh+13*18+0)  = ke_m_sh[5][0];
    *(pk_sh+13*18+1)  = ke_m_sh[5][1];
    *(pk_sh+12*18+6)  = ke_m_sh[4][2];
    *(pk_sh+12*18+7)  = ke_m_sh[4][3];
    *(pk_sh+13*18+6)  = ke_m_sh[5][2];
    *(pk_sh+13*18+7)  = ke_m_sh[5][3];
    *(pk_sh+12*18+12) = ke_m_sh[4][4];
    *(pk_sh+12*18+13) = ke_m_sh[4][5];
    *(pk_sh+13*18+12) = ke_m_sh[5][4];
    *(pk_sh+13*18+13) = ke_m_sh[5][5];

    // Pass control to stiff_b_sh function
    stiffe_b_sh (&ke_b_sh[0][0], pemod, pnu, pxlocal, pthick, pfarea, pslength, ptr, n);

    /* Add contribution of linear-elastic plate bending stiffness matrix to element
       stiffness matrix */
    *(pk_sh+2*18+2)   = ke_b_sh[0][0];
    *(pk_sh+2*18+3)   = ke_b_sh[0][1];
    *(pk_sh+2*18+4)   = ke_b_sh[0][2];
    *(pk_sh+3*18+2)   = ke_b_sh[1][0];
    *(pk_sh+3*18+3)   = ke_b_sh[1][1];
    *(pk_sh+3*18+4)   = ke_b_sh[1][2];
    *(pk_sh+4*18+2)   = ke_b_sh[2][0];
    *(pk_sh+4*18+3)   = ke_b_sh[2][1];
    *(pk_sh+4*18+4)   = ke_b_sh[2][2];
    *(pk_sh+2*18+8)   = ke_b_sh[0][3];
    *(pk_sh+2*18+9)   = ke_b_sh[0][4];
    *(pk_sh+2*18+10)  = ke_b_sh[0][5];
    *(pk_sh+3*18+8)   = ke_b_sh[1][3];
    *(pk_sh+3*18+9)   = ke_b_sh[1][4];
    *(pk_sh+3*18+10)  = ke_b_sh[1][5];
    *(pk_sh+4*18+8)   = ke_b_sh[2][3];
    *(pk_sh+4*18+9)   = ke_b_sh[2][4];
    *(pk_sh+4*18+10)  = ke_b_sh[2][5];
    *(pk_sh+2*18+14)  = ke_b_sh[0][6];
    *(pk_sh+2*18+15)  = ke_b_sh[0][7];
    *(pk_sh+2*18+16)  = ke_b_sh[0][8];
    *(pk_sh+3*18+14)  = ke_b_sh[1][6];
    *(pk_sh+3*18+15)  = ke_b_sh[1][7];
    *(pk_sh+3*18+16)  = ke_b_sh[1][8];
    *(pk_sh+4*18+14)  = ke_b_sh[2][6];
    *(pk_sh+4*18+15)  = ke_b_sh[2][7];
    *(pk_sh+4*18+16)  = ke_b_sh[2][8];
    *(pk_sh+8*18+2)   = ke_b_sh[3][0];
    *(pk_sh+8*18+3)   = ke_b_sh[3][1];
    *(pk_sh+8*18+4)   = ke_b_sh[3][2];
    *(pk_sh+9*18+2)   = ke_b_sh[4][0];
    *(pk_sh+9*18+3)   = ke_b_sh[4][1];
    *(pk_sh+9*18+4)   = ke_b_sh[4][2];
    *(pk_sh+10*18+2)  = ke_b_sh[5][0];
    *(pk_sh+10*18+3)  = ke_b_sh[5][1];
    *(pk_sh+10*18+4)  = ke_b_sh[5][2];
    *(pk_sh+8*18+8)   = ke_b_sh[3][3];
    *(pk_sh+8*18+9)   = ke_b_sh[3][4];
    *(pk_sh+8*18+10)  = ke_b_sh[3][5];
    *(pk_sh+9*18+8)   = ke_b_sh[4][3];
    *(pk_sh+9*18+9)   = ke_b_sh[4][4];
    *(pk_sh+9*18+10)  = ke_b_sh[4][5];
    *(pk_sh+10*18+8)  = ke_b_sh[5][3];
    *(pk_sh+10*18+9)  = ke_b_sh[5][4];
    *(pk_sh+10*18+10) = ke_b_sh[5][5];
    *(pk_sh+8*18+14)  = ke_b_sh[3][6];
    *(pk_sh+8*18+15)  = ke_b_sh[3][7];
    *(pk_sh+8*18+16)  = ke_b_sh[3][8];
    *(pk_sh+9*18+14)  = ke_b_sh[4][6];
    *(pk_sh+9*18+15)  = ke_b_sh[4][7];
    *(pk_sh+9*18+16)  = ke_b_sh[4][8];
    *(pk_sh+10*18+14) = ke_b_sh[5][6];
    *(pk_sh+10*18+15) = ke_b_sh[5][7];
    *(pk_sh+10*18+16) = ke_b_sh[5][8];
    *(pk_sh+14*18+2)  = ke_b_sh[6][0];
    *(pk_sh+14*18+3)  = ke_b_sh[6][1];
    *(pk_sh+14*18+4)  = ke_b_sh[6][2];
    *(pk_sh+15*18+2)  = ke_b_sh[7][0];
    *(pk_sh+15*18+3)  = ke_b_sh[7][1];
    *(pk_sh+15*18+4)  = ke_b_sh[7][2];
    *(pk_sh+16*18+2)  = ke_b_sh[8][0];
    *(pk_sh+16*18+3)  = ke_b_sh[8][1];
    *(pk_sh+16*18+4)  = ke_b_sh[8][2];
    *(pk_sh+14*18+8)  = ke_b_sh[6][3];
    *(pk_sh+14*18+9)  = ke_b_sh[6][4];
    *(pk_sh+14*18+10) = ke_b_sh[6][5];
    *(pk_sh+15*18+8)  = ke_b_sh[7][3];
    *(pk_sh+15*18+9)  = ke_b_sh[7][4];
    *(pk_sh+15*18+10) = ke_b_sh[7][5];
    *(pk_sh+16*18+8)  = ke_b_sh[8][3];
    *(pk_sh+16*18+9)  = ke_b_sh[8][4];
    *(pk_sh+16*18+10) = ke_b_sh[8][5];
    *(pk_sh+14*18+14) = ke_b_sh[6][6];
    *(pk_sh+14*18+15) = ke_b_sh[6][7];
    *(pk_sh+14*18+16) = ke_b_sh[6][8];
    *(pk_sh+15*18+14) = ke_b_sh[7][6];
    *(pk_sh+15*18+15) = ke_b_sh[7][7];
    *(pk_sh+15*18+16) = ke_b_sh[7][8];
    *(pk_sh+16*18+14) = ke_b_sh[8][6];
    *(pk_sh+16*18+15) = ke_b_sh[8][7];
    *(pk_sh+16*18+16) = ke_b_sh[8][8];

    // Assign in-plane rotation stiffness to element stiffness matrix
    *(pk_sh+5*18+5)   = ke_b_sh[1][1] / 10000;
    *(pk_sh+11*18+11) = ke_b_sh[4][4] / 10000;
    *(pk_sh+17*18+17) = ke_b_sh[7][7] / 10000;
}

void stiffe_m_sh (double *pke_m_sh, double *pemod, double *pnu, double *pxlocal,
    double *pthick, double *pfarea, long ptr, long n)
{
    // Initialize function variables
    int i, j, k;
    double sum;
    double Bm[3][6]; // Membrane strain-displacement matrix
    double C[3][3]; // Plane stress constitutive matrix
    double Bm_C[6][3];

    // Define plane stress constitutive matrix
    C[0][2] = C[1][2] = C[2][0] = C[2][1] = 0;
    C[0][0] = C[1][1] = *(pemod+ptr+n) / (1 - pow(*(pnu+n),2));
    C[0][1] = C[1][0] = *(pemod+ptr+n) / (1 - pow(*(pnu+n),2)) * (*(pnu+n));
    C[2][2] = *(pemod+ptr+n) / (1 - pow(*(pnu+n),2)) * (1 - *(pnu+n)) / 2;

    // Compute membrane strain-displacement matrix
    Bm[0][1] = Bm[0][3] = Bm[0][4] = Bm[0][5] = Bm[1][0] = Bm[1][2] = Bm[1][4] =
        Bm[2][5] = 0;
    Bm[0][0] = Bm[2][1] = -(*(pxlocal+n*3+2) / (2 * (*(pfarea+n))));
    Bm[0][2] = Bm[2][3] = *(pxlocal+n*3+2) / (2 * (*(pfarea+n)));
    Bm[1][1] = Bm[2][0] = (*(pxlocal+n*3+1) - *(pxlocal+n*3)) / (2 * (*(pfarea+n)));
    Bm[1][3] = Bm[2][2] = -(*(pxlocal+n*3+1) / (2 * (*(pfarea+n))));
    Bm[1][5] = Bm[2][4] = *(pxlocal+n*3) / (2 * (*(pfarea+n)));

    // Compute linear-elastic membrane element stiffness matrix
    for (i = 0; i < 6; ++i) {
        for (j = 0; j < 3; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += Bm[k][i] * C[k][j];
            }
            Bm_C[i][j] = sum;
        }
    }
    for (i = 0; i < 6; ++i) {
        for (j = 0; j < 6; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += Bm_C[i][k] * Bm[k][j];
            }
            *(pke_m_sh+i*6+j) = *(pthick+n) * (*(pfarea+n)) * sum;
        }
    }
}

void stiffe_b_sh (double *pke_b_sh, double *pemod, double *pnu, double *pxlocal,
    double *pthick, double *pfarea, double *pslength, long ptr, long n)
{
    // Initialize function variables
    int i, j, k;
    double sum;
    // Coefficients for computing plate bending strain-displacement matrix
    double x23, l12, l23, l31, p4, p5, p6, t4, t5, q4, q5, r4, r5;
    double E1, E2, E3, E4; // Plate bending "constitutive matrix" coefficients
    double Q[9][9], b1, b2, b3;

    // Define plate bending "constitutive matrix"
    E1 = E3 = *(pemod+ptr+n) * pow(*(pthick+n),3) / (12 * (1 - pow(*(pnu+n),2)));
    E2 = *(pemod+ptr+n) * pow(*(pthick+n),3) / (12 * (1 - pow(*(pnu+n),2))) * (*(pnu+n));
    E4 = *(pemod+ptr+n) * pow(*(pthick+n),3) /
        (12 * (1 - pow(*(pnu+n),2))) * (1 - *(pnu+n)) / 2;

    // Define local coordinate coefficients for convenience in defining p, t, q, and r
    x23 = *(pxlocal+n*3) - *(pxlocal+n*3+1);
    l12 = pow(*(pslength+n*3),2);
    l23 = pow(*(pslength+n*3+1),2);
    l31 = pow(*(pslength+n*3+2),2);

    // Compute coefficients for assembly of plate bending strain-displacement matrix
    p4 = -6 * x23 / l23;
    p5 = -6 * (*(pxlocal+n*3+1)) / l31;
    p6 = 6 * (*(pxlocal+n*3)) / l12;
    t4 = 6 * (*(pxlocal+n*3+2)) / l23;
    t5 = -6 * (*(pxlocal+n*3+2)) / l31;
    q4 = -3 * x23 * (*(pxlocal+n*3+2)) / l23;
    q5 = 3 * (*(pxlocal+n*3+1)) * (*(pxlocal+n*3+2)) / l31;
    r4 = 3 * pow(*(pxlocal+n*3+2),2) / l23;
    r5 = 3 * pow(*(pxlocal+n*3+2),2) / l31;

    /* Compute the TRANSPOSE of the alpha matrix; computing the transpose eliminates a
       nested loop operation */
    double alpha_T[9][9] =
    {
        {*(pxlocal+n*3+2) * p6, -(*(pxlocal+n*3+2) * p6), *(pxlocal+n*3+2) * p5,
            -(*(pxlocal+n*3) * t5), 0, x23 * t5,
            -(*(pxlocal+n*3+1) * p6) - *(pxlocal+n*3) * p5, -x23 * p6,
            x23 * p5 + *(pxlocal+n*3+2) * t5},
        {0, 0, -(*(pxlocal+n*3+2) * q5), x23 + *(pxlocal+n*3) * r5, x23, x23 * (1 - r5),
            *(pxlocal+n*3) * q5 + *(pxlocal+n*3+2), *(pxlocal+n*3+2),
            -x23 * q5 + *(pxlocal+n*3+2) * (1 - r5)},
        {-4 * (*(pxlocal+n*3+2)), 2 * (*(pxlocal+n*3+2)), *(pxlocal+n*3+2) * (2 - r5),
            -(*(pxlocal+n*3) * q5), 0, x23 * q5, -4 * x23 + *(pxlocal+n*3) * r5, 2 * x23,
            x23 * (2 - r5) + *(pxlocal+n*3+2) * q5},
        {-(*(pxlocal+n*3+2) * p6), *(pxlocal+n*3+2) * p6, *(pxlocal+n*3+2) * p4, 0,
            *(pxlocal+n*3) * t4, -(*(pxlocal+n*3+1) * t4), *(pxlocal+n*3+1) * p6,
            x23 * p6 + *(pxlocal+n*3) * p4,
            -(*(pxlocal+n*3+1) * p4) + *(pxlocal+n*3+2) * t4},
        {0, 0, *(pxlocal+n*3+2) * q4, *(pxlocal+n*3+1),
            *(pxlocal+n*3+1) + *(pxlocal+n*3) * r4, *(pxlocal+n*3+1) * (1 - r4),
            -(*(pxlocal+n*3+2)), -(*(pxlocal+n*3+2)) + *(pxlocal+n*3) * q4,
            *(pxlocal+n*3+2) * (r4 - 1) - *(pxlocal+n*3+1) * q4},
        {-2 * (*(pxlocal+n*3+2)), 4 * (*(pxlocal+n*3+2)), *(pxlocal+n*3+2) * (r4 - 2), 0,
            -(*(pxlocal+n*3) * q4), *(pxlocal+n*3+1) * q4, 2 * (*(pxlocal+n*3+1)),
            -4 * (*(pxlocal+n*3+1)) + *(pxlocal+n*3) * r4,
            *(pxlocal+n*3+1) * (2 - r4) - *(pxlocal+n*3+2) * q4},
        {0, 0, -(*(pxlocal+n*3+2) * (p4 + p5)), *(pxlocal+n*3) * t5,
            -(*(pxlocal+n*3) * t4), -x23 * t5 + *(pxlocal+n*3+1) * t4,
            *(pxlocal+n*3) * p5, -(*(pxlocal+n*3) * p4),
            -x23 * p5 + *(pxlocal+n*3+1) * p4 - *(pxlocal+n*3+2) * (t4 + t5)},
        {0, 0, *(pxlocal+n*3+2) * (q4 - q5), *(pxlocal+n*3) * (r5 - 1),
            *(pxlocal+n*3) * (r4 - 1),
            -x23 * r5 - *(pxlocal+n*3+1) * r4 - *(pxlocal+n*3), *(pxlocal+n*3) * q5,
            *(pxlocal+n*3) * q4,
            -x23 * q5 - *(pxlocal+n*3+1) * q4 + *(pxlocal+n*3+2) * (r4 - r5)},
        {0, 0, *(pxlocal+n*3+2) * (r4 - r5), -(*(pxlocal+n*3) * q5),
            -(*(pxlocal+n*3) * q4), *(pxlocal+n*3+1) * q4 + x23 * q5,
            *(pxlocal+n*3) * (r5 - 2), *(pxlocal+n*3) * (r4 - 2),
            -x23 * r5 - *(pxlocal+n*3+1) * r4 + 4 * (*(pxlocal+n*3)) +
            *(pxlocal+n*3+2) * (q5 - q4)}
    };

    // Compute linear-elastic plate bending element stiffness matrix
    for (i = 0; i < 3; ++i) {
        b1 = b2 = b3 = 0;
        for (j = 0; j < 3; ++j) {
            b1 += E1 * alpha_T[i][j] + E2 * alpha_T[i][j + 3];
            b2 += E2 * alpha_T[i][j] + E3 * alpha_T[i][j + 3];
            b3 += E4 * alpha_T[i][j + 6];
        }
        for (j = 0; j < 3; ++j) {
            Q[i][j] = (E1 * alpha_T[i][j] + E2 * alpha_T[i][j + 3] + b1) / 24;
            Q[i][j + 3] = (E2 * alpha_T[i][j] + E3 * alpha_T[i][j + 3] + b2) / 24;
            Q[i][j + 6] = (E4 * alpha_T[i][j + 6] + b3) / 24;
        }
    }
    for (i = 3; i < 6; ++i) {
        b1 = b2 = b3 = 0;
        for (j = 3; j < 6; ++j) {
            b1 += E1 * alpha_T[i][j - 3] + E2 * alpha_T[i][j];
            b2 += E2 * alpha_T[i][j - 3] + E3 * alpha_T[i][j];
            b3 += E4 * alpha_T[i][j + 3];
        }
        for (j = 3; j < 6; ++j) {
            Q[i][j - 3] = (E1 * alpha_T[i][j - 3] + E2 * alpha_T[i][j] + b1) / 24;
            Q[i][j] = (E2 * alpha_T[i][j - 3] + E3 * alpha_T[i][j] + b2) / 24;
            Q[i][j + 3] = (E4 * alpha_T[i][j + 3] + b3) / 24;
        }
    }
    for (i = 6; i < 9; ++i) {
        b1 = b2 = b3 = 0;
        for (j = 6; j < 9; ++j) {
            b1 += E1 * alpha_T[i][j - 6] + E2 * alpha_T[i][j - 3];
            b2 += E2 * alpha_T[i][j - 6] + E3 * alpha_T[i][j - 3];
            b3 += E4 * alpha_T[i][j];
        }
        for (j = 6; j < 9; ++j) {
            Q[i][j - 6] = (E1 * alpha_T[i][j - 6] + E2 * alpha_T[i][j - 3] + b1) / 24;
            Q[i][j - 3] = (E2 * alpha_T[i][j - 6] + E3 * alpha_T[i][j - 3] + b2) / 24;
            Q[i][j] = (E4 * alpha_T[i][j] + b3) / 24;
        }
    }
    for (i = 0; i < 9; ++i) {
        for (j = 0; j < 9; ++j) {
            sum = 0;
            for (k = 0; k < 9; ++k) {
                sum += Q[i][k] * alpha_T[j][k];
            }
            *(pke_b_sh+i*9+j) = sum / (2 * (*(pfarea+n)));
        }
    }
}

void stiffg_sh (double *pk_sh, double *pemod, double *pnu, double *pxlocal,
    double *pthick, double *pdeffarea_ip, double *pdm, long ptr, long n)
{
    // Initialize function variables
    int i, j, k;
    double sum;
    double C[3][3]; // Constitutive matrix with plane stress coefficient
    double Bm[3][6]; // Membrane strain-displacement matrix
    double Nm[3]; // Element internal membrane forces
    double C_Bm[6][6], Bnl_N[9][6];
    double kg_sh[9][9]; // Membrane stiffness matrix

    // Define constitutive matrix with plane stress coefficient
    C[0][2] = C[1][2] = C[2][0] = C[2][1] = 0;
    C[0][0] = C[1][1] = *(pemod+ptr+n) / (1 - pow(*(pnu+n),2));
    C[0][1] = C[1][0] = *(pemod+ptr+n) / (1 - pow(*(pnu+n),2)) * (*(pnu+n));
    C[2][2] = *(pemod+ptr+n) / (1 - pow(*(pnu+n),2)) * (1 - *(pnu+n)) / 2;

    // Compute membrane strain-displacement matrix
    Bm[0][1] = Bm[0][3] = Bm[0][4] = Bm[0][5] = Bm[1][0] = Bm[1][2] = Bm[1][4] =
        Bm[2][5] = 0;
    Bm[0][0] = Bm[2][1] = -(*(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n))));
    Bm[0][2] = Bm[2][3] = *(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n)));
    Bm[1][1] = Bm[2][0] = (*(pxlocal+n*3+1) - *(pxlocal+n*3)) /
        (2 * (*(pdeffarea_ip+n)));
    Bm[1][3] = Bm[2][2] = -(*(pxlocal+n*3+1) / (2 * (*(pdeffarea_ip+n))));
    Bm[1][5] = Bm[2][4] = *(pxlocal+n*3) / (2 * (*(pdeffarea_ip+n)));

    // Compute element internal membrane forces
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 6; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += C[i][k] * Bm[k][j];
            }
            C_Bm[i][j] = sum;
        }
    }
    for (i = 0; i < 3; ++i) {
        sum = 0;
        for (j = 0; j < 6; ++j) {
            sum += C_Bm[i][j] * (*(pdm+j));
        }
        Nm[i] = *(pthick+n) * sum;
    }

    // Assemble matrix of membrane force components
    double N[6][6] =
    {
        {Nm[0], Nm[2], 0, 0, 0, 0},
        {Nm[2], Nm[1], 0, 0, 0, 0},
        {0, 0, Nm[0], Nm[2], 0, 0},
        {0, 0, Nm[2], Nm[1], 0, 0},
        {0, 0, 0, 0, Nm[0], Nm[2]},
        {0, 0, 0, 0, Nm[2], Nm[1]}
    };

    // Compute the nonlinear membrane strain-displacement matrix
    double Bnl[6][9] =
    {
        {-(*(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n)))), 0, 0,
            *(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n))), 0, 0, 0, 0, 0},
        {(*(pxlocal+n*3+1) - *(pxlocal+n*3)) / (2 * (*(pdeffarea_ip+n))), 0, 0,
            -(*(pxlocal+n*3+1) / (2 * (*(pdeffarea_ip+n)))), 0, 0,
            *(pxlocal+n*3) / (2 * (*(pdeffarea_ip+n))), 0, 0},
        {0, -(*(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n)))), 0, 0,
            *(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n))), 0, 0, 0, 0},
        {0, (*(pxlocal+n*3+1) - *(pxlocal+n*3)) / (2 * (*(pdeffarea_ip+n))), 0, 0,
            -(*(pxlocal+n*3+1) / (2 * (*(pdeffarea_ip+n)))), 0, 0,
            *(pxlocal+n*3) / (2 * (*(pdeffarea_ip+n))), 0},
        {0, 0, -(*(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n)))), 0, 0,
            *(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n))), 0, 0, 0},
        {0, 0, (*(pxlocal+n*3+1) - *(pxlocal+n*3)) / (2 * (*(pdeffarea_ip+n))), 0, 0,
            -(*(pxlocal+n*3+1) / (2 * (*(pdeffarea_ip+n)))), 0, 0,
            *(pxlocal+n*3) / (2 * (*(pdeffarea_ip+n)))}
    };

    // Compute geometrically nonlinear membrane element stiffness matrix
    for (i = 0; i < 9; ++i) {
        for (j = 0; j < 6; ++j) {
            sum = 0;
            for (k = 0; k < 6; ++k) {
                sum += Bnl[k][i] * N[k][j];
            }
            Bnl_N[i][j] = sum;
        }
    }
    for (i = 0; i < 9; ++i) {
        for (j = 0; j < 9; ++j) {
            sum = 0;
            for (k = 0; k < 6; ++k) {
                sum += Bnl_N[i][k] * Bnl[k][j];
            }
            kg_sh[i][j] = *(pdeffarea_ip+n) * sum;
        }
    }

    /* Add contribution of geometric nonlinearity from membrane to element stiffness
       matrix */
    *(pk_sh+0*18+0)   += kg_sh[0][0];
    *(pk_sh+0*18+1)   += kg_sh[0][1];
    *(pk_sh+0*18+2)   += kg_sh[0][2];
    *(pk_sh+1*18+0)   += kg_sh[1][0];
    *(pk_sh+1*18+1)   += kg_sh[1][1];
    *(pk_sh+1*18+2)   += kg_sh[1][2];
    *(pk_sh+2*18+0)   += kg_sh[2][0];
    *(pk_sh+2*18+1)   += kg_sh[2][1];
    *(pk_sh+2*18+2)   += kg_sh[2][2];
    *(pk_sh+0*18+6)   += kg_sh[0][3];
    *(pk_sh+0*18+7)   += kg_sh[0][4];
    *(pk_sh+0*18+8)   += kg_sh[0][5];
    *(pk_sh+1*18+6)   += kg_sh[1][3];
    *(pk_sh+1*18+7)   += kg_sh[1][4];
    *(pk_sh+1*18+8)   += kg_sh[1][5];
    *(pk_sh+2*18+6)   += kg_sh[2][3];
    *(pk_sh+2*18+7)   += kg_sh[2][4];
    *(pk_sh+2*18+8)   += kg_sh[2][5];
    *(pk_sh+0*18+12)  += kg_sh[0][6];
    *(pk_sh+0*18+13)  += kg_sh[0][7];
    *(pk_sh+0*18+14)  += kg_sh[0][8];
    *(pk_sh+1*18+12)  += kg_sh[1][6];
    *(pk_sh+1*18+13)  += kg_sh[1][7];
    *(pk_sh+1*18+14)  += kg_sh[1][8];
    *(pk_sh+2*18+12)  += kg_sh[2][6];
    *(pk_sh+2*18+13)  += kg_sh[2][7];
    *(pk_sh+2*18+14)  += kg_sh[2][8];
    *(pk_sh+6*18+0)   += kg_sh[3][0];
    *(pk_sh+6*18+1)   += kg_sh[3][1];
    *(pk_sh+6*18+2)   += kg_sh[3][2];
    *(pk_sh+7*18+0)   += kg_sh[4][0];
    *(pk_sh+7*18+1)   += kg_sh[4][1];
    *(pk_sh+7*18+2)   += kg_sh[4][2];
    *(pk_sh+8*18+0)   += kg_sh[5][0];
    *(pk_sh+8*18+1)   += kg_sh[5][1];
    *(pk_sh+8*18+2)   += kg_sh[5][2];
    *(pk_sh+6*18+6)   += kg_sh[3][3];
    *(pk_sh+6*18+7)   += kg_sh[3][4];
    *(pk_sh+6*18+8)   += kg_sh[3][5];
    *(pk_sh+7*18+6)   += kg_sh[4][3];
    *(pk_sh+7*18+7)   += kg_sh[4][4];
    *(pk_sh+7*18+8)   += kg_sh[4][5];
    *(pk_sh+8*18+6)   += kg_sh[5][3];
    *(pk_sh+8*18+7)   += kg_sh[5][4];
    *(pk_sh+8*18+8)   += kg_sh[5][5];
    *(pk_sh+6*18+12)  += kg_sh[3][6];
    *(pk_sh+6*18+13)  += kg_sh[3][7];
    *(pk_sh+6*18+14)  += kg_sh[3][8];
    *(pk_sh+7*18+12)  += kg_sh[4][6];
    *(pk_sh+7*18+13)  += kg_sh[4][7];
    *(pk_sh+7*18+14)  += kg_sh[4][8];
    *(pk_sh+8*18+12)  += kg_sh[5][6];
    *(pk_sh+8*18+13)  += kg_sh[5][7];
    *(pk_sh+8*18+14)  += kg_sh[5][8];
    *(pk_sh+12*18+0)  += kg_sh[6][0];
    *(pk_sh+12*18+1)  += kg_sh[6][1];
    *(pk_sh+12*18+2)  += kg_sh[6][2];
    *(pk_sh+13*18+0)  += kg_sh[7][0];
    *(pk_sh+13*18+1)  += kg_sh[7][1];
    *(pk_sh+13*18+2)  += kg_sh[7][2];
    *(pk_sh+14*18+0)  += kg_sh[8][0];
    *(pk_sh+14*18+1)  += kg_sh[8][1];
    *(pk_sh+14*18+2)  += kg_sh[8][2];
    *(pk_sh+12*18+6)  += kg_sh[6][3];
    *(pk_sh+12*18+7)  += kg_sh[6][4];
    *(pk_sh+12*18+8)  += kg_sh[6][5];
    *(pk_sh+13*18+6)  += kg_sh[7][3];
    *(pk_sh+13*18+7)  += kg_sh[7][4];
    *(pk_sh+13*18+8)  += kg_sh[7][5];
    *(pk_sh+14*18+6)  += kg_sh[8][3];
    *(pk_sh+14*18+7)  += kg_sh[8][4];
    *(pk_sh+14*18+8)  += kg_sh[8][5];
    *(pk_sh+12*18+12) += kg_sh[6][6];
    *(pk_sh+12*18+13) += kg_sh[6][7];
    *(pk_sh+12*18+14) += kg_sh[6][8];
    *(pk_sh+13*18+12) += kg_sh[7][6];
    *(pk_sh+13*18+13) += kg_sh[7][7];
    *(pk_sh+13*18+14) += kg_sh[7][8];
    *(pk_sh+14*18+12) += kg_sh[8][6];
    *(pk_sh+14*18+13) += kg_sh[8][7];
    *(pk_sh+14*18+14) += kg_sh[8][8];
}

void stiffm_sh (double *pk_sh, double *pemod, double *pnu, double *pxlocal,
    double *pthick, double *pdeffarea_ip, double *pdefslen_ip, double *pyield,
    double *pefN_temp, double *pefM_temp, double *pchi_temp, int yv, double *pNo,
    double *palpha, double *pMe, double *pNbar, double *pMbar, double *pMNbar,
    double *pq_fact, double *pr_fact, double *ps_fact, int *ph_fact, long ptr,
    long n)
{
    // Initialize function variables
    int j, k;
    double sum, sum2;
    // Membrane, plate bending, and coupled stiffness matrices
    double km_m_sh[6][6], km_b_sh[9][9], km_mb_sh[6][9];
    double C[3][3]; // Plane stress constitutive matrix
    double c_fact, g_fact, d_fact; // Factors for computation of plastic flow directions
    double gradN_Nbar[3], gradM_Mbar[3]; // Gradients of quadratic stress intensities
    double fn[3], fm[3]; // Plastic flow directions
    // Factors for computation of elasto-plastic modular matrices
    double fn_C[3], fm_C[3];
    double j_fact, k_fact, B_fact;
    double df_da, da_dchi; // Derivatives involved with pseudo hardening parameter

    // Compute plane stress constitutive matrix
    C[0][2] = C[1][2] = C[2][0] = C[2][1] = 0;
    C[0][0] = C[1][1] = *(pemod+ptr+n) / (1 - pow(*(pnu+n),2));
    C[0][1] = C[1][0] = *(pemod+ptr+n) / (1 - pow(*(pnu+n),2)) * (*(pnu+n));
    C[2][2] = *(pemod+ptr+n) / (1 - pow(*(pnu+n),2)) * (1 - *(pnu+n)) / 2;

    // Compute factors for computation of plastic flow directions
    c_fact = 1 / pow(*pNo,2) - *pMbar / (4 * (*pq_fact)) +
        *ps_fact * pow(*pMe,2) / (4 * pow(*pq_fact,2));
    if (*ph_fact == 1) {
        g_fact = *pMNbar * (1 / (4 * (*pq_fact)) + 1 / (*pNo * (*pr_fact)));
        d_fact = 1  / (2 * pow(*pMe,2)) - *pNbar / (4 * (*pq_fact)) +
            0.12 * pow(*pNo,2) * (*ps_fact) / pow(*pq_fact,2) +
            *pMbar * (*pNo) / (2 * pow(*pMe,2) * (*pr_fact));
    } else {
        g_fact = *pMNbar / (4 * (*pq_fact));
        d_fact = 1  / (2 * pow(*pMe,2)) - *pNbar / (4 * (*pq_fact)) +
            0.12 * pow(*pNo,2) * (*ps_fact) / pow(*pq_fact,2);
    }

    // Compute gradients of quadratic stress intensities
    gradN_Nbar[0] = 2 * (*(pefN_temp+n*9+yv*3)) - *(pefN_temp+n*9+yv*3+1);
    gradN_Nbar[1] = 2 * (*(pefN_temp+n*9+yv*3+1)) - *(pefN_temp+n*9+yv*3);
    gradN_Nbar[2] = 6 * (*(pefN_temp+n*9+yv*3+2));
    gradM_Mbar[0] = 2 * (*(pefM_temp+n*9+yv*3)) - *(pefM_temp+n*9+yv*3+1);
    gradM_Mbar[1] = 2 * (*(pefM_temp+n*9+yv*3+1)) - *(pefM_temp+n*9+yv*3);
    gradM_Mbar[2] = 6 * (*(pefM_temp+n*9+yv*3+2));

    // Compute plastic flow directions
    for (j = 0; j < 3; ++j) {
        fn[j] = c_fact * gradN_Nbar[j] + g_fact * gradM_Mbar[j];
        fm[j] = g_fact * gradN_Nbar[j] + d_fact * gradM_Mbar[j];
    }

    /* Compute factors for computation of elasto-plastic modular matrices in stiffm_m_sh,
       stiffm_b_sh, and stiffm_mb_sh functions */
    for (j = 0; j < 3; ++j) {
        sum = sum2 = 0;
        for (k = 0; k < 3; ++k) {
            sum += fn[k] * C[k][j];
            sum2 += fm[k] * C[k][j];
        }
        fn_C[j] = sum;
        fm_C[j] = sum2;
    }
    sum = sum2 = 0;
    for (j = 0; j < 3; ++j) {
        sum += fn_C[j] * fn[j];
        sum2 += fm_C[j] * fm[j];
    }
    j_fact = *(pthick+n) * sum;
    k_fact = pow(*(pthick+n),3) * sum2 / 12;
    B_fact = 2 * sqrt(pow(g_fact,2) * (*pNbar) + pow(d_fact,2) * (*pMbar) +
        2 * d_fact * g_fact * (*pMNbar));

    // Compute derivatives involved with pseudo hardening parameter
    if (*ph_fact == 1) {
        df_da = -(*pMbar / (*palpha * pow(*pMe,2))) +
            *ps_fact * (*pNbar) * pow(*pMe,2) / (2 * pow(*pq_fact,2) * (*palpha)) -
            *pr_fact / (*palpha * pow(*pMe,2) * (*pNo)) + 2 * pow(*pMNbar,2) /
            (*palpha * (*pNo) * (*pr_fact));
    } else {
        df_da = -(*pMbar / (*palpha * pow(*pMe,2))) +
            *ps_fact * (*pNbar) * pow(*pMe,2) / (2 * pow(*pq_fact,2) * (*palpha));
    }
    if (*(pchi_temp+n*3+yv) >= 1e-6) {
        da_dchi = 0.52 * (*(pemod+ptr+n)) * (*(pthick+n)) *
            exp(-2.6 * sqrt(*(pchi_temp+n*3+yv))) /
            (3 * (*(pyield+ptr+n)) * sqrt(*(pchi_temp+n*3+yv)));
    } else {
        da_dchi = 0;
    }

    // Pass control to stiffm_m_sh function
    stiffm_m_sh (&km_m_sh[0][0], pxlocal, pthick, pdeffarea_ip, &C[0][0], fn, &j_fact,
        &k_fact, &B_fact, &df_da, &da_dchi, n);

    /* Add contribution of elasto-plastic membrane stiffness matrix to element stiffness
       matrix */
    *(pk_sh+0*18+0)   = km_m_sh[0][0];
    *(pk_sh+0*18+1)   = km_m_sh[0][1];
    *(pk_sh+1*18+0)   = km_m_sh[1][0];
    *(pk_sh+1*18+1)   = km_m_sh[1][1];
    *(pk_sh+0*18+6)   = km_m_sh[0][2];
    *(pk_sh+0*18+7)   = km_m_sh[0][3];
    *(pk_sh+1*18+6)   = km_m_sh[1][2];
    *(pk_sh+1*18+7)   = km_m_sh[1][3];
    *(pk_sh+0*18+12)  = km_m_sh[0][4];
    *(pk_sh+0*18+13)  = km_m_sh[0][5];
    *(pk_sh+1*18+12)  = km_m_sh[1][4];
    *(pk_sh+1*18+13)  = km_m_sh[1][5];
    *(pk_sh+6*18+0)   = km_m_sh[2][0];
    *(pk_sh+6*18+1)   = km_m_sh[2][1];
    *(pk_sh+7*18+0)   = km_m_sh[3][0];
    *(pk_sh+7*18+1)   = km_m_sh[3][1];
    *(pk_sh+6*18+6)   = km_m_sh[2][2];
    *(pk_sh+6*18+7)   = km_m_sh[2][3];
    *(pk_sh+7*18+6)   = km_m_sh[3][2];
    *(pk_sh+7*18+7)   = km_m_sh[3][3];
    *(pk_sh+6*18+12)  = km_m_sh[2][4];
    *(pk_sh+6*18+13)  = km_m_sh[2][5];
    *(pk_sh+7*18+12)  = km_m_sh[3][4];
    *(pk_sh+7*18+13)  = km_m_sh[3][5];
    *(pk_sh+12*18+0)  = km_m_sh[4][0];
    *(pk_sh+12*18+1)  = km_m_sh[4][1];
    *(pk_sh+13*18+0)  = km_m_sh[5][0];
    *(pk_sh+13*18+1)  = km_m_sh[5][1];
    *(pk_sh+12*18+6)  = km_m_sh[4][2];
    *(pk_sh+12*18+7)  = km_m_sh[4][3];
    *(pk_sh+13*18+6)  = km_m_sh[5][2];
    *(pk_sh+13*18+7)  = km_m_sh[5][3];
    *(pk_sh+12*18+12) = km_m_sh[4][4];
    *(pk_sh+12*18+13) = km_m_sh[4][5];
    *(pk_sh+13*18+12) = km_m_sh[5][4];
    *(pk_sh+13*18+13) = km_m_sh[5][5];

    // Pass control to stiffm_b_sh function
    stiffm_b_sh (&km_b_sh[0][0], pxlocal, pthick, pdeffarea_ip, pdefslen_ip, &C[0][0],
        fm, &j_fact, &k_fact, &B_fact, &df_da, &da_dchi, n);

    /* Add contribution of elasto-plastic plate bending stiffness matrix to element
       stiffness matrix */
    *(pk_sh+2*18+2)   = km_b_sh[0][0];
    *(pk_sh+2*18+3)   = km_b_sh[0][1];
    *(pk_sh+2*18+4)   = km_b_sh[0][2];
    *(pk_sh+3*18+2)   = km_b_sh[1][0];
    *(pk_sh+3*18+3)   = km_b_sh[1][1];
    *(pk_sh+3*18+4)   = km_b_sh[1][2];
    *(pk_sh+4*18+2)   = km_b_sh[2][0];
    *(pk_sh+4*18+3)   = km_b_sh[2][1];
    *(pk_sh+4*18+4)   = km_b_sh[2][2];
    *(pk_sh+2*18+8)   = km_b_sh[0][3];
    *(pk_sh+2*18+9)   = km_b_sh[0][4];
    *(pk_sh+2*18+10)  = km_b_sh[0][5];
    *(pk_sh+3*18+8)   = km_b_sh[1][3];
    *(pk_sh+3*18+9)   = km_b_sh[1][4];
    *(pk_sh+3*18+10)  = km_b_sh[1][5];
    *(pk_sh+4*18+8)   = km_b_sh[2][3];
    *(pk_sh+4*18+9)   = km_b_sh[2][4];
    *(pk_sh+4*18+10)  = km_b_sh[2][5];
    *(pk_sh+2*18+14)  = km_b_sh[0][6];
    *(pk_sh+2*18+15)  = km_b_sh[0][7];
    *(pk_sh+2*18+16)  = km_b_sh[0][8];
    *(pk_sh+3*18+14)  = km_b_sh[1][6];
    *(pk_sh+3*18+15)  = km_b_sh[1][7];
    *(pk_sh+3*18+16)  = km_b_sh[1][8];
    *(pk_sh+4*18+14)  = km_b_sh[2][6];
    *(pk_sh+4*18+15)  = km_b_sh[2][7];
    *(pk_sh+4*18+16)  = km_b_sh[2][8];
    *(pk_sh+8*18+2)   = km_b_sh[3][0];
    *(pk_sh+8*18+3)   = km_b_sh[3][1];
    *(pk_sh+8*18+4)   = km_b_sh[3][2];
    *(pk_sh+9*18+2)   = km_b_sh[4][0];
    *(pk_sh+9*18+3)   = km_b_sh[4][1];
    *(pk_sh+9*18+4)   = km_b_sh[4][2];
    *(pk_sh+10*18+2)  = km_b_sh[5][0];
    *(pk_sh+10*18+3)  = km_b_sh[5][1];
    *(pk_sh+10*18+4)  = km_b_sh[5][2];
    *(pk_sh+8*18+8)   = km_b_sh[3][3];
    *(pk_sh+8*18+9)   = km_b_sh[3][4];
    *(pk_sh+8*18+10)  = km_b_sh[3][5];
    *(pk_sh+9*18+8)   = km_b_sh[4][3];
    *(pk_sh+9*18+9)   = km_b_sh[4][4];
    *(pk_sh+9*18+10)  = km_b_sh[4][5];
    *(pk_sh+10*18+8)  = km_b_sh[5][3];
    *(pk_sh+10*18+9)  = km_b_sh[5][4];
    *(pk_sh+10*18+10) = km_b_sh[5][5];
    *(pk_sh+8*18+14)  = km_b_sh[3][6];
    *(pk_sh+8*18+15)  = km_b_sh[3][7];
    *(pk_sh+8*18+16)  = km_b_sh[3][8];
    *(pk_sh+9*18+14)  = km_b_sh[4][6];
    *(pk_sh+9*18+15)  = km_b_sh[4][7];
    *(pk_sh+9*18+16)  = km_b_sh[4][8];
    *(pk_sh+10*18+14) = km_b_sh[5][6];
    *(pk_sh+10*18+15) = km_b_sh[5][7];
    *(pk_sh+10*18+16) = km_b_sh[5][8];
    *(pk_sh+14*18+2)  = km_b_sh[6][0];
    *(pk_sh+14*18+3)  = km_b_sh[6][1];
    *(pk_sh+14*18+4)  = km_b_sh[6][2];
    *(pk_sh+15*18+2)  = km_b_sh[7][0];
    *(pk_sh+15*18+3)  = km_b_sh[7][1];
    *(pk_sh+15*18+4)  = km_b_sh[7][2];
    *(pk_sh+16*18+2)  = km_b_sh[8][0];
    *(pk_sh+16*18+3)  = km_b_sh[8][1];
    *(pk_sh+16*18+4)  = km_b_sh[8][2];
    *(pk_sh+14*18+8)  = km_b_sh[6][3];
    *(pk_sh+14*18+9)  = km_b_sh[6][4];
    *(pk_sh+14*18+10) = km_b_sh[6][5];
    *(pk_sh+15*18+8)  = km_b_sh[7][3];
    *(pk_sh+15*18+9)  = km_b_sh[7][4];
    *(pk_sh+15*18+10) = km_b_sh[7][5];
    *(pk_sh+16*18+8)  = km_b_sh[8][3];
    *(pk_sh+16*18+9)  = km_b_sh[8][4];
    *(pk_sh+16*18+10) = km_b_sh[8][5];
    *(pk_sh+14*18+14) = km_b_sh[6][6];
    *(pk_sh+14*18+15) = km_b_sh[6][7];
    *(pk_sh+14*18+16) = km_b_sh[6][8];
    *(pk_sh+15*18+14) = km_b_sh[7][6];
    *(pk_sh+15*18+15) = km_b_sh[7][7];
    *(pk_sh+15*18+16) = km_b_sh[7][8];
    *(pk_sh+16*18+14) = km_b_sh[8][6];
    *(pk_sh+16*18+15) = km_b_sh[8][7];
    *(pk_sh+16*18+16) = km_b_sh[8][8];

    // Pass control to stiffm_b function
    stiffm_mb_sh (&km_mb_sh[0][0], pxlocal, pthick, pdeffarea_ip, pdefslen_ip, &C[0][0],
        fn, fm, &j_fact, &k_fact, &B_fact, &df_da, &da_dchi, n);

    /* Add contribution of elasto-plastic coupled stiffness matrix to element stiffness
       matrix */
    *(pk_sh+0*18+2)   = *(pk_sh+2*18+0)   = km_mb_sh[0][0];
    *(pk_sh+0*18+3)   = *(pk_sh+3*18+0)   = km_mb_sh[0][1];
    *(pk_sh+0*18+4)   = *(pk_sh+4*18+0)   = km_mb_sh[0][2];
    *(pk_sh+0*18+8)   = *(pk_sh+8*18+0)   = km_mb_sh[0][3];
    *(pk_sh+0*18+9)   = *(pk_sh+9*18+0)   = km_mb_sh[0][4];
    *(pk_sh+0*18+10)  = *(pk_sh+10*18+0)  = km_mb_sh[0][5];
    *(pk_sh+0*18+14)  = *(pk_sh+14*18+0)  = km_mb_sh[0][6];
    *(pk_sh+0*18+15)  = *(pk_sh+15*18+0)  = km_mb_sh[0][7];
    *(pk_sh+0*18+16)  = *(pk_sh+16*18+0)  = km_mb_sh[0][8];
    *(pk_sh+1*18+2)   = *(pk_sh+2*18+1)   = km_mb_sh[1][0];
    *(pk_sh+1*18+3)   = *(pk_sh+3*18+1)   = km_mb_sh[1][1];
    *(pk_sh+1*18+4)   = *(pk_sh+4*18+1)   = km_mb_sh[1][2];
    *(pk_sh+1*18+8)   = *(pk_sh+8*18+1)   = km_mb_sh[1][3];
    *(pk_sh+1*18+9)   = *(pk_sh+9*18+1)   = km_mb_sh[1][4];
    *(pk_sh+1*18+10)  = *(pk_sh+10*18+1)  = km_mb_sh[1][5];
    *(pk_sh+1*18+14)  = *(pk_sh+14*18+1)  = km_mb_sh[1][6];
    *(pk_sh+1*18+15)  = *(pk_sh+15*18+1)  = km_mb_sh[1][7];
    *(pk_sh+1*18+16)  = *(pk_sh+16*18+1)  = km_mb_sh[1][8];
    *(pk_sh+6*18+2)   = *(pk_sh+2*18+6)   = km_mb_sh[2][0];
    *(pk_sh+6*18+3)   = *(pk_sh+3*18+6)   = km_mb_sh[2][1];
    *(pk_sh+6*18+4)   = *(pk_sh+4*18+6)   = km_mb_sh[2][2];
    *(pk_sh+6*18+8)   = *(pk_sh+8*18+6)   = km_mb_sh[2][3];
    *(pk_sh+6*18+9)   = *(pk_sh+9*18+6)   = km_mb_sh[2][4];
    *(pk_sh+6*18+10)  = *(pk_sh+10*18+6)  = km_mb_sh[2][5];
    *(pk_sh+6*18+14)  = *(pk_sh+14*18+6)  = km_mb_sh[2][6];
    *(pk_sh+6*18+15)  = *(pk_sh+15*18+6)  = km_mb_sh[2][7];
    *(pk_sh+6*18+16)  = *(pk_sh+16*18+6)  = km_mb_sh[2][8];
    *(pk_sh+7*18+2)   = *(pk_sh+2*18+7)   = km_mb_sh[3][0];
    *(pk_sh+7*18+3)   = *(pk_sh+3*18+7)   = km_mb_sh[3][1];
    *(pk_sh+7*18+4)   = *(pk_sh+4*18+7)   = km_mb_sh[3][2];
    *(pk_sh+7*18+8)   = *(pk_sh+8*18+7)   = km_mb_sh[3][3];
    *(pk_sh+7*18+9)   = *(pk_sh+9*18+7)   = km_mb_sh[3][4];
    *(pk_sh+7*18+10)  = *(pk_sh+10*18+7)  = km_mb_sh[3][5];
    *(pk_sh+7*18+14)  = *(pk_sh+14*18+7)  = km_mb_sh[3][6];
    *(pk_sh+7*18+15)  = *(pk_sh+15*18+7)  = km_mb_sh[3][7];
    *(pk_sh+7*18+16)  = *(pk_sh+16*18+7)  = km_mb_sh[3][8];
    *(pk_sh+12*18+2)  = *(pk_sh+2*18+12)  = km_mb_sh[4][0];
    *(pk_sh+12*18+3)  = *(pk_sh+3*18+12)  = km_mb_sh[4][1];
    *(pk_sh+12*18+4)  = *(pk_sh+4*18+12)  = km_mb_sh[4][2];
    *(pk_sh+12*18+8)  = *(pk_sh+8*18+12)  = km_mb_sh[4][3];
    *(pk_sh+12*18+9)  = *(pk_sh+9*18+12)  = km_mb_sh[4][4];
    *(pk_sh+12*18+10) = *(pk_sh+10*18+12) = km_mb_sh[4][5];
    *(pk_sh+12*18+14) = *(pk_sh+14*18+12) = km_mb_sh[4][6];
    *(pk_sh+12*18+15) = *(pk_sh+15*18+12) = km_mb_sh[4][7];
    *(pk_sh+12*18+16) = *(pk_sh+16*18+12) = km_mb_sh[4][8];
    *(pk_sh+13*18+2)  = *(pk_sh+2*18+13)  = km_mb_sh[5][0];
    *(pk_sh+13*18+3)  = *(pk_sh+3*18+13)  = km_mb_sh[5][1];
    *(pk_sh+13*18+4)  = *(pk_sh+4*18+13)  = km_mb_sh[5][2];
    *(pk_sh+13*18+8)  = *(pk_sh+8*18+13)  = km_mb_sh[5][3];
    *(pk_sh+13*18+9)  = *(pk_sh+9*18+13)  = km_mb_sh[5][4];
    *(pk_sh+13*18+10) = *(pk_sh+10*18+13) = km_mb_sh[5][5];
    *(pk_sh+13*18+14) = *(pk_sh+14*18+13) = km_mb_sh[5][6];
    *(pk_sh+13*18+15) = *(pk_sh+15*18+13) = km_mb_sh[5][7];
    *(pk_sh+13*18+16) = *(pk_sh+16*18+13) = km_mb_sh[5][8];

    // Assign in-plane rotation stiffness to element stiffness matrix
    *(pk_sh+5*18+5)   = km_b_sh[1][1] / 10000;
    *(pk_sh+11*18+11) = km_b_sh[4][4] / 10000;
    *(pk_sh+17*18+17) = km_b_sh[7][7] / 10000;
}

void stiffm_m_sh (double *pkm_m_sh, double *pxlocal, double *pthick,
    double *pdeffarea_ip, double *pC, double *pfn, double *pj_fact, double *pk_fact,
    double *pB_fact, double *pdf_da, double *pda_dchi, long n)
{
    // Initialize function variables
    int i, j, k;
    double sum;
    // Factor for computation of membrane elasto-plastic modular matrix
    double xi = *(pthick+n) / (*pj_fact + *pk_fact - *pB_fact * (*pdf_da) * (*pda_dchi));
    double C_star[3][3]; // Membrane elasto-plastic modular matrix
    double Bm[3][6]; // Membrane strain-displacement matrix
    double N[3][3], xi_N_C[3][3], Bm_Cstar[6][3];

    // Compute membrane elasto-plastic modular matrix
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            N[i][j] = *(pfn+i) * (*(pfn+j));
        }
        for (j = 0; j < 3; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += N[i][k] * (*(pC+k*3+j));
            }
            xi_N_C[i][j] = xi * sum;
        }
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += (*(pC+i*3+k)) * xi_N_C[k][j];
            }
            C_star[i][j] = *(pthick+n) * (*(pC+i*3+j) - sum);
        }
    }

    // Compute membrane strain-displacement matrix
    Bm[0][1] = Bm[0][3] = Bm[0][4] = Bm[0][5] = Bm[1][0] = Bm[1][2] = Bm[1][4] =
        Bm[2][5] = 0;
    Bm[0][0] = Bm[2][1] = -(*(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n))));
    Bm[0][2] = Bm[2][3] = *(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n)));
    Bm[1][1] = Bm[2][0] = (*(pxlocal+n*3+1) - *(pxlocal+n*3)) /
        (2 * (*(pdeffarea_ip+n)));
    Bm[1][3] = Bm[2][2] = -(*(pxlocal+n*3+1) / (2 * (*(pdeffarea_ip+n))));
    Bm[1][5] = Bm[2][4] = *(pxlocal+n*3) / (2 * (*(pdeffarea_ip+n)));

    // Compute elasto-plastic membrane element stiffness matrix
    for (i = 0; i < 6; ++i) {
        for (j = 0; j < 3; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += Bm[k][i] * C_star[k][j];
            }
            Bm_Cstar[i][j] = sum;
        }
    }
    for (i = 0; i < 6; ++i) {
        for (j = 0; j < 6; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += Bm_Cstar[i][k] * Bm[k][j];
            }
            *(pkm_m_sh+i*6+j) = *(pdeffarea_ip+n) * sum;
        }
    }
}

void stiffm_b_sh (double *pkm_b_sh, double *pxlocal, double *pthick,
    double *pdeffarea_ip, double *pdefslen_ip, double *pC, double *pfm, double *pj_fact,
    double *pk_fact, double *pB_fact, double *pdf_da, double *pda_dchi, long n)
{
    // Initialize function variables
    int i, j, k;
    double sum;
    // Factor for computation of plate bending elasto-plastic modular matrix
    double xi = pow(*(pthick+n),3) /
        (12 * (*pj_fact + *pk_fact - *pB_fact * (*pdf_da) * (*pda_dchi)));
    double D_star[3][3]; // Plate bending elasto-plastic modular matrix
    // Coefficients for computing plate bending strain-displacement matrix
    double x23, l12, l23, l31, p4, p5, p6, t4, t5, q4, q5, r4, r5;
    double Q[9][9], b[3];
    double M[3][3], xi_M_C[3][3];

    // Compute plate bending elasto-plastic modular matrix
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            M[i][j] = *(pfm+i) * (*(pfm+j));
        }
        for (j = 0; j < 3; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += M[i][k] * (*(pC+k*3+j));
            }
            xi_M_C[i][j] = xi * sum;
        }
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += (*(pC+i*3+k)) * xi_M_C[k][j];
            }
            D_star[i][j] = pow(*(pthick+n),3) * (*(pC+i*3+j) - sum) / 12;
        }
    }

    // Define local coordinate coefficients for convenience in defining p, t, q, and r
    x23 = *(pxlocal+n*3) - *(pxlocal+n*3+1);
    l12 = pow(*(pdefslen_ip+n*3),2);
    l23 = pow(*(pdefslen_ip+n*3+1),2);
    l31 = pow(*(pdefslen_ip+n*3+2),2);

    // Compute coefficients for assembly of plate bending strain-displacement matrix
    p4 = -6 * x23 / l23;
    p5 = -6 * (*(pxlocal+n*3+1)) / l31;
    p6 = 6 * (*(pxlocal+n*3)) / l12;
    t4 = 6 * (*(pxlocal+n*3+2)) / l23;
    t5 = -6 * (*(pxlocal+n*3+2)) / l31;
    q4 = -3 * x23 * (*(pxlocal+n*3+2)) / l23;
    q5 = 3 * (*(pxlocal+n*3+1)) * (*(pxlocal+n*3+2)) / l31;
    r4 = 3 * pow(*(pxlocal+n*3+2),2) / l23;
    r5 = 3 * pow(*(pxlocal+n*3+2),2) / l31;

    /* Compute the TRANSPOSE of the alpha matrix; computing the transpose eliminates a
       nested loop operation */
    double alpha_T[9][9] =
    {
        {*(pxlocal+n*3+2) * p6, -(*(pxlocal+n*3+2) * p6), *(pxlocal+n*3+2) * p5,
            -(*(pxlocal+n*3) * t5), 0, x23 * t5,
            -(*(pxlocal+n*3+1) * p6) - *(pxlocal+n*3) * p5, -x23 * p6,
            x23 * p5 + *(pxlocal+n*3+2) * t5},
        {0, 0, -(*(pxlocal+n*3+2) * q5), x23 + *(pxlocal+n*3) * r5, x23, x23 * (1 - r5),
            *(pxlocal+n*3) * q5 + *(pxlocal+n*3+2), *(pxlocal+n*3+2),
            -x23 * q5 + *(pxlocal+n*3+2) * (1 - r5)},
        {-4 * (*(pxlocal+n*3+2)), 2 * (*(pxlocal+n*3+2)), *(pxlocal+n*3+2) * (2 - r5),
            -(*(pxlocal+n*3) * q5), 0, x23 * q5, -4 * x23 + *(pxlocal+n*3) * r5, 2 * x23,
            x23 * (2 - r5) + *(pxlocal+n*3+2) * q5},
        {-(*(pxlocal+n*3+2) * p6), *(pxlocal+n*3+2) * p6, *(pxlocal+n*3+2) * p4, 0,
            *(pxlocal+n*3) * t4, -(*(pxlocal+n*3+1) * t4), *(pxlocal+n*3+1) * p6,
            x23 * p6 + *(pxlocal+n*3) * p4,
            -(*(pxlocal+n*3+1) * p4) + *(pxlocal+n*3+2) * t4},
        {0, 0, *(pxlocal+n*3+2) * q4, *(pxlocal+n*3+1),
            *(pxlocal+n*3+1) + *(pxlocal+n*3) * r4, *(pxlocal+n*3+1) * (1 - r4),
            -(*(pxlocal+n*3+2)), -(*(pxlocal+n*3+2)) + *(pxlocal+n*3) * q4,
            *(pxlocal+n*3+2) * (r4 - 1) - *(pxlocal+n*3+1) * q4},
        {-2 * (*(pxlocal+n*3+2)), 4 * (*(pxlocal+n*3+2)), *(pxlocal+n*3+2) * (r4 - 2), 0,
            -(*(pxlocal+n*3) * q4), *(pxlocal+n*3+1) * q4, 2 * (*(pxlocal+n*3+1)),
            -4 * (*(pxlocal+n*3+1)) + *(pxlocal+n*3) * r4,
            *(pxlocal+n*3+1) * (2 - r4) - *(pxlocal+n*3+2) * q4},
        {0, 0, -(*(pxlocal+n*3+2) * (p4 + p5)), *(pxlocal+n*3) * t5,
            -(*(pxlocal+n*3) * t4), -x23 * t5 + *(pxlocal+n*3+1) * t4,
            *(pxlocal+n*3) * p5, -(*(pxlocal+n*3) * p4),
            -x23 * p5 + *(pxlocal+n*3+1) * p4 - *(pxlocal+n*3+2) * (t4 + t5)},
        {0, 0, *(pxlocal+n*3+2) * (q4 - q5), *(pxlocal+n*3) * (r5 - 1),
            *(pxlocal+n*3) * (r4 - 1),
            -x23 * r5 - *(pxlocal+n*3+1) * r4 - *(pxlocal+n*3), *(pxlocal+n*3) * q5,
            *(pxlocal+n*3) * q4,
            -x23 * q5 - *(pxlocal+n*3+1) * q4 + *(pxlocal+n*3+2) * (r4 - r5)},
        {0, 0, *(pxlocal+n*3+2) * (r4 - r5), -(*(pxlocal+n*3) * q5),
            -(*(pxlocal+n*3) * q4), *(pxlocal+n*3+1) * q4 + x23 * q5,
            *(pxlocal+n*3) * (r5 - 2), *(pxlocal+n*3) * (r4 - 2),
            -x23 * r5 - *(pxlocal+n*3+1) * r4 + 4 * (*(pxlocal+n*3)) +
            *(pxlocal+n*3+2) * (q5 - q4)}
    };

    // Compute elasto-plastic plate bending element stiffness matrix
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            b[j] = 0;
            for (k = 0; k < 3; ++k) {
                b[j] += D_star[j][0] * alpha_T[i][k] + D_star[j][1] * alpha_T[i][k + 3] +
                    D_star[j][2] * alpha_T[i][k + 6];
            }
        }
        for (j = 0; j < 3; ++j) {
            Q[i][j] = (D_star[0][0] * alpha_T[i][j] + D_star[1][0] * alpha_T[i][j + 3] +
                D_star[2][0] * alpha_T[i][j + 6] + b[0]) / 24;
            Q[i][j + 3] = (D_star[0][1] * alpha_T[i][j] +
                D_star[1][1] * alpha_T[i][j + 3] + D_star[2][1] * alpha_T[i][j + 6] +
                b[1]) / 24;
            Q[i][j + 6] = (D_star[0][2] * alpha_T[i][j] +
                D_star[1][2] * alpha_T[i][j + 3] + D_star[2][2] * alpha_T[i][j + 6] +
                b[2]) / 24;
        }
    }
    for (i = 3; i < 6; ++i) {
        for (j = 0; j < 3; ++j) {
            b[j] = 0;
            for (k = 3; k < 6; ++k) {
                b[j] += D_star[j][0] * alpha_T[i][k - 3] + D_star[j][1] * alpha_T[i][k] +
                    D_star[j][2] * alpha_T[i][k + 3];
            }
        }
        for (j = 3; j < 6; ++j) {
            Q[i][j - 3] = (D_star[0][0] * alpha_T[i][j - 3] +
                D_star[1][0] * alpha_T[i][j] + D_star[2][0] * alpha_T[i][j + 3] + b[0])
                / 24;
            Q[i][j] = (D_star[0][1] * alpha_T[i][j - 3] +
                D_star[1][1] * alpha_T[i][j] + D_star[2][1] * alpha_T[i][j + 3] + b[1])
                / 24;
            Q[i][j + 3] = (D_star[0][2] * alpha_T[i][j - 3] +
                D_star[1][2] * alpha_T[i][j] + D_star[2][2] * alpha_T[i][j + 3] + b[2])
                / 24;
        }
    }
    for (i = 6; i < 9; ++i) {
        for (j = 0; j < 3; ++j) {
            b[j] = 0;
            for (k = 6; k < 9; ++k) {
                b[j] += D_star[j][0] * alpha_T[i][k - 6] +
                     D_star[j][1] * alpha_T[i][k - 3] + D_star[j][2] * alpha_T[i][k];
            }
        }
        for (j = 6; j < 9; ++j) {
            Q[i][j - 6] = (D_star[0][0] * alpha_T[i][j - 6] +
                D_star[1][0] * alpha_T[i][j - 3] + D_star[2][0] * alpha_T[i][j] + b[0])
                / 24;
            Q[i][j - 3] = (D_star[0][1] * alpha_T[i][j - 6] +
                D_star[1][1] * alpha_T[i][j - 3] + D_star[2][1] * alpha_T[i][j] + b[1])
                / 24;
            Q[i][j] = (D_star[0][2] * alpha_T[i][j - 6] +
                D_star[1][2] * alpha_T[i][j - 3] + D_star[2][2] * alpha_T[i][j] + b[2])
                / 24;
        }
    }
    for (i = 0; i < 9; ++i) {
        for (j = 0; j < 9; ++j) {
            sum = 0;
            for (k = 0; k < 9; ++k) {
                sum += Q[i][k] * alpha_T[j][k];
            }
            *(pkm_b_sh+i*9+j) = sum / (2 * (*(pdeffarea_ip+n)));
        }
    }
}

void stiffm_mb_sh (double *pkm_mb_sh, double *pxlocal, double *pthick,
    double *pdeffarea_ip, double *pdefslen_ip, double *pC, double *pfn, double *pfm,
    double *pj_fact, double *pk_fact, double *pB_fact, double *pdf_da, double *pda_dchi,
    long n)
{
    // Initialize function variables
    int i, j, k;
    double sum;
    // Factor for computation of coupled elasto-plastic modular matrix
    double xi = -pow(*(pthick+n),4) /
        (12 * (*pj_fact + *pk_fact - *pB_fact * (*pdf_da) * (*pda_dchi)));
    double cd[3][3]; // Coupled elasto-plastic modular matrix
    double Bm[3][6]; // Membrane strain-displacement matrix
    // Coefficients for computing plate bending strain-displacement matrix
    double x23, l12, l23, l31, p4, p5, p6, t4, t5, q4, q5, r4, r5;
    double NM[3][3], xi_NM_C[3][3], cdL[3][9], Bm_cdL[6][9];

    // Compute coupled elasto-plastic modular matrix
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            NM[i][j] = *(pfn+i) * (*(pfm+j));
        }
        for (j = 0; j < 3; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += NM[i][k] * (*(pC+k*3+j));
            }
            xi_NM_C[i][j] = xi * sum;
        }
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += (*(pC+i*3+k)) * xi_NM_C[k][j];
            }
            cd[i][j] = sum;
        }
    }

    // Compute membrane strain-displacement matrix
    Bm[0][1] = Bm[0][3] = Bm[0][4] = Bm[0][5] = Bm[1][0] = Bm[1][2] = Bm[1][4] =
        Bm[2][5] = 0;
    Bm[0][0] = Bm[2][1] = -(*(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n))));
    Bm[0][2] = Bm[2][3] = *(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n)));
    Bm[1][1] = Bm[2][0] = (*(pxlocal+n*3+1) - *(pxlocal+n*3)) /
        (2 * (*(pdeffarea_ip+n)));
    Bm[1][3] = Bm[2][2] = -(*(pxlocal+n*3+1) / (2 * (*(pdeffarea_ip+n))));
    Bm[1][5] = Bm[2][4] = *(pxlocal+n*3) / (2 * (*(pdeffarea_ip+n)));

    // Define local coordinate coefficients for convenience in defining p, t, q, and r
    x23 = *(pxlocal+n*3) - *(pxlocal+n*3+1);
    l12 = pow(*(pdefslen_ip+n*3),2);
    l23 = pow(*(pdefslen_ip+n*3+1),2);
    l31 = pow(*(pdefslen_ip+n*3+2),2);

    // Compute coefficients for assembly of plate bending strain-displacement matrix
    p4 = -6 * x23 / l23;
    p5 = -6 * (*(pxlocal+n*3+1)) / l31;
    p6 = 6 * (*(pxlocal+n*3)) / l12;
    t4 = 6 * (*(pxlocal+n*3+2)) / l23;
    t5 = -6 * (*(pxlocal+n*3+2)) / l31;
    q4 = -3 * x23 * (*(pxlocal+n*3+2)) / l23;
    q5 = 3 * (*(pxlocal+n*3+1)) * (*(pxlocal+n*3+2)) / l31;
    r4 = 3 * pow(*(pxlocal+n*3+2),2) / l23;
    r5 = 3 * pow(*(pxlocal+n*3+2),2) / l31;

    /* Compute the TRANSPOSE of the alpha matrix; computing the transpose eliminates a
       nested loop operation */
    double alpha_T[9][9] =
    {
        {*(pxlocal+n*3+2) * p6, -(*(pxlocal+n*3+2) * p6), *(pxlocal+n*3+2) * p5,
            -(*(pxlocal+n*3) * t5), 0, x23 * t5,
            -(*(pxlocal+n*3+1) * p6) - *(pxlocal+n*3) * p5, -x23 * p6,
            x23 * p5 + *(pxlocal+n*3+2) * t5},
        {0, 0, -(*(pxlocal+n*3+2) * q5), x23 + *(pxlocal+n*3) * r5, x23, x23 * (1 - r5),
            *(pxlocal+n*3) * q5 + *(pxlocal+n*3+2), *(pxlocal+n*3+2),
            -x23 * q5 + *(pxlocal+n*3+2) * (1 - r5)},
        {-4 * (*(pxlocal+n*3+2)), 2 * (*(pxlocal+n*3+2)), *(pxlocal+n*3+2) * (2 - r5),
            -(*(pxlocal+n*3) * q5), 0, x23 * q5, -4 * x23 + *(pxlocal+n*3) * r5, 2 * x23,
            x23 * (2 - r5) + *(pxlocal+n*3+2) * q5},
        {-(*(pxlocal+n*3+2) * p6), *(pxlocal+n*3+2) * p6, *(pxlocal+n*3+2) * p4, 0,
            *(pxlocal+n*3) * t4, -(*(pxlocal+n*3+1) * t4), *(pxlocal+n*3+1) * p6,
            x23 * p6 + *(pxlocal+n*3) * p4,
            -(*(pxlocal+n*3+1) * p4) + *(pxlocal+n*3+2) * t4},
        {0, 0, *(pxlocal+n*3+2) * q4, *(pxlocal+n*3+1),
            *(pxlocal+n*3+1) + *(pxlocal+n*3) * r4, *(pxlocal+n*3+1) * (1 - r4),
            -(*(pxlocal+n*3+2)), -(*(pxlocal+n*3+2)) + *(pxlocal+n*3) * q4,
            *(pxlocal+n*3+2) * (r4 - 1) - *(pxlocal+n*3+1) * q4},
        {-2 * (*(pxlocal+n*3+2)), 4 * (*(pxlocal+n*3+2)), *(pxlocal+n*3+2) * (r4 - 2), 0,
            -(*(pxlocal+n*3) * q4), *(pxlocal+n*3+1) * q4, 2 * (*(pxlocal+n*3+1)),
            -4 * (*(pxlocal+n*3+1)) + *(pxlocal+n*3) * r4,
            *(pxlocal+n*3+1) * (2 - r4) - *(pxlocal+n*3+2) * q4},
        {0, 0, -(*(pxlocal+n*3+2) * (p4 + p5)), *(pxlocal+n*3) * t5,
            -(*(pxlocal+n*3) * t4), -x23 * t5 + *(pxlocal+n*3+1) * t4,
            *(pxlocal+n*3) * p5, -(*(pxlocal+n*3) * p4),
            -x23 * p5 + *(pxlocal+n*3+1) * p4 - *(pxlocal+n*3+2) * (t4 + t5)},
        {0, 0, *(pxlocal+n*3+2) * (q4 - q5), *(pxlocal+n*3) * (r5 - 1),
            *(pxlocal+n*3) * (r4 - 1),
            -x23 * r5 - *(pxlocal+n*3+1) * r4 - *(pxlocal+n*3), *(pxlocal+n*3) * q5,
            *(pxlocal+n*3) * q4,
            -x23 * q5 - *(pxlocal+n*3+1) * q4 + *(pxlocal+n*3+2) * (r4 - r5)},
        {0, 0, *(pxlocal+n*3+2) * (r4 - r5), -(*(pxlocal+n*3) * q5),
            -(*(pxlocal+n*3) * q4), *(pxlocal+n*3+1) * q4 + x23 * q5,
            *(pxlocal+n*3) * (r5 - 2), *(pxlocal+n*3) * (r4 - 2),
            -x23 * r5 - *(pxlocal+n*3+1) * r4 + 4 * (*(pxlocal+n*3)) +
            *(pxlocal+n*3+2) * (q5 - q4)}
    };

    // Compute coupled elasto-plastic element stiffness matrix
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            cdL[i][j * 3] = cdL[i][j * 3 + 1] = cdL[i][j * 3 + 2] = cd[i][j] / 6;
        }
    }
    for (i = 0; i < 6; ++i) {
        for (j = 0; j < 9; ++j) {
            sum = 0;
            for (k = 0; k < 3; ++k) {
                sum += Bm[k][i] * cdL[k][j];
            }
            Bm_cdL[i][j] = sum;
        }
    }
    for (i = 0; i < 6; ++i) {
        for (j = 0; j < 9; ++j) {
            sum = 0;
            for (k = 0; k < 9; ++k) {
                sum += Bm_cdL[i][k] * alpha_T[j][k];
            }
            *(pkm_mb_sh+i*9+j) = sum;
        }
    }
}

void mass_sh (double *psm, double *pcarea, double *pdens, double *pthick,
              double *pfarea, double *pslength, 
              double *px, long *pminc, long *pmcode, double *pjac)
{
	long i, j, k, l, m, ie, je, ptr, ptr2, ptr3;

    ptr = NE_TR + NE_FR;
    ptr2 = NE_TR * 2 + NE_FR * 2;
    ptr3 = NE_TR + NE_FR * 3;
    
    double el12[3], el23[3], el31[3], normal[3], Mtot;
    
	double m_sh[18][18]; // General element mass matrix 
	
	for (i = 0; i < NE_SH; ++i) {
		
        // Compute element side-lengths
        j = *(pminc+ptr2+i*3+0) - 1;
        k = *(pminc+ptr2+i*3+1) - 1;
        l = *(pminc+ptr2+i*3+2) - 1;
        
        /* Note that the order is reversed to be from Vertex-1 to Vertex-3 so that localz
         is properly oriented */
        for (m = 0; m < 3; ++m) {
            el23[m] = *(px+l*3+m) - *(px+k*3+m);
            el31[m] = *(px+l*3+m) - *(px+j*3+m);
            el12[m] = *(px+k*3+m) - *(px+j*3+m);
        }
        *(pslength+i*3+1) = sqrt(dot(el23,el23,3));
        *(pslength+i*3+2) = sqrt(dot(el31,el31,3));
        *(pslength+i*3) = sqrt(dot(el12,el12,3));
        
        // Compute element face area
        cross(el12,el31,normal,0); // Compute vector normal to shell surface
        *(pfarea+i) = 0.5 * sqrt(dot(normal,normal,3));
        
		// Initialize element mass array to zero
		for (j = 0; j < 18; ++j) {
			for (k = 0; k < 18; k++) {
				m_sh[j][k] = 0;
			}
		}
        
        // Assemble lumped element mass matrix (ADINA theory manual Sec 2.7.7)
        
        // Translational degrees of freedom
        Mtot = (*(pdens+i) * (*(pfarea+i)) * (*(pthick+i)));
        m_sh[0][0] = m_sh[1][1] = m_sh[2][2] = Mtot/3;
        m_sh[6][6] = m_sh[7][7] = m_sh[8][8] = Mtot/3;
        m_sh[12][12] = m_sh[13][13] = m_sh[14][14] = Mtot/3;
        
        // Rotational degrees of freedom
        m_sh[3][3] = m_sh[4][4] = m_sh[5][5] = Mtot/3*pow(*(pthick+i),2)/12;
        m_sh[9][9] = m_sh[10][10] = m_sh[11][11] = Mtot/3*pow(*(pthick+i),2)/12;
        m_sh[15][15] = m_sh[16][16] = m_sh[17][17] = Mtot/3*pow(*(pthick+i),2)/12;
        
		// Assemble system mass matrix - full order [NEQ][NEQ]
        if (SLVFLAG == 0) {
            /* Initialize index and then assign element mass components to structure mass array by index and mcode */
            for (ie = 0; ie < 18; ++ie) {
                for (je = 0; je < 18; ++je) {
                    j = *(pmcode+i*18+ie);
                    k = *(pmcode+i*18+je);
                    /* Add current element mass component to previous elements' components to the given DOFs */
                    if (j != 0) {
                        if (j == k) {
                            *(psm+j-1) += m_sh[ie][je];
                        }
                    }
                }
            }
        }
        else {
            /*Build the full order (i.e. [NEQ][NEQ] mass mastrix using mcode*/
            for (ie = 0; ie < 18; ++ie) {
                for (je = 0; je < 18; ++je) {
                    j = *(pmcode+i*18+ie);
                    k = *(pmcode+i*18+je);
                    
                    if ((j != 0) && (k != 0)) {
                        *(psm+(j-1)*NEQ+k-1) += m_sh[je][ie];
                    }
                }
            }
        }
    }
}

int forces_sh (double *pf_temp, double *pef_ip, double *pef_i, double *pefN_temp,
    double *pefM_temp, double *pdd, double *pd_temp, double *pchi_temp, double *px_temp,
    double *px_ip, double *pemod, double *pnu, double *pxlocal, double *pthick,
    double *pfarea, double *pdeffarea_ip, double *pslength, double *pdefslen_ip,
    double *pyield, double *pc1_ip, double *pc2_ip, double *pc3_ip, double *pc1_i,
    double *pc2_i, double *pc3_i, long *pminc, long *pmcode, long *pjcode)
{
    // Initialize function variables
    long i, j, k, n, ptr, ptr2, ptr3, ptr4, ptr5;
    double sum, sum2;
    int flagm[6] = {0, 1, 6, 7, 12, 13}; // Array of membrane DOFs
    int flagb[9] = {2, 3, 4, 8, 9, 10, 14, 15, 16}; // Array of plate bending DOFs

    // Temporary element force vector in updated local coordinate system
    double ef_temp[18];
    double def[18]; // Element incremental force vector in local coordinate system

    // Transformation matrices
    double T_ip[18][18]; // From global to previous configuration
    double T_i[18][18]; // From global to current configuration

    /* General element stiffness matrix in local coordinate system; may serve as ke or
       kep + kg */
    double k_sh[18][18];
    double ke_m_sh[6][6], ke_b_sh[9][9]; // Membrane and plate bending stiffness matrices

    // Various displacement measures
    // Element local nodal coordinates from current and previous configurations
    double xlocal_i[3], xlocal_ip[3];
    double D[18]; // Element total nodal displacements in global coordinate system
    double d[18]; // Element total nodal displacements in local coordinate system
    double DD[18]; // Element incremental nodal displacements in global coordinate system
    double dd[18]; // Element total nodal displacements in local coordinate system
    /* Element total membrane nodal displacements in body-attached coordinate
       system */
    double dm[6];
    /* Element incremental membrane nodal displacements in body-attached
       coordinate system */
    double ddm[6];
    /* Element incremental plate bending nodal displacements in local coordinate
       system */
    double ddb[9];

    // Ivanov's yield criteria parameters
    double No; // Uniaxial yield force per unit width
    double alpha[3]; // Partial plastification factor
    double Me[3]; // Modified uniaxial yield moment per unit width
    double Nbar[3], Mbar[3], MNbar[3]; // Quadratic stress intensities
    // Factors for computation of Ivanov's yield criteria
    double q_fact[3], r_fact[3], s_fact[3];
    int h_fact[3];
    double phi[3]; // Values of yield function at Vertices 1, 2, and 3
    int yv; // Yielded vertex; "0" indicates no yielded vertex
    // Factors for computation of plastic flow directions
    double c_fact, g_fact, d_fact;
    double gradN_Nbar[3], gradM_Mbar[3]; // Gradients of quadratic stress intensities
    double C[3][3]; // Plane stress constitutive matrix
    double fn[3], fm[3]; // Plastic flow directions
    // Factors for computation of plastic strain rate multiplier
    double fn_C[3], fm_C[3];
    double j_fact, k_fact, B_fact;
    double df_da, da_dchi; // Derivatives involved with pseudo hardening parameter
    double fnC_strn, fmC_curv;
    double lambda; // Plastic strain rate multiplier
    // Total (elastic plus plastic) strain and curvature increments
    double strn[3], curv[3][3];

    ptr = NE_TR + NE_FR * 3;
    ptr2 = NE_TR + NE_FR;
    ptr3 = NE_TR * 6 + NE_FR * 14;
    ptr4 = NE_TR * 2 + NE_FR * 14;
    ptr5 = NE_TR * 2 + NE_FR * 2;
    for (n = 0; n < NE_SH; ++n) {
        // Initialize all elements to zero
        for (i = 0; i < 18; ++i) {
            for (j = 0; j < 18; ++j) {
                k_sh[i][j] = 0;
                T_ip[i][j] = T_i[i][j] = 0;
            }
            ef_temp[i] = def[i] = 0;
        }

        // Assign non-zero elements of coordinate transformation matrix
        T_ip[0][0] = T_ip[3][3] = T_ip[6][6] = T_ip[9][9] = T_ip[12][12] =
            T_ip[15][15] = *(pc1_ip+ptr+n*3);
        T_ip[0][1] = T_ip[3][4] = T_ip[6][7] = T_ip[9][10] = T_ip[12][13] =
            T_ip[15][16] = *(pc1_ip+ptr+n*3+1);
        T_ip[0][2] = T_ip[3][5] = T_ip[6][8] = T_ip[9][11] = T_ip[12][14] =
            T_ip[15][17] = *(pc1_ip+ptr+n*3+2);
        T_ip[1][0] = T_ip[4][3] = T_ip[7][6] = T_ip[10][9] = T_ip[13][12] =
            T_ip[16][15] = *(pc2_ip+ptr+n*3);
        T_ip[1][1] = T_ip[4][4] = T_ip[7][7] = T_ip[10][10] = T_ip[13][13] =
            T_ip[16][16] = *(pc2_ip+ptr+n*3+1);
        T_ip[1][2] = T_ip[4][5] = T_ip[7][8] = T_ip[10][11] = T_ip[13][14] =
            T_ip[16][17] = *(pc2_ip+ptr+n*3+2);
        T_ip[2][0] = T_ip[5][3] = T_ip[8][6] = T_ip[11][9] = T_ip[14][12] =
            T_ip[17][15] = *(pc3_ip+ptr+n*3);
        T_ip[2][1] = T_ip[5][4] = T_ip[8][7] = T_ip[11][10] = T_ip[14][13] =
            T_ip[17][16] = *(pc3_ip+ptr+n*3+1);
        T_ip[2][2] = T_ip[5][5] = T_ip[8][8] = T_ip[11][11] = T_ip[14][14] =
            T_ip[17][17] = *(pc3_ip+ptr+n*3+2);

        if (ANAFLAG == 1) {
            // Pass control to stiffe_sh function
            stiffe_sh (&k_sh[0][0], pemod, pnu, pxlocal, pthick, pfarea, pslength, ptr2,
                n);

            /* Retrieve element total nodal displacements from generalized total nodal
               displacement vector */
            for (i = 0; i < 18; ++i) {
                D[i] = 0;
                j = *(pmcode+ptr3+n*18+i);
                if (j != 0) {
                    D[i] = *(pd_temp+j-1);
                }
            }

            /* Transform element total nodal displacement vector from global into
               local coordinate system */
            for (i = 0; i < 18; ++i) {
                sum = 0;
                for (j = 0; j < 18; ++j) {
                    sum += T_ip[i][j] * D[j];
                }
                d[i] = sum;
            }

            // Compute element force vector
            for (i = 0; i < 18; ++i) {
                sum = 0;
                for (j = 0; j < 18; ++j) {
                    sum += k_sh[i][j] * d[j];
                }
                *(pef_i+ptr4+n*18+i) = sum;
            }
        } else if (ANAFLAG == 2) {
            // Pass control to stiffe_m_sh function
            stiffe_m_sh (&ke_m_sh[0][0], pemod, pnu, pxlocal, pthick, pfarea, ptr2, n);

            // Pass control to stiffe_b_sh function
            stiffe_b_sh (&ke_b_sh[0][0], pemod, pnu, pxlocal, pthick, pfarea, pslength,
                ptr2, n);

            // Pass control to mem_coord function
            mem_coord (dm, n, 0, *(pminc+ptr5+n*3) - 1, *(pminc+ptr5+n*3+1) - 1,
                *(pminc+ptr5+n*3+2) - 1, px_temp, pc1_i, pc2_i, pc3_i, ptr);

            // Assign element membrane nodal displacements
            dm[5] = dm[2] - *(pxlocal+n*3+2);
            dm[4] = dm[1] - *(pxlocal+n*3+1);
            dm[2] = dm[0] - *(pxlocal+n*3);
            dm[0] = dm[1] = dm[3] = 0;

            /* Retrieve element incremental nodal displacements from generalized
               incremental nodal displacement vectors */
            for (i = 0; i < 18; ++i) {
                DD[i] = 0;
                j = *(pmcode+ptr3+n*18+i);
                if (j != 0) {
                    DD[i] = *(pdd+j-1);
                }
            }

            /* Transform element incremental nodal displacement vector from global into
               local coordinate system and assign element incremental plate bending nodal
               displacement vector */
            for (i = 0; i < 9; ++i) {
                sum = 0;
                for (j = 0; j < 18; ++j) {
                    sum += T_ip[flagb[i]][j] * DD[j];
                }
                ddb[i] = sum;
            }

            // Compute membrane contribution to element total internal force vector
            for (i = 0; i < 6; ++i) {
                sum = 0;
                for (j = 0; j < 6; ++j) {
                    sum += ke_m_sh[i][j] * dm[j];
                }
                *(pef_ip+ptr4+n*18+flagm[i]) = 0;
                ef_temp[flagm[i]] = sum;
            }

            /* Compute plate bending contribution to element incremental internal force
               vector */
            for (i = 0; i < 9; ++i) {
                sum = 0;
                for (j = 0; j < 9; ++j) {
                    sum += ke_b_sh[i][j] * ddb[j];
                }
                def[flagb[i]] = sum;
            }
        } else {
            // Pass control to mem_coord function
            mem_coord (xlocal_i, n, 0, *(pminc+ptr5+n*3) - 1, *(pminc+ptr5+n*3+1) - 1,
                *(pminc+ptr5+n*3+2) - 1, px_temp, pc1_i, pc2_i, pc3_i, ptr);

            // Pass control to mem_coord function
            mem_coord (xlocal_ip, n, 0, *(pminc+ptr5+n*3) - 1, *(pminc+ptr5+n*3+1) - 1,
                *(pminc+ptr5+n*3+2) - 1, px_ip, pc1_ip, pc2_ip, pc3_ip, ptr);

            // Assign incremental membrane displacement vector
            ddm[0] = ddm[1] = ddm[3] = 0;
            ddm[2] = xlocal_i[0] - xlocal_ip[0];
            ddm[4] = xlocal_i[1] - xlocal_ip[1];
            ddm[5] = xlocal_i[2] - xlocal_ip[2];

            /* Retrieve element incremental nodal displacements from generalized
               incremental nodal displacement vectors */
            for (i = 0; i < 18; ++i) {
                DD[i] = 0;
                j = *(pmcode+ptr3+n*18+i);
                if (j != 0) {
                    DD[i] = *(pdd+j-1);
                }
            }

            /* Transform element incremental nodal displacement vector from global into
               local coordinate system and assign element incremental plate bending nodal
               displacement vector */
            for (i = 0; i < 9; ++i) {
                sum = 0;
                for (j = 0; j < 18; ++j) {
                    sum += T_ip[flagb[i]][j] * DD[j];
                }
                ddb[i] = sum;
            }

            // Pass control to strn_curv function
            strn_curv (strn, &curv[0][0], ddm, ddb, pxlocal, pdeffarea_ip, pdefslen_ip,
                n);

            // Compute uniaxial yield force per unit width
            No = *(pyield+ptr2+n) * (*(pthick+n));
            yv = 0; // Reset yielded vertex flag
            // Compute plane stress constitutive matrix
            C[0][2] = C[1][2] = C[2][0] = C[2][1] = 0;
            C[0][0] = C[1][1] = *(pemod+ptr2+n) / (1 - pow(*(pnu+n),2));
            C[0][1] = C[1][0] = *(pemod+ptr2+n) / (1 - pow(*(pnu+n),2)) * (*(pnu+n));
            C[2][2] = *(pemod+ptr2+n) / (1 - pow(*(pnu+n),2)) * (1 - *(pnu+n)) / 2;

            /* Compute plastic strain rate multiplier and Ivanov's yield criteria at
               Vertices 1, 2, and 3 */
            for (i = 0; i < 3; ++i) {
                // Compute modified uniaxial yield moment per unit width
                alpha[i] = 1.0 - 0.4 * exp(-2.6 * sqrt(*(pchi_temp+n*3+i)));
                Me[i] = alpha[i] * 0.25 * (*(pyield+ptr2+n)) * pow(*(pthick+n),2);

                // Compute quadratic stress intensities
                Nbar[i] = pow(*(pefN_temp+n*9+i*3),2) + pow(*(pefN_temp+n*9+i*3+1),2) -
                    *(pefN_temp+n*9+i*3) * (*(pefN_temp+n*9+i*3+1)) +
                    3 * pow(*(pefN_temp+n*9+i*3+2),2);
                Mbar[i] = pow(*(pefM_temp+n*9+i*3),2) + pow(*(pefM_temp+n*9+i*3+1),2) -
                    *(pefM_temp+n*9+i*3) * (*(pefM_temp+n*9+i*3+1)) +
                    3 * pow(*(pefM_temp+n*9+i*3+2),2);
                MNbar[i] = *(pefM_temp+n*9+i*3) * (*(pefN_temp+n*9+i*3)) +
                    *(pefM_temp+n*9+i*3+1) * (*(pefN_temp+n*9+i*3+1)) -
                    0.5 * (*(pefM_temp+n*9+i*3)) * (*(pefN_temp+n*9+i*3+1)) -
                    0.5 * (*(pefM_temp+n*9+i*3+1)) * (*(pefN_temp+n*9+i*3)) +
                    3 * (*(pefM_temp+n*9+i*3+2)) * (*(pefN_temp+n*9+i*3+2));

                // Compute factors for computation of Ivanov's yield criteria
                q_fact[i] = Nbar[i] * pow(Me[i],2) + 0.48 * Mbar[i] * pow(No,2);
                if (q_fact[i] >= 1e-4) {
                    r_fact[i] = sqrt(pow(No,2) * pow(Mbar[i],2) +
                        4 * pow(Me[i],2) * pow(MNbar[i],2));
                    if (r_fact[i] / (2 * pow(Me[i],2) * No) >= 1e-4) {
                        h_fact[i] = 1;
                    } else {
                        h_fact[i] = 0;
                    }
                    s_fact[i] = Nbar[i] * Mbar[i] - pow(MNbar[i],2);
                    if (h_fact[i] == 1) {
                        phi[i] = Nbar[i] / pow(No,2) + 0.5 * Mbar[i] / pow(Me[i],2) -
                            0.25 * s_fact[i] / q_fact[i] +
                            r_fact[i] / (2 * pow(Me[i],2) * No);
                    } else {
                        phi[i] = Nbar[i] / pow(No,2) + 0.5 * Mbar[i] / pow(Me[i],2) -
                            0.25 * s_fact[i] / q_fact[i];
                    }

                    if (phi[i] >= 1 - phitol) {
                        // Compute factors for computation of plastic flow directions
                        c_fact = 1 / pow(No,2) - Mbar[i] / (4 * q_fact[i]) +
                            s_fact[i] * pow(Me[i],2) / (4 * pow(q_fact[i],2));
                        if (h_fact[i] == 1) {
                            g_fact =
                                MNbar[i] * (1 / (4 * q_fact[i]) + 1 / (No * r_fact[i]));
                            d_fact = 1  / (2 * pow(Me[i],2)) -
                                Nbar[i] / (4 * q_fact[i]) +
                                0.12 * pow(No,2) * s_fact[i] / pow(q_fact[i],2) +
                                Mbar[i] * No / (2 * pow(Me[i],2) * r_fact[i]);
                        } else {
                            g_fact = MNbar[i] / (4 * q_fact[i]);
                            d_fact = 1  / (2 * pow(Me[i],2)) -
                                Nbar[i] / (4 * q_fact[i]) +
                                0.12 * pow(No,2) * s_fact[i] / pow(q_fact[i],2);
                        }

                        // Compute gradients of quadratic stress intensities
                        gradN_Nbar[0] = 2 * (*(pefN_temp+n*9+i*3)) -
                            *(pefN_temp+n*9+i*3+1);
                        gradN_Nbar[1] = 2 * (*(pefN_temp+n*9+i*3+1)) -
                            *(pefN_temp+n*9+i*3);
                        gradN_Nbar[2] = 6 * (*(pefN_temp+n*9+i*3+2));
                        gradM_Mbar[0] = 2 * (*(pefM_temp+n*9+i*3)) -
                            *(pefM_temp+n*9+i*3+1);
                        gradM_Mbar[1] = 2 * (*(pefM_temp+n*9+i*3+1)) -
                            *(pefM_temp+n*9+i*3);
                        gradM_Mbar[2] = 6 * (*(pefM_temp+n*9+i*3+2));

                        // Compute plastic flow directions
                        for (j = 0; j < 3; ++j) {
                            fn[j] = c_fact * gradN_Nbar[j] + g_fact * gradM_Mbar[j];
                            fm[j] = g_fact * gradN_Nbar[j] + d_fact * gradM_Mbar[j];
                        }

                        /* Compute factors for computation of plastic strain rate
                           multiplier */
                        for (j = 0; j < 3; ++j) {
                            sum = sum2 = 0;
                            for (k = 0; k < 3; ++k) {
                                sum += fn[k] * C[k][j];
                                sum2 += fm[k] * C[k][j];
                            }
                            fn_C[j] = sum;
                            fm_C[j] = sum2;
                        }
                        sum = sum2 = 0;
                        for (j = 0; j <= 2; ++j) {
                            sum += fn_C[j] * fn[j];
                            sum2 += fm_C[j] * fm[j];
                        }
                        j_fact = *(pthick+n) * sum;
                        k_fact = pow(*(pthick+n),3) * sum2 / 12;
                        B_fact = 2 * sqrt(pow(g_fact,2) * Nbar[i] +
                            pow(d_fact,2) * Mbar[i] + 2 * d_fact * g_fact * MNbar[i]);

                        /* Compute derivatives involved with partial plasticifation
                           parameter */
                        if (h_fact[i] == 1) {
                            df_da = -Mbar[i] / (alpha[i] * pow(Me[i],2)) +
                                s_fact[i] * Nbar[i] * pow(Me[i],2) /
                                (2 * pow(q_fact[i],2) * alpha[i]) -
                                r_fact[i] / (alpha[i] * pow(Me[i],2) * No) +
                                2 * pow(MNbar[i],2) / (alpha[i] * No * r_fact[i]);
                        } else {
                            df_da = -Mbar[i] / (alpha[i] * pow(Me[i],2)) +
                                s_fact[i] * Nbar[i] * pow(Me[i],2) /
                                (2 * pow(q_fact[i],2) * alpha[i]);
                        }
                        if (*(pchi_temp+n*3+i) >= 1e-6) {
                            da_dchi = 0.52 * (*(pemod+ptr2+n)) * (*(pthick+n)) *
                                exp(-2.6 * sqrt(*(pchi_temp+n*3+i))) /
                                (3 * (*(pyield+ptr2+n)) * sqrt(*(pchi_temp+n*3+i)));
                        } else {
                            da_dchi = 0;
                        }

                        // Compute plastic strain rate multiplier
                        fnC_strn = fmC_curv = 0;
                        for (j = 0; j <= 2; ++j) {
                            fnC_strn += fn_C[j] * strn[j];
                            fmC_curv += fm_C[j] * curv[i][j];
                        }
                        fnC_strn *= (*(pthick+n));
                        fmC_curv *= pow(*(pthick+n),3) / 12;
                        lambda = (fnC_strn + fmC_curv) /
                            (j_fact + k_fact - B_fact * df_da * da_dchi);

                        // Compute element equivalent plastic curvature
                        *(pchi_temp+n*3+i) +=
                            sqrt(pow((*(pemod+ptr2+n) * (*(pthick+n))) /
                            (3 * (*(pyield+ptr2+n))),2) * pow(B_fact * lambda,2));

                        // Compute element internal membrane forces and bending moments
                        for (j = 0; j <= 2; ++j) {
                            sum = sum2 = 0;
                            for (k = 0; k <= 2; ++k) {
                                sum += C[j][k] * (strn[k] - lambda * fn[k]);
                                sum2 += C[j][k] * (curv[i][k] - lambda * fm[k]);
                            }
                            *(pefN_temp+n*9+i*3+j) += *(pthick+n) * sum;
                            *(pefM_temp+n*9+i*3+j) += pow(*(pthick+n),3) * sum2 / 12;
                        }
                    } else {
                        // Compute element internal membrane forces and bending moments
                        for (j = 0; j <= 2; ++j) {
                            sum = sum2 = 0;
                            for (k = 0; k <= 2; ++k) {
                                sum += C[j][k] * strn[k];
                                sum2 += C[j][k] * curv[i][k];
                            }
                            *(pefN_temp+n*9+i*3+j) += *(pthick+n) * sum;
                            *(pefM_temp+n*9+i*3+j) += pow(*(pthick+n),3) * sum2 / 12;
                        }
                    }
                } else {
                    // Compute element internal membrane forces and bending moments
                    for (j = 0; j <= 2; ++j) {
                        sum = sum2 = 0;
                        for (k = 0; k <= 2; ++k) {
                            sum += C[j][k] * strn[k];
                            sum2 += C[j][k] * curv[i][k];
                        }
                        *(pefN_temp+n*9+i*3+j) += *(pthick+n) * sum;
                        *(pefM_temp+n*9+i*3+j) += pow(*(pthick+n),3) * sum2 / 12;
                    }
                }

                // Compute modified uniaxial yield moment per unit width
                alpha[i] = 1.0 - 0.4 * exp(-2.6 * sqrt(*(pchi_temp+n*3+i)));
                Me[i] = alpha[i] * 0.25 * (*(pyield+ptr2+n)) * pow(*(pthick+n),2);

                // Compute quadratic stress intensities
                Nbar[i] = pow(*(pefN_temp+n*9+i*3),2) + pow(*(pefN_temp+n*9+i*3+1),2) -
                    *(pefN_temp+n*9+i*3) * (*(pefN_temp+n*9+i*3+1)) +
                    3 * pow(*(pefN_temp+n*9+i*3+2),2);
                Mbar[i] = pow(*(pefM_temp+n*9+i*3),2) + pow(*(pefM_temp+n*9+i*3+1),2) -
                    *(pefM_temp+n*9+i*3) * (*(pefM_temp+n*9+i*3+1)) +
                    3 * pow(*(pefM_temp+n*9+i*3+2),2);
                MNbar[i] = *(pefM_temp+n*9+i*3) * (*(pefN_temp+n*9+i*3)) +
                    *(pefM_temp+n*9+i*3+1) * (*(pefN_temp+n*9+i*3+1)) -
                    0.5 * (*(pefM_temp+n*9+i*3)) * (*(pefN_temp+n*9+i*3+1)) -
                    0.5 * (*(pefM_temp+n*9+i*3+1)) * (*(pefN_temp+n*9+i*3)) +
                    3 * (*(pefM_temp+n*9+i*3+2)) * (*(pefN_temp+n*9+i*3+2));

                // Compute factors for computation of Ivanov's yield criteria
                q_fact[i] = Nbar[i] * pow(Me[i],2) + 0.48 * Mbar[i] * pow(No,2);
                if (q_fact[i] >= 1e-4) {
                    r_fact[i] = sqrt(pow(No,2) * pow(Mbar[i],2) +
                        4 * pow(Me[i],2) * pow(MNbar[i],2));
                    if (r_fact[i] / (2 * pow(Me[i],2) * No) >= 1e-4) {
                        h_fact[i] = 1;
                    } else {
                        h_fact[i] = 0;
                    }
                    s_fact[i] = Nbar[i] * Mbar[i] - pow(MNbar[i],2);

                    /* Compute Ivanov's yield criteria and if value is greater than one,
                       bring element internal membrane forces and bending moments back to
                       yield surface */
                    if (h_fact[i] == 1) {
                        phi[i] = Nbar[i] / pow(No,2) + 0.5 * Mbar[i] / pow(Me[i],2) -
                            0.25 * s_fact[i] / q_fact[i] +
                            r_fact[i] / (2 * pow(Me[i],2) * No);
                    } else {
                        phi[i] = Nbar[i] / pow(No,2) + 0.5 * Mbar[i] / pow(Me[i],2) -
                            0.25 * s_fact[i] / q_fact[i];
                    }
                    if (phi[i] > 1 + 10 * phitol) {
                        return 1;
                    } else if (phi[i] > 1 + phitol) {
                        do {
                            // Compute factors for computation of plastic flow directions
                            c_fact = 1 / pow(No,2) - Mbar[i] / (4 * q_fact[i]) +
                                s_fact[i] * pow(Me[i],2) / (4 * pow(q_fact[i],2));
                            if (h_fact[i] == 1) {
                                g_fact = MNbar[i] * (1 / (4 * q_fact[i]) +
                                    1 / (No * r_fact[i]));
                                d_fact = 1  / (2 * pow(Me[i],2)) -
                                    Nbar[i] / (4 * q_fact[i]) +
                                    0.12 * pow(No,2) * s_fact[i] / pow(q_fact[i],2) +
                                    Mbar[i] * No / (2 * pow(Me[i],2) * r_fact[i]);
                            } else {
                                g_fact = MNbar[i] / (4 * q_fact[i]);
                                d_fact = 1  / (2 * pow(Me[i],2)) -
                                    Nbar[i] / (4 * q_fact[i]) +
                                    0.12 * pow(No,2) * s_fact[i] / pow(q_fact[i],2);
                            }

                            // Compute gradients of quadratic stress intensities
                            gradN_Nbar[0] = 2 * (*(pefN_temp+n*9+i*3)) -
                                *(pefN_temp+n*9+i*3+1);
                            gradN_Nbar[1] = 2 * (*(pefN_temp+n*9+i*3+1)) -
                                *(pefN_temp+n*9+i*3);
                            gradN_Nbar[2] = 6 * (*(pefN_temp+n*9+i*3+2));
                            gradM_Mbar[0] = 2 * (*(pefM_temp+n*9+i*3)) -
                                *(pefM_temp+n*9+i*3+1);
                            gradM_Mbar[1] = 2 * (*(pefM_temp+n*9+i*3+1)) -
                                *(pefM_temp+n*9+i*3);
                            gradM_Mbar[2] = 6 * (*(pefM_temp+n*9+i*3+2));

                            // Compute plastic flow directions
                            for (j = 0; j < 3; ++j) {
                                fn[j] = c_fact * gradN_Nbar[j] + g_fact * gradM_Mbar[j];
                                fm[j] = g_fact * gradN_Nbar[j] + d_fact * gradM_Mbar[j];
                            }

                            /* Compute factors for computation of plastic strain rate
                               multiplier */
                            for (j = 0; j < 3; ++j) {
                                sum = sum2 = 0;
                                for (k = 0; k < 3; ++k) {
                                    sum += fn[k] * C[k][j];
                                    sum2 += fm[k] * C[k][j];
                                }
                                fn_C[j] = sum;
                                fm_C[j] = sum2;
                            }
                            sum = sum2 = 0;
                            for (j = 0; j < 3; ++j) {
                                sum += fn_C[j] * fn[j];
                                sum2 += fm_C[j] * fm[j];
                            }
                            j_fact = *(pthick+n) * sum;
                            k_fact = pow(*(pthick+n),3) * sum2 / 12;
                            B_fact = 2 * sqrt(pow(g_fact,2) * Nbar[i] +
                                pow(d_fact,2) * Mbar[i] +
                                2 * d_fact * g_fact * MNbar[i]);

                            /* Compute derivatives involved with partial plasticifation
                               parameter */
                            if (h_fact[i] == 1) {
                                df_da = -Mbar[i] / (alpha[i] * pow(Me[i],2)) +
                                    s_fact[i] * Nbar[i] * pow(Me[i],2) /
                                    (2 * pow(q_fact[i],2) * alpha[i]) -
                                    r_fact[i] / (alpha[i] * pow(Me[i],2) * No) +
                                    2 * pow(MNbar[i],2) / (alpha[i] * No * r_fact[i]);
                            } else {
                                df_da = -Mbar[i] / (alpha[i] * pow(Me[i],2)) +
                                    s_fact[i] * Nbar[i] * pow(Me[i],2) /
                                    (2 * pow(q_fact[i],2) * alpha[i]);
                            }
                            if (*(pchi_temp+n*3+i) >= 1e-6) {
                                da_dchi = 0.52 * (*(pemod+ptr2+n)) * (*(pthick+n)) *
                                    exp(-2.6 * sqrt(*(pchi_temp+n*3+i))) /
                                    (3 * (*(pyield+ptr2+n)) * sqrt(*(pchi_temp+n*3+i)));
                            } else {
                                da_dchi = 0;
                            }

                            /* Compute plastic strain rate multiplier required to return
                               element internal membrane forces and bending moments to
                               yield surface */
                            lambda = (phi[i] - 1) /
                                (j_fact + k_fact - B_fact * df_da * da_dchi);

                            /* Return element internal membrane forces and bending
                               moments to yield surface */
                            for (j = 0; j < 3; ++j) {
                                sum = sum2 = 0;
                                for (k = 0; k < 3; ++k) {
                                    sum += C[j][k] * (-lambda * fn[k]);
                                    sum2 += C[j][k] * (-lambda * fm[k]);
                                }
                                *(pefN_temp+n*9+i*3+j) += *(pthick+n) * sum;
                                *(pefM_temp+n*9+i*3+j) += pow(*(pthick+n),3) * sum2 / 12;
                            }

                            // Compute modified uniaxial yield moment per unit width
                            alpha[i] = 1.0 - 0.4 * exp(-2.6 * sqrt(*(pchi_temp+n*3+i)));
                            Me[i] = alpha[i] * 0.25 * (*(pyield+ptr2+n)) *
                                pow(*(pthick+n),2);

                            // Compute quadratic stress intensities
                            Nbar[i] = pow(*(pefN_temp+n*9+i*3),2) +
                                pow(*(pefN_temp+n*9+i*3+1),2) -
                                *(pefN_temp+n*9+i*3) * (*(pefN_temp+n*9+i*3+1)) +
                                3 * pow(*(pefN_temp+n*9+i*3+2),2);
                            Mbar[i] = pow(*(pefM_temp+n*9+i*3),2) +
                                pow(*(pefM_temp+n*9+i*3+1),2) -
                                *(pefM_temp+n*9+i*3) * (*(pefM_temp+n*9+i*3+1)) +
                                3 * pow(*(pefM_temp+n*9+i*3+2),2);
                            MNbar[i] = *(pefM_temp+n*9+i*3) * (*(pefN_temp+n*9+i*3)) +
                                *(pefM_temp+n*9+i*3+1) * (*(pefN_temp+n*9+i*3+1)) -
                                0.5 * (*(pefM_temp+n*9+i*3)) * (*(pefN_temp+n*9+i*3+1)) -
                                0.5 * (*(pefM_temp+n*9+i*3+1)) * (*(pefN_temp+n*9+i*3)) +
                                3 * (*(pefM_temp+n*9+i*3+2)) * (*(pefN_temp+n*9+i*3+2));

                            // Compute factors for computation of Ivanov's yield criteria
                            q_fact[i] = Nbar[i] * pow(Me[i],2) +
                                0.48 * Mbar[i] * pow(No,2);
                            if (q_fact[i] >= 1e-4) {
                                r_fact[i] = sqrt(pow(No,2) * pow(Mbar[i],2) +
                                    4 * pow(Me[i],2) * pow(MNbar[i],2));
                                if (r_fact[i] / (2 * pow(Me[i],2) * No) >= 1e-4) {
                                    h_fact[i] = 1;
                                } else {
                                    h_fact[i] = 0;
                                }
                                s_fact[i] = Nbar[i] * Mbar[i] - pow(MNbar[i],2);

                                // Compute Ivanov's yield criteria
                                if (h_fact[i] == 1) {
                                    phi[i] = Nbar[i] / pow(No,2) +
                                        0.5 * Mbar[i] / pow(Me[i],2) -
                                        0.25 * s_fact[i] / q_fact[i] +
                                        r_fact[i] / (2 * pow(Me[i],2) * No);
                                } else {
                                    phi[i] = Nbar[i] / pow(No,2) +
                                        0.5 * Mbar[i] / pow(Me[i],2) -
                                        0.25 * s_fact[i] / q_fact[i];
                                }
                            }
                        } while (phi[i] > 1 + phitol);

                        if (yv == 0) {
                            yv = i + 1;
                        } else {
                            if (phi[i] < phi[yv - 1]) {
                                yv = i + 1;
                            }
                        }
                    } else if (phi[i] >= 1 - phitol) {
                        if (yv == 0) {
                            yv = i + 1;
                        } else {
                            if (phi[i] < phi[yv - 1]) {
                                yv = i + 1;
                            }
                        }
                    }
                } else {
                    phi[i] = 0;
                }
            }

            /* Check if yielding has occured at any vertex and compute appropriate
               element stiffness matrix */
            if (yv == 0) {
                // Pass control to stiffe_m_sh function
                stiffe_m_sh (&ke_m_sh[0][0], pemod, pnu, pxlocal, pthick, pfarea, ptr2,
                    n);

                // Pass control to stiffe_b_sh function
                stiffe_b_sh (&ke_b_sh[0][0], pemod, pnu, pxlocal, pthick, pfarea,
                    pslength, ptr2, n);

                // Pass control to mem_coord function
                mem_coord (dm, n, 0, *(pminc+ptr5+n*3) - 1, *(pminc+ptr5+n*3+1) - 1,
                    *(pminc+ptr5+n*3+2) - 1, px_temp, pc1_i, pc2_i, pc3_i, ptr);

                // Assign element membrane nodal displacements
                dm[5] = dm[2] - *(pxlocal+n*3+2);
                dm[4] = dm[1] - *(pxlocal+n*3+1);
                dm[2] = dm[0] - *(pxlocal+n*3);
                dm[0] = dm[1] = dm[3] = 0;

                /* Retrieve element incremental nodal displacements from generalized
                   incremental nodal displacement vectors */
                for (i = 0; i < 18; ++i) {
                    DD[i] = 0;
                    j = *(pmcode+ptr3+n*18+i);
                    if (j != 0) {
                        DD[i] = *(pdd+j-1);
                    }
                }

                /* Transform element incremental nodal displacement vector from global
                   into local coordinate system */
                for (i = 0; i < 18; ++i) {
                    sum = 0;
                    for (j = 0; j < 18; ++j) {
                        sum += T_ip[i][j] * DD[j];
                    }
                    dd[i] = sum;
                }

                // Assign element incremental plate bending nodal displacement vector
                for (i = 0; i < 9; ++i) {
                    ddb[i] = dd[flagb[i]];
                }

                // Compute membrane contribution to element total internal force vector
                for (i = 0; i < 6; ++i) {
                    sum = 0;
                    for (j = 0; j < 6; ++j) {
                        sum += ke_m_sh[i][j] * dm[j];
                    }
                    *(pef_ip+ptr4+n*18+flagm[i]) = 0;
                    ef_temp[flagm[i]] = sum;
                }

                /* Compute plate bending contribution to element incremental internal
                   force vector */
                for (i = 0; i < 9; ++i) {
                    sum = 0;
                    for (j = 0; j < 9; ++j) {
                        sum += ke_b_sh[i][j] * ddb[j];
                    }
                    def[flagb[i]] = sum;
                }
            } else {
                // Pass control to stiffm_sh function
                stiffm_sh (&k_sh[0][0], pemod, pnu, pxlocal, pthick, pdeffarea_ip,
                    pdefslen_ip, pyield, pefN_temp, pefM_temp, pchi_temp, yv-1, &No,
                    &alpha[yv-1], &Me[yv-1], &Nbar[yv-1], &Mbar[yv-1], &MNbar[yv-1],
                    &q_fact[yv-1], &r_fact[yv-1], &s_fact[yv-1], &h_fact[yv-1], ptr2, n);

                /* Transform element total nodal displacement vector from global into
                   local coordinate system */
                for (i = 0; i < 18; ++i) {
                    sum = 0;
                    for (j = 0; j < 18; ++j) {
                        sum += T_ip[i][j] * DD[j];
                    }
                    dd[i] = sum;
                }

                // Compute element force vector
                for (i = 0; i < 18; ++i) {
                    sum = 0;
                    for (j = 0; j < 18; ++j) {
                        sum += k_sh[i][j] * dd[j];
                    }
                    def[i] = sum;
                }
            }
        }

        // Assign non-zero elements of coordinate transformation matrix
        T_i[0][0] = T_i[3][3] = T_i[6][6] = T_i[9][9] = T_i[12][12] = T_i[15][15] =
            *(pc1_i+ptr+n*3);
        T_i[0][1] = T_i[3][4] = T_i[6][7] = T_i[9][10] = T_i[12][13] = T_i[15][16] =
            *(pc1_i+ptr+n*3+1);
        T_i[0][2] = T_i[3][5] = T_i[6][8] = T_i[9][11] = T_i[12][14] = T_i[15][17] =
            *(pc1_i+ptr+n*3+2);
        T_i[1][0] = T_i[4][3] = T_i[7][6] = T_i[10][9] = T_i[13][12] = T_i[16][15] =
            *(pc2_i+ptr+n*3);
        T_i[1][1] = T_i[4][4] = T_i[7][7] = T_i[10][10] = T_i[13][13] = T_i[16][16] =
            *(pc2_i+ptr+n*3+1);
        T_i[1][2] = T_i[4][5] = T_i[7][8] = T_i[10][11] = T_i[13][14] = T_i[16][17] =
            *(pc2_i+ptr+n*3+2);
        T_i[2][0] = T_i[5][3] = T_i[8][6] = T_i[11][9] = T_i[14][12] = T_i[17][15] =
            *(pc3_i+ptr+n*3);
        T_i[2][1] = T_i[5][4] = T_i[8][7] = T_i[11][10] = T_i[14][13] = T_i[17][16] =
            *(pc3_i+ptr+n*3+1);
        T_i[2][2] = T_i[5][5] = T_i[8][8] = T_i[11][11] = T_i[14][14] = T_i[17][17] =
            *(pc3_i+ptr+n*3+2);

        // Store current element force vector
        if (ANAFLAG == 2) {
            /* Construct coordinate transformation matrix which transforms previous
               configuration to current configuration, i.e. T_i * T_ip^T */
            double Ti_Tip[18][18];
            for (i = 0; i < 18; ++i) {
                for (j = 0; j < 18; ++j) {
                    sum = 0;
                    for (k = 0; k < 18; ++k) {
                        sum += T_i[i][k] * T_ip[j][k];
                    }
                    Ti_Tip[i][j] = sum;
                }
            }

            /* Add contribution of incremental element force vector to previous element
               force vector and update reference configuration */
            for (i = 0; i < 18; ++i) {
                sum = 0;
                for (j = 0; j < 18; ++j) {
                    sum += Ti_Tip[i][j] * (def[j] + *(pef_ip+ptr4+n*18+j));
                }
                *(pef_i+ptr4+n*18+i) = ef_temp[i] + sum;
            }
        } else if (ANAFLAG == 3) {
            /* Construct coordinate transformation matrix which transforms previous
               configuration to current configuration, i.e. T_i * T_ip^T */
            double Ti_Tip[18][18];
            for (i = 0; i < 18; ++i) {
                for (j = 0; j < 18; ++j) {
                    sum = 0;
                    for (k = 0; k < 18; ++k) {
                        sum += T_i[i][k] * T_ip[j][k];
                    }
                    Ti_Tip[i][j] = sum;
                }
            }

            if (yv == 0) {
                /* Add contribution of incremental element force vector to previous
                   element force vector and update reference configuration */
                for (i = 0; i < 18; ++i) {
                    sum = 0;
                    for (j = 0; j < 18; ++j) {
                        sum += Ti_Tip[i][j] * (def[j] + *(pef_ip+ptr4+n*18+j));
                    }
                    *(pef_i+ptr4+n*18+i) = ef_temp[i] + sum;
                }
            } else {
                /* Add contribution of incremental element force vector to previous
                   element force vector and update reference configuration */
                for (i = 0; i < 18; ++i) {
                    sum = 0;
                    for (j = 0; j < 18; ++j) {
                        sum += Ti_Tip[i][j] * (def[j] + *(pef_ip+ptr4+n*18+j));
                    }
                    *(pef_i+ptr4+n*18+i) = sum;
                }
            }
        }

        /* Transform element force vector from local into global coordinate system and
           add element contribution to generalized internal force vector */
        for (i = 0; i < 18; ++i) {
            sum = 0;
            for (j = 0; j < 18; ++j) {
                sum += T_i[j][i] * (*(pef_i+ptr4+n*18+j));
            }
            j = *(pmcode+ptr3+n*18+i);
            if (j != 0) {
                *(pf_temp+j-1) += sum;
            }
        }
    }
    return 0;
}

void mem_coord (double *pxlocal, long n, long i, long j, long k, long l, double *px,
    double *pc1, double *pc2, double *pc3, long ptr)
{
    // Initialize function variables
    int m;
    double sum;
    double X[3][3]; // Matrix of translated element global coordinates
    double x[3][3]; // Matrix of element local coordinates
    double T[3][3]; // Coordinate coordinate transformation matrix

    /* Store element global coordinates in matrix, translating Vertex-1 to the
       origin of the local coordinate system */
    for (m = 0; m < 3; ++m) {
        X[m][0] = 0;
        X[m][1] = *(px+k*3+m) - *(px+j*3+m);
        X[m][2] = *(px+l*3+m) - *(px+j*3+m);
    }

    // Assign elements of coordinate coordinate transformation matrix
    T[0][0] = *(pc1+ptr+n*3);
    T[0][1] = *(pc1+ptr+n*3+1);
    T[0][2] = *(pc1+ptr+n*3+2);
    T[1][0] = *(pc2+ptr+n*3);
    T[1][1] = *(pc2+ptr+n*3+1);
    T[1][2] = *(pc2+ptr+n*3+2);
    T[2][0] = *(pc3+ptr+n*3);
    T[2][1] = *(pc3+ptr+n*3+1);
    T[2][2] = *(pc3+ptr+n*3+2);

    // Compute initial element local coordinates
    for (j = 0; j < 3; ++j) {
        for (k = 0; k < 3; ++k) {
            sum = 0;
            for (l = 0; l < 3; ++l) {
                sum += T[j][l] * X[l][k];
            }
            x[j][k] = sum;
        }
    }

    // Assign local membrane coordinates
    *(pxlocal+i*3) = x[0][1];
    *(pxlocal+i*3+1) = x[0][2];
    *(pxlocal+i*3+2) = x[1][2];
}

void strn_curv (double *pstrn, double *pcurv, double *pddm, double *pddb,
    double *pxlocal, double *pdeffarea_ip, double *pdefslen_ip, long n)
{
    // Initialize function variables
    int i, j, k;
    double sum;
    double Bm[3][6]; // Membrane strain-displacement matrix
    // Coefficients for computing plate bending strain-displacement matrix
    double x23, l12, l23, l31, p4, p5, p6, t4, t5, q4, q5, r4, r5;

    // Compute membrane strain-displacement matrix
    Bm[0][1] = Bm[0][3] = Bm[0][4] = Bm[0][5] = Bm[1][0] = Bm[1][2] = Bm[1][4] =
        Bm[2][5] = 0;
    Bm[0][0] = Bm[2][1] = -(*(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n))));
    Bm[0][2] = Bm[2][3] = *(pxlocal+n*3+2) / (2 * (*(pdeffarea_ip+n)));
    Bm[1][1] = Bm[2][0] = (*(pxlocal+n*3+1) - *(pxlocal+n*3)) /
        (2 * (*(pdeffarea_ip+n)));
    Bm[1][3] = Bm[2][2] = -(*(pxlocal+n*3+1) / (2 * (*(pdeffarea_ip+n))));
    Bm[1][5] = Bm[2][4] = *(pxlocal+n*3) / (2 * (*(pdeffarea_ip+n)));

    // Compute total strain increment
    for (i = 0; i < 3; ++i) {
        sum = 0;
        for (j = 0; j < 6; ++j) {
            sum += Bm[i][j] * (*(pddm+j));
        }
        *(pstrn+i) = sum;
    }

    // Define local coordinate coefficients for convenience in defining p, t, q, and r
    x23 = *(pxlocal+n*3) - *(pxlocal+n*3+1);
    l12 = pow(*(pdefslen_ip+n*3),2);
    l23 = pow(*(pdefslen_ip+n*3+1),2);
    l31 = pow(*(pdefslen_ip+n*3+2),2);

    // Compute coefficients for assembly of plate bending strain-displacement matrix
    p4 = -6 * x23 / l23;
    p5 = -6 * (*(pxlocal+n*3+1)) / l31;
    p6 = 6 * (*(pxlocal+n*3)) / l12;
    t4 = 6 * (*(pxlocal+n*3+2)) / l23;
    t5 = -6 * (*(pxlocal+n*3+2)) / l31;
    q4 = -3 * x23 * (*(pxlocal+n*3+2)) / l23;
    q5 = 3 * (*(pxlocal+n*3+1)) * (*(pxlocal+n*3+2)) / l31;
    r4 = 3 * pow(*(pxlocal+n*3+2),2) / l23;
    r5 = 3 * pow(*(pxlocal+n*3+2),2) / l31;

    /* Define matrices representing area coordinate matrices, evaluated at each vertex,
       multiplied by alpha matrix */
    double LL_alpha[3][3][9] =
    {
        {
            {*(pxlocal+n*3+2) * p6, 0, -4 * (*(pxlocal+n*3+2)), -(*(pxlocal+n*3+2) * p6),
                0, -2 * (*(pxlocal+n*3+2)), 0, 0, 0},
            {-(*(pxlocal+n*3) * t5), x23 + *(pxlocal+n*3) * r5, -(*(pxlocal+n*3) * q5),
                0, *(pxlocal+n*3+1), 0, *(pxlocal+n*3) * t5, *(pxlocal+n*3) * (r5 - 1),
                -(*(pxlocal+n*3) * q5)},
            {-(*(pxlocal+n*3+1) * p6) - *(pxlocal+n*3) * p5,
                *(pxlocal+n*3) * q5 + *(pxlocal+n*3+2), -4 * x23 + *(pxlocal+n*3) * r5,
                *(pxlocal+n*3+1) * p6, -(*(pxlocal+n*3+2)), 2 * (*(pxlocal+n*3+1)),
                *(pxlocal+n*3) * p5, *(pxlocal+n*3) * q5, *(pxlocal+n*3) * (r5 - 2)}
        },
        {
            {-(*(pxlocal+n*3+2) * p6), 0, 2 * (*(pxlocal+n*3+2)), *(pxlocal+n*3+2) * p6,
                0, 4 * (*(pxlocal+n*3+2)), 0, 0, 0},
            {0, x23, 0, *(pxlocal+n*3) * t4, *(pxlocal+n*3+1) + *(pxlocal+n*3) * r4,
                -(*(pxlocal+n*3) * q4), -(*(pxlocal+n*3) * t4),
                *(pxlocal+n*3) * (r4 - 1), -(*(pxlocal+n*3) * q4)},
            {-x23 * p6, *(pxlocal+n*3+2), 2 * x23, x23 * p6 + *(pxlocal+n*3) * p4,
                -(*(pxlocal+n*3+2)) + *(pxlocal+n*3) * q4,
                -4 * (*(pxlocal+n*3+1)) + *(pxlocal+n*3) * r4, -(*(pxlocal+n*3) * p4),
                *(pxlocal+n*3) * q4, *(pxlocal+n*3) * (r4 - 2)}
        },
        {
            {*(pxlocal+n*3+2) * p5, -(*(pxlocal+n*3+2) * q5),
                *(pxlocal+n*3+2) * (2 - r5), *(pxlocal+n*3+2) * p4,
                *(pxlocal+n*3+2) * q4, *(pxlocal+n*3+2) * (r4 - 2),
                -(*(pxlocal+n*3+2) * (p4 + p5)), *(pxlocal+n*3+2) * (q4 - q5),
                *(pxlocal+n*3+2) * (r4 - r5)},
            {x23 * t5, x23 * (1 - r5), x23 * q5, -(*(pxlocal+n*3+1) * t4),
                *(pxlocal+n*3+1) * (1 - r4), *(pxlocal+n*3+1) * q4,
                -x23 * t5 + *(pxlocal+n*3+1) * t4,
                -x23 * r5 - *(pxlocal+n*3+1) * r4 - *(pxlocal+n*3),
                *(pxlocal+n*3+1) * q4 + x23 * q5},
            {x23 * p5 + *(pxlocal+n*3+2) * t5, -x23 * q5 + *(pxlocal+n*3+2) * (1 - r5),
                x23 * (2 - r5) + *(pxlocal+n*3+2) * q5,
                -(*(pxlocal+n*3+1) * p4) + *(pxlocal+n*3+2) * t4,
                *(pxlocal+n*3+2) * (r4 - 1) - *(pxlocal+n*3+1) * q4,
                *(pxlocal+n*3+1) * (2 - r4) - *(pxlocal+n*3+2) * q4,
                -x23 * p5 + *(pxlocal+n*3+1) * p4 - *(pxlocal+n*3+2) * (t4 + t5),
                -x23 * q5 - *(pxlocal+n*3+1) * q4 + *(pxlocal+n*3+2) * (r4 - r5),
                -x23 * r5 - *(pxlocal+n*3+1) * r4 + 4 * (*(pxlocal+n*3)) +
                *(pxlocal+n*3+2) * (q5 - q4)}
        }
    };

    // Compute total curvature increment
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            sum = 0;
            for (k = 0; k < 9; ++k) {
                sum += LL_alpha[i][j][k] * (*(pddb+k));
            }
            *(pcurv+i*3+j) = sum / (2 * (*(pdeffarea_ip+n)));
        }
    }
}
