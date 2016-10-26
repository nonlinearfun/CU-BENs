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

extern long NJ, NE_TR, NEQ;
extern int ANAFLAG, ALGFLAG, OPTFLAG, SLVFLAG;
extern FILE *IFP[2], *OFP[5];

void prop_tr (double *px, double *pemod, double *pcarea, double *pdens, double *pllength,
    double *pyield, double *pc1, double *pc2, double *pc3, long *pminc)
{
    // Initialize function variables
    long i, j, k, l;
    double el[3];

    fprintf(OFP[0], "\nTruss Element Properties:\n\tElement\t\tModulus\t\tArea\t\tDensity\t\t");
    fprintf(OFP[0], "Length\t\tYield Stress\n");

    for (i = 0; i < NE_TR; ++i) {
        // Compute element lengths
        j = *(pminc+i*2) - 1;
        k = *(pminc+i*2+1) - 1;
        for (l = 0; l < 3; ++l) {
            el[l] = *(px+k*3+l) - *(px+j*3+l);
        }
        *(pllength+i) = sqrt(dot(el,el,3));

        // Compute direction cosines
        *(pc1+i) = el[0] / *(pllength+i);
        *(pc2+i) = el[1] / *(pllength+i);
        *(pc3+i) = el[2] / *(pllength+i);

        // Read in element properties from input file
        fscanf(IFP[0], "%lf,%lf,%lf,%lf\n", pemod+i, pcarea+i, pdens+i, pyield+i);

        // Write out element properties to output files
        fprintf(OFP[0], "\t%ld\t\t%lf\t%lf\t%lf\t%lf\t%lf\n", i + 1, *(pemod+i), *(pcarea+i), *(pdens+i),
            *(pllength+i), *(pyield+i));

    }

	if (OPTFLAG == 2) {
		for (i = 0; i < NE_TR; ++i) {
			fprintf(IFP[1], "%e,%e,%e\n", *(pemod+i), *(pcarea+i),
				*(pyield+i));
		}
	}
}

void stiff_tr (double *pss, double *pemod, double *pcarea, double *pllength,
    double *pdefllen_ip, double *pyield, double *pc1_ip, double *pc2_ip, double *pc3_ip,
    double *pef_ip, long *pmaxa, long *pmcode)
{
    // Initialize primary function variables
    long i, j, k, n, je, ie;
    double sum;

    double k_tr[2][2]; // General element stiffness matrix in global coordinate system
    double T_ip[2][6]; // Coordinate transformation matrix
    double K_tr[6][6]; // Total element stiffness matrix in global coordinate system

    /* Result of matrix multiplication for transformation of stiffness matrix from local
       to global coordinate system */
    double T_k[6][2];

    // Yield surface parameters
    double Py; // Squash load
    double phitol = 1e-4; // Allowable +/- deviation from 1.0 of phi

    for (n = 0; n < NE_TR; ++n) {
        // Assemble element stiffness matrix
        if (ANAFLAG == 1) {
            k_tr[0][0] = k_tr[1][1] = *(pemod+n) * (*(pcarea+n)) / *(pllength+n);
            k_tr[0][1] = k_tr[1][0] = -(*(pemod+n) * (*(pcarea+n)) / *(pllength+n));
        } else if (ANAFLAG == 2) {
            k_tr[0][0] = k_tr[1][1] = *(pemod+n) * (*(pcarea+n)) *
                pow(*(pdefllen_ip+n),2) / pow(*(pllength+n),3);
            k_tr[0][1] = k_tr[1][0] = -(*(pemod+n) * (*(pcarea+n)) *
                pow(*(pdefllen_ip+n),2) / pow(*(pllength+n),3));
        } else {
            k_tr[0][0] = k_tr[1][1] = *(pemod+n) * (*(pcarea+n)) *
                pow(*(pdefllen_ip+n),2) / pow(*(pllength+n),3);
            k_tr[0][1] = k_tr[1][0] = -(*(pemod+n) * (*(pcarea+n)) *
                pow(*(pdefllen_ip+n),2) / pow(*(pllength+n),3));

            Py = *(pcarea+n) * (*(pyield+n));
            // Pass control to stiffm_tr function
            stiffm_tr (&k_tr[0][0], pef_ip, &Py, n);
        }

        // Assign elements of coordinate transformation matrix
        T_ip[0][3] = T_ip[0][4] = T_ip[0][5] = T_ip[1][0] = T_ip[1][1] = T_ip[1][2] = 0;
        T_ip[0][0] = *(pc1_ip+n);
        T_ip[0][1] = *(pc2_ip+n);
        T_ip[0][2] = *(pc3_ip+n);
        T_ip[1][3] = *(pc1_ip+n);
        T_ip[1][4] = *(pc2_ip+n);
        T_ip[1][5] = *(pc3_ip+n);

        /* Transform the element stiffness matrix from local into global coordinate
           system */
        for (i = 0; i < 6; ++i) {
            for (j = 0; j < 2; ++j) {
                sum = 0;
                for (k = 0; k < 2; ++k) {
                    sum += T_ip[k][i] * k_tr[k][j];
                }
                T_k[i][j] = sum;
            }
        }
        for (i = 0; i < 6; ++i) {
            for (j = 0; j < 6; ++j) {
                sum = 0;
                for (k = 0; k < 2; ++k) {
                    sum += T_k[i][k] * T_ip[k][j];
                }
                K_tr[i][j] = sum;
            }
        }

        /* If 2nd order analysis is being performed, add contribution to global element
           stiffness matrix; this is independent of element orientation (Bathe, 1996) */
        if (ANAFLAG == 2 ||
            (ANAFLAG == 3 && pow(*(pef_ip+n*2) / Py, 2) >= 1 - phitol)) {
                for (i = 0; i < 6; ++i) {
                    K_tr[i][i] += *(pef_ip+n*2) / *(pdefllen_ip+n);
                }
                K_tr[0][3] -= *(pef_ip+n*2) / *(pdefllen_ip+n);
                K_tr[1][4] -= *(pef_ip+n*2) / *(pdefllen_ip+n);
                K_tr[2][5] -= *(pef_ip+n*2) / *(pdefllen_ip+n);
                K_tr[3][0] -= *(pef_ip+n*2) / *(pdefllen_ip+n);
                K_tr[4][1] -= *(pef_ip+n*2) / *(pdefllen_ip+n);
                K_tr[5][2] -= *(pef_ip+n*2) / *(pdefllen_ip+n);
        }
		if (SLVFLAG == 0) {
			/* Initialize index and then assign element tangent stiffness coefficients of
			 element n to the structure stiffness matrix by index, mcode, and maxa */
			for (je = 0; je < 6; ++je) {
				j = *(pmcode+n*6+je);
				if (j != 0) {
					// Check mcode above current entry to find rank of "j"
					for (ie = 0; ie <= je; ++ie) {
						i = *(pmcode+n*6+ie);
						if (i != 0) {
							if (i > j) { // Find element address as diagonal address + delta
								k = *(pmaxa+i-1) + (i - j);
							} else {
								k = *(pmaxa+j-1) + (j - i);
							}
							/* Add current element stiffness to previous elements'
							 contributions to the given DOFs */
							*(pss+k-1) += K_tr[ie][je];
						}
					}
				}
			}
		}
		else {
			/* Build the full order (i.e. [NEQ][NEQ] system stiffness matrix using mcode */
			for (ie = 0; ie < 6; ++ie) {
				for (je = 0; je < 6; ++je) {
					i = *(pmcode+n*6+ie);
					j = *(pmcode+n*6+je);
					if ((i != 0) && (j != 0)) {
						*(pss+(i-1)*NEQ+j-1) += K_tr[je][ie];
					}
                    
				}
			}
		}
    }
}

void stiffm_tr (double *pk_tr, double *pef_ip, double *pPy, long n)
{
    // Yield surface parameters
    double phitol = 1e-4; // Allowable tolerance on phi

    double G[2]; // Vector of yield surface gradients
    double kt_G[2], G_kt_G;

    if (pow(*(pef_ip+n*2) / *pPy, 2) > 1 + phitol) {
        // Assign elements of vector of yield surface gradients
        G[0] = 2 * (*(pef_ip+n*2)) / pow(*pPy,2);
        G[1] = 2 * (*(pef_ip+n*2+1)) / pow(*pPy,2);

        /* The following contains all operations required to compute the element plastic
           reduction matrix */
        kt_G[0] = *(pk_tr) * G[0] + *(pk_tr+1) * G[1];
        kt_G[1] = *(pk_tr+2) * G[0] + *(pk_tr+3) * G[1];
        G_kt_G = kt_G[0] * G[0] + kt_G[1] * G[1];
        *(pk_tr) -= pow(kt_G[0],2) / G_kt_G;
        *(pk_tr+1) -= kt_G[0] * kt_G[1] / G_kt_G;
        *(pk_tr+2) -= kt_G[1] * kt_G[0] / G_kt_G;
        *(pk_tr+3) -= pow(kt_G[1],2) / G_kt_G;
    }
}

void forces_tr (double *pf_temp, double *pef_i, double *pd, double *pemod,
    double *pcarea, double *pllength, double *pdefllen_i, double *pyield, double *pc1_i,
    double *pc2_i, double *pc3_i, long *pmcode)
{
    // Initialize function variables
    long n, i, j, k;
    double strain, sum;

    // Yield surface parameters
    double Py; // Squash load
    double phitol = 1e-4; // Allowable +/- deviation from 1.0 of phi

    if (ANAFLAG == 1) {
        // Total element nodal displacements
        double D[6]; // Global coordinate system
        double d[2]; //  Local coordinate system

        // Element stiffness matrix
        double kk[2][2]; // Local coordinate system
        double K[6][6]; // Global coordinate system

        double T[2][6]; // Transformation matrix (local to global)

        /* Result of matrix multiplication for transformation of stiffness matrix from
           local to global coordinate system */
        double T_kk[6][2];

        for (n = 0; n < NE_TR; ++n) {
            /* Retrieve element incremental nodal displacements from generalized nodal
               displacement vector */
            for (i = 0; i < 6; ++i) {
                D[i] = 0;
                j = *(pmcode+n*6+i);
                if (j != 0) {
                    D[i] = *(pd+j-1);
                }
            }

            // Assign elements of transformation matrix
            T[0][3] = T[0][4] = T[0][5] = T[1][0] = T[1][1] = T[1][2] = 0;
            T[0][0] = T[1][3] = *(pc1_i+n);
            T[0][1] = T[1][4] = *(pc2_i+n);
            T[0][2] = T[1][5] = *(pc3_i+n);

            /* Transform element incremental nodal displacement vector from global into
               local coordinate system */
            for (i = 0; i < 2; ++i) {
                sum = 0;
                for (j = 0; j < 6; ++j) {
                    sum += T[i][j] * D[j];
                }
                d[i] = sum;
            }

            // Assemble element stiffness matrix
            kk[0][0] = kk[1][1] = *(pcarea+n) * (*(pemod+n)) / *(pllength+n);
            kk[0][1] = kk[1][0] = -(*(pcarea+n) * (*(pemod+n)) / *(pllength+n));

            // Compute element force vector
            for (i = 0; i < 2; ++i) {
                sum = 0;
                for (j = 0; j < 2; ++j) {
                    sum += kk[i][j] * d[j];
                }
                *(pef_i+n*2+i) = sum;
            }

            /* Transform element force vector from local into global coordinate system
               and add element contribution to generalized internal force vector */
           for (j = 0; j < 6; ++j) {
                k = *(pmcode+n*6+j);
                if (k != 0) {
                    switch(j) {
                        case(0):
                            *(pf_temp+k-1) -= *(pef_i+n*2) * (*(pc1_i+n));
                            break;
                        case(1):
                            *(pf_temp+k-1) -= *(pef_i+n*2) * (*(pc2_i+n));
                            break;
                        case(2):
                            *(pf_temp+k-1) -= *(pef_i+n*2) * (*(pc3_i+n));
                            break;
                        case(3):
                            *(pf_temp+k-1) -= *(pef_i+n*2+1) * (*(pc1_i+n));
                            break;
                        case(4):
                            *(pf_temp+k-1) -= *(pef_i+n*2+1) * (*(pc2_i+n));
                            break;
                        case(5):
                            *(pf_temp+k-1) -= *(pef_i+n*2+1) * (*(pc3_i+n));
                            break;
                    }
                }
            }
        }
    } else {
        for (n = 0; n < NE_TR; ++n) {
            // Compute element internal force vector
            strain = (*(pdefllen_i+n) - *(pllength+n)) / (*(pllength+n));
            *(pef_i+n*2) = *(pemod+n) * (*(pcarea+n)) * (strain + 0.5 * pow(strain,2)) *
                (*(pdefllen_i+n)) / (*(pllength+n));
            *(pef_i+n*2+1) = -(*(pef_i+n*2));

            /* If element internal force vector exceeds the yield stress, set equal to
               yield stress */
            if (ANAFLAG == 3) {
                Py = *(pyield+n) * (*(pcarea+n));
                if (pow(*(pef_i+n*2) / Py, 2) >= 1 - phitol) {
                    if (*(pef_i+n*2) < 0) {
                        *(pef_i+n*2) = -Py;
                        *(pef_i+n*2+1) = Py;
                    } else {
                        *(pef_i+n*2) = Py;
                        *(pef_i+n*2+1) = -Py;
                    }
                }
            }

            /* Transform element force vector from local into global coordinate system
               and add element contribution to generalized internal force vector */
           for (j = 0; j < 6; ++j) {
                k = *(pmcode+n*6+j);
                if (k != 0) {
                    switch(j) {
                        case(0):
                            *(pf_temp+k-1) -= *(pef_i+n*2) * (*(pc1_i+n));
                            break;
                        case(1):
                            *(pf_temp+k-1) -= *(pef_i+n*2) * (*(pc2_i+n));
                            break;
                        case(2):
                            *(pf_temp+k-1) -= *(pef_i+n*2) * (*(pc3_i+n));
                            break;
                        case(3):
                            *(pf_temp+k-1) -= *(pef_i+n*2+1) * (*(pc1_i+n));
                            break;
                        case(4):
                            *(pf_temp+k-1) -= *(pef_i+n*2+1) * (*(pc2_i+n));
                            break;
                        case(5):
                            *(pf_temp+k-1) -= *(pef_i+n*2+1) * (*(pc3_i+n));
                            break;
                    }
                }
            }
        }
    }
}


void mass_tr (double *psm, double *pcarea, double *pllength, double *pdens, double *px, long *pminc, long *pmcode, double *pjac)
{
	long i, j, k, l,ie, je;
	double el[3];
    
	double m_tr[6][6]; // General element mass matrix 
	
	for (i = 0; i < NE_TR; ++i) {
		
        // Compute element lengths
        j = *(pminc+i*2) - 1;
        k = *(pminc+i*2+1) - 1;
        for (l = 0; l < 3; ++l) {
            el[l] = *(px+k*3+l) - *(px+j*3+l);
        }
        *(pllength+i) = sqrt(dot(el,el,3));
        
		// Initialize element mass array to zero
		for (j = 0; j < 6; ++j) {
			for (k = 0; k < 6; k++) {
				m_tr[j][k] = 0;
			}
		}
        
        // Assemble element mass matrix (Bathe "FE Procedures" (1996) eqn 4.25 - Consistent mass matrix)
        m_tr[0][0] = (*(pdens+i) * (*(pcarea+i)) * (*(pllength+i)))/3;
        m_tr[1][1] = m_tr[2][2] = m_tr[3][3] = m_tr[4][4] = m_tr[5][5] = m_tr[0][0];
        
        m_tr[3][0] = (*(pdens+i) * (*(pcarea+i)) * (*(pllength+i)))/6;
        m_tr[4][1] = m_tr[5][2] = m_tr[0][3] = m_tr[1][4] = m_tr[2][5] = m_tr[3][0];
        
		
		// Assemble system mass array - size [lss] - for use in skyline solver or system mass matrix - full order [NEQ][NEQ]
        if (SLVFLAG == 0) {
            /* Initialize index and then assign element mass components to structure mass array by index and mcode */
            for (ie = 0; ie < 6; ++ie) {
                for (je = 0; je < 6; ++je) {
                    j = *(pmcode+i*6+ie);
                    k = *(pmcode+i*6+je);
                    /* Add current element mass component to previous elements' components to the given DOFs by summing across rows of mass matrix */
                    if ((j != 0) && (k != 0)) {
                        *(psm+j-1) += m_tr[ie][je];
                    }
                }
            }
        }
        else {
            /*Build the full order (i.e. [NEQ][NEQ] mass mastrix using mcode*/
            for (ie = 0; ie < 6; ++ie) {
                for (je = 0; je < 6; ++je) {
                    j = *(pmcode+i*6+ie);
                    k = *(pmcode+i*6+je);
                    
                    if ((j != 0) && (k != 0)) {
                        *(psm+(j-1)*NEQ+k-1) += m_tr[je][ie];
                    }
                }
            }
        }
    }
}
