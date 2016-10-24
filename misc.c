//********************************************************************************
//**																			**
//**  Pertains to CU-BEN ver 3.0												**
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

extern long NJ, NE_TR, NE_FR, NE_SH, NEQ;
extern int ANAFLAG, ALGFLAG, OPTFLAG;
extern FILE *IFP[2], *OFP[7];

void transform (double *pk, double *pT, double *pK, int n)
{
    // Initialize function variables
    int i, j, k;
    double temp[n][n];
    double sum;

    // Multiply transpose of transformation matrix by element stiffness matrix
    for (i = 0; i < n; i ++) {
        for (j = 0; j < n; j ++) {
            sum = 0;
            for (k = 0; k < n; k ++) {
                sum += (*(pT+k*n+i)) * (*(pk+k*n+j));
            }
            temp[i][j] = sum;
        }
    }

    // Multiply above result by transformation matrix
    for (i = 0; i < n; i ++) {
        for (j = 0; j < n; j ++) {
            sum = 0;
            for (k = 0; k < n; k ++) {
                sum += temp[i][k] * (*(pT+k*n+j));
            }
            *(pK+i*n+j) = sum;
        }
    }
}

void updatc (double *px_temp, double *px_ip, double *pxfr_temp, double *pdd,
    double *pdefllen_i, double *pdeffarea_i, double *pdefslen_i, double *poffset,
    int *posflag, double *pauxpt, double *pc1_i, double *pc2_i, double *pc3_i,
    long *pminc, long *pjcode)
{
    // Initialize function variables
    long i, j, k, l, m, ptr, ptr2;
    double el[3];
    double temp[3], localx[3], localy[3], localz[3];
    double el12[3], el23[3], el31[3], normal[3];

    // Update nodal coordinates during each increment
    for (i = 0; i < NJ; ++i) {
        for (j = 0; j < 3; ++j) {
            k = *(pjcode+i*7+j);
            if (k != 0) {
                *(px_ip+i*3+j) = *(px_temp+i*3+j);
                *(px_temp+i*3+j) += *(pdd+k-1);
            } else {
                *(px_ip+i*3+j) = *(px_temp+i*3+j);
            }
        }
    }

    /* Update element lengths and direction cosines (if applicable) during each increment
       for truss elements */
    for (i = 0; i < NE_TR; ++i) {
        j = *(pminc+i*2) - 1;
        k = *(pminc+i*2+1) - 1;
        el[0] = *(px_temp+k*3) - *(px_temp+j*3);
        el[1] = *(px_temp+k*3+1) - *(px_temp+j*3+1);
        el[2] = *(px_temp+k*3+2) - *(px_temp+j*3+2);
        *(pdefllen_i+i) = sqrt(dot(el,el,3));
        // Update direction cosines
        *(pc1_i+i) = el[0] / *(pdefllen_i+i);
        *(pc2_i+i) = el[1] / *(pdefllen_i+i);
        *(pc3_i+i) = el[2] / *(pdefllen_i+i);
    }

    /* Update element lengths and direction cosines (if applicable) during each increment
       for frame elements */
    ptr = NE_TR * 2;
    for (i = 0; i < NE_FR; ++i) {
        // Update member lengths
        j = *(pminc+ptr+i*2) - 1;
        k = *(pminc+ptr+i*2+1) - 1;
        if (*(posflag+i) == 0) {
            for (l = 0; l < 3; ++l) {
                *(pxfr_temp+i*6+l) = *(px_temp+j*3+l);
                *(pxfr_temp+i*6+3+l) = *(px_temp+k*3+l);
            }
        } else {
            for (l = 0; l < 3; ++l) {
                *(pxfr_temp+i*6+l) = *(px_temp+j*3+l) + *(poffset+i*6+l);
                *(pxfr_temp+i*6+3+l) = *(px_temp+k*3+l) + *(poffset+i*6+3+l);
            }
        }
        el[0] = *(pxfr_temp+i*6+3) - *(pxfr_temp+i*6);
        el[1] = *(pxfr_temp+i*6+3+1) - *(pxfr_temp+i*6+1);
        el[2] = *(pxfr_temp+i*6+3+2) - *(pxfr_temp+i*6+2);
        *(pdefllen_i+NE_TR+i) = sqrt(dot(el,el,3));

        // Update direction cosines
        localx[0] = el[0] / *(pdefllen_i+NE_TR+i);
        localx[1] = el[1] / *(pdefllen_i+NE_TR+i);
        localx[2] = el[2] / *(pdefllen_i+NE_TR+i);
        temp[0] = *(pauxpt+i*3) - *(pxfr_temp+i*6);
        temp[1] = *(pauxpt+i*3+1) - *(pxfr_temp+i*6+1);
        temp[2] = *(pauxpt+i*3+2) - *(pxfr_temp+i*6+2);
        cross(localx,temp,localz,1);
        cross(localz,localx,localy,1);
        for (j = 0; j < 3; ++j) {
            *(pc1_i+NE_TR+i*3+j) = localx[j];
            *(pc2_i+NE_TR+i*3+j) = localy[j];
            *(pc3_i+NE_TR+i*3+j) = localz[j];
        }
    }

    /* Update element side lengths, face areas, and direction cosines (if applicable)
       during each increment for shell elements */
    ptr = NE_TR * 2 + NE_FR * 2;
    ptr2 = NE_TR + NE_FR * 3;
    for (i = 0; i < NE_SH; ++i) {
        // Update element side-lengths
        j = *(pminc+ptr+i*3) - 1;
        k = *(pminc+ptr+i*3+1) - 1;
        l = *(pminc+ptr+i*3+2) - 1;
        /* Note that the order is reversed to be from Vertex-1 to Vertex-3 so that
           localz is properly oriented */
        for (m = 0; m < 3; ++m) {
            el23[m] = *(px_temp+l*3+m) - *(px_temp+k*3+m);
            el31[m] = *(px_temp+l*3+m) - *(px_temp+j*3+m);
            el12[m] = *(px_temp+k*3+m) - *(px_temp+j*3+m);
        }
        *(pdefslen_i+i*3+1) = sqrt(dot(el23,el23,3));
        *(pdefslen_i+i*3+2) = sqrt(dot(el31,el31,3));
        *(pdefslen_i+i*3) = sqrt(dot(el12,el12,3));

        cross(el12,el31,normal,0); // Compute vector normal to shell surface
        // Update element face area
        *(pdeffarea_i+i) = 0.5 * sqrt(dot(normal,normal,3));

        // Update direction cosines
        for (m = 0; m < 3; ++m) {
            localx[m] = el12[m] / *(pdefslen_i+i*3);
            localz[m] = normal[m] / (2 * (*(pdeffarea_i+i)));
        }
        cross(localz,localx,localy,1);
        for (m = 0; m < 3; ++m) {
            *(pc1_i+ptr2+i*3+m) = localx[m];
            *(pc2_i+ptr2+i*3+m) = localy[m];
            *(pc3_i+ptr2+i*3+m) = localz[m];
        }
    }
}

int test (double *pd_temp, double *pdd, double *pf_temp, double *pfp, double *pqtot,
    double *pf_ip, double *pintener1, int *pconvchk, double *ptoldisp, double *ptolforc,
    double *ptolener)
{
    // Initialize function variables
    long i;
    double deltad = 0, totald = 0;
    double unbfi = 0, unbfp = 0;
    double inteneri = 0;
    double c;
    *pconvchk = 0;

    if (*ptoldisp < 1) {
        // Compute norm of displacements
        deltad = dot(pdd, pdd, NEQ);
        totald = dot(pd_temp, pd_temp, NEQ);
        // Check against displacement tolerance
        if (totald != 0) {
            c = sqrt(deltad) / sqrt(totald);
            if (c > *ptoldisp) {
                *pconvchk += 10;
            }
        } else {
            fprintf(OFP[0], "\n***ERROR*** Displacements are zero\n");
            return 1;
        }
    }

    if (*ptolforc < 1) {
        // Compute norm of unbalanced forces
        for (i = 0; i < NEQ; ++i) {
            unbfi += (*(pqtot+i) - *(pf_temp+i)) * (*(pqtot+i) - *(pf_temp+i));
            unbfp += (*(pqtot+i) - *(pfp+i)) * (*(pqtot+i) - *(pfp+i));
        }
        // Check against force tolerance
        if (unbfp != 0) {
            c = sqrt(unbfi) / sqrt(unbfp);
            if (c > *ptolforc) {
               *pconvchk += 100;
            }
        } else {
            fprintf(OFP[0], "\n***ERROR*** Force increment is zero\n");
            return 1;
        }
    }

    if (*ptolener < 1) {
        // Compute incremental internal energy during each iteration
        for (i = 0; i < NEQ; ++i) {
            inteneri += *(pdd+i) * (*(pqtot+i) - *(pf_ip+i));
        }
        // Check against tolerances
        if (*pintener1 != 0) {
            c = fabs(inteneri / *pintener1);
            if (c > *ptolener) {
                *pconvchk += 1000;
            }
        } else {
            fprintf(OFP[0], "\n***ERROR*** Energy increment is zero\n");
            return 1;
        }
    }
    return 0;
}

double dot (double *pa, double *pb, int n)
{
    // Initialize function variables
    int i;
    double dp = 0;

    for (i = 0; i < n; ++i) {
        dp += *(pa+i) * (*(pb+i));
    }
    return dp;
}

void cross (double *pa, double *pb, double *pc, int flag)
{
    // Initialize function variables
    int i;
    double length;

    // Compute cross product
    *(pc) = *(pa+1) * (*(pb+2)) - *(pa+2) * (*(pb+1));
    *(pc+1) = *(pa+2) * (*(pb)) - *(pa) * (*(pb+2));
    *(pc+2) = *(pa) * (*(pb+1)) - *(pa+1) * (*(pb));

    if (flag == 1) {
        // Compute length and normalize cross product
        length = sqrt(*(pc) * (*(pc)) + *(pc+1) * (*(pc+1)) + *(pc+2) * (*(pc+2)));
        for (i = 0; i < 3; ++i) {
            *(pc+i) /= length;
        }
    }
}

void inverse (double *pGT_k_G, int n)
{
    // Initialize function variables
    int i, j, k;
    double aug[n][2 * n];
    double temp[n];
    double m;

    // Augment matrix to include identity matrix on right hand side
    for (i = 0; i < n; ++i) {
        for (j = 0; j < 2 * n; ++j) {
            if (j < n) {
                aug[i][j] = *(pGT_k_G+n*i+j);
            } else {
                if (j - n == i) {
                    aug[i][j] = 1;
                } else {
                    aug[i][j] = 0;
                }
            }
        }
    }

    // Begin Gauss-Jordan elimination procedure
    for (i = 0; i < n; ++i) {
        // Find first non-zero element in Column "i"
        for (j = i; j < n; ++j) {
            if (aug[i][j] != 0) {
                k = j;
                break;
            }
        }
        /* If first non-zero element in Column "i" is not in Row "i", swap Row "i" for
           Row "k" */
        if (k != i) {
            for (j = 0; j < 2 * n; ++j) {
                temp[j] = aug[i][j];
                aug[i][j] = aug[k][j];
                aug[k][j] = temp[j];
            }
        }
        // Zero out elements in Column "i" above and below Row "i"
        for (j = 0; j < n; ++j) {
            if (j != i) {
                m = aug[j][i] / aug[i][i];
                for (k = 0; k < 2 * n; ++k) {
                    aug[j][k] -= m * aug[i][k];
                }
            }
        }
    }

    /* Output inverted matrix; left hand side of augmented matrix must be identity
       matrix, hence division by "augmented[i][i]" */
    for (i = 0; i < n; ++i) {
        for (j = n; j < 2 * n; ++j) {
            *(pGT_k_G+n*i+j-n) = aug[i][j] / aug[i][i];
        }
    }
}

void output (double *plpf, int *pitecnt, double *pd, double *pef, int flag)
{
    // Initialize function variables
    long i, j, ptr;
    
    if (flag == 0) {
        // Print layout for output of displacement results
        fprintf(OFP[1], "Model Displacements:\n\tLambda\t\tIterations");
        for (i = 0; i < NEQ; ++i) {
            if (i + 1 <= 1000) {
                fprintf(OFP[1], "\tDOF %ld\t", i + 1);
            } else {
                fprintf(OFP[1], "DOF %ld\t", i + 1);
            }
        }

        if (NE_TR > 0) {
            // Print layout for output of truss element force results
            fprintf(OFP[2], "Truss Element Forces (averaged):\n\tLambda\t\tIterations");
            fprintf(OFP[2], "\t");
            for (i = 0; i < NE_TR; ++i) {
                fprintf(OFP[2], "Element %ld\t", i + 1);
            }
        }

        if (NE_FR > 0) {
            // Print layout for output of frame element force results
            fprintf(OFP[3], "Frame Element Forces (averaged):\n\tLambda\t\tIterations");
            fprintf(OFP[2], "\t");
            for (i = 0; i < NE_FR; ++i) {
                fprintf(OFP[3], "Element %ld\t\t\t\t\t\t\t\t\t\t\t\t\t", i + 1);
            }
            fprintf(OFP[3], "\n\t\t\t\t\t");
            for (i = 0; i < NE_FR; ++i) {
                fprintf(OFP[3], "X-Force\t\tY-Force\t\tZ-Force\t\tX-Moment\tY-Moment\t");
                fprintf(OFP[3], "Z-Moment\tBi-Moment\t");
            }
        }

        if (NE_SH > 0) {
            // Print layout for output of shell element force results
            fprintf(OFP[4], "Shell Element Forces (averaged):\n\tLambda\t\tIterations");
            fprintf(OFP[2], "\t");
            for (i = 0; i < NE_SH; ++i) {
                fprintf(OFP[4], "Element %ld\t\t\t\t\t\t\t\t\t\t\t", i + 1);
            }
            fprintf(OFP[4], "\n\t\t\t\t\t");
            for (i = 0; i < NE_SH; ++i) {
                fprintf(OFP[4], "X-Force\t\tY-Force\t\tZ-Force\t\tX-Moment\tY-Moment\t");
                fprintf(OFP[4], "Z-Moment\t");
            }
        }
    } else {
		if (ANAFLAG == 1) {
			*pitecnt = 0;
		}

        // Output iteration completion in terminal window
        if (ALGFLAG == 4 || ALGFLAG == 5) {
            printf("Time = %e complete\n", *plpf);
        }
        else if (ANAFLAG == 1){
            printf("Analysis complete\n");
        }
        else {
            printf("LPF = %e complete (%d)\n", *plpf, *pitecnt);
        }
        
        // Output displacement response
        fprintf(OFP[1], "\n\t%e\t%d\t", *plpf, *pitecnt);
        for (i = 0; i < NEQ; ++i) {
            fprintf(OFP[1], "\t%e", *(pd+i));
        }

        if (NE_TR > 0) {
            // Output truss element forces
            fprintf(OFP[2], "\n\t%e\t%d\t", *plpf, *pitecnt);
            for (i = 0; i < NE_TR; ++i) {
                fprintf(OFP[2], "\t%e", (fabs(*(pef+i*2)) + fabs(*(pef+i*2+1))) / 2);
            }
        }

        if (NE_FR > 0) {
            // Output frame element forces
            ptr = NE_TR * 2;
            fprintf(OFP[3], "\n\t%e\t%d\t", *plpf, *pitecnt);
            for (i = 0; i < NE_FR; ++i) {
                for (j = 0; j < 7; ++j) {
                    fprintf(OFP[3], "\t%e", (fabs(*(pef+ptr+i*14+j)) +
                        fabs(*(pef+ptr+i*14+7+j))) / 2);
                }
            }
        }

        if (NE_SH > 0) {
            // Output shell element forces
            ptr = NE_TR * 2 + NE_FR * 14;
            fprintf(OFP[4], "\n\t%e\t%d\t", *plpf, *pitecnt);
            for (i = 0; i < NE_SH; ++i) {
                for (j = 0; j < 6; ++j) {
                    fprintf(OFP[4], "\t%e", (fabs(*(pef+ptr+i*18+j)) +
                        fabs(*(pef+ptr+i*18+6+j)) + fabs(*(pef+ptr+i*18+12+j))) / 3);
                }
            }
        }
    }
}

int closeio (int flag)
{
    // Initialize function variables
    int i;

    if (flag == 0) {
        // Close the I/O
        if (OPTFLAG == 1) {
			fclose(IFP[0]);
			for (i = 0; i < 5; ++i) {
				fclose(OFP[i]);
			}
        } else {
        	fclose(IFP[0]);
        	fclose(IFP[1]);
			for (i = 0; i < 5; ++i) {
				fclose(OFP[i]);
			}
        }
    } else {
        // Notify user that solution failed
        fprintf(OFP[0], "\nSolution failed\n");

        // Close the I/O
        if (OPTFLAG == 1) {
			fclose(IFP[0]);
			for (i = 0; i < 5; ++i) {
				fclose(OFP[i]);
			}
        } else {
        	fclose(IFP[0]);
        	fclose(IFP[1]);
			for (i = 0; i < 5; ++i) {
				fclose(OFP[i]);
			}
        }
    }
    return 0;
}
