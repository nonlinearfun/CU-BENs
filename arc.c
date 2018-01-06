//********************************************************************************
//**																			**
//**  Pertains to CU-BEN ver 3.1415												**
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

extern long NEQ;
extern int OPTFLAG;
extern FILE *IFP[4], *OFP[8];

int msal (double *pdk, long *pdkdof, long *pjnt, long *pjcode)
{
    // Initialize function variables
    int dir;
    long jt;

    // Read in joint load from input file
    fscanf(IFP[0], "%ld,%d,%lf\n", &jt, &dir, pdk);
    jt = *(pjnt+jt-1) + 1;

    // Write out joint load to output file
    fprintf(OFP[0], "\nInitial Prescribed Displacement:\n\tGlobal Joint\tDirection\t");
    fprintf(OFP[0], "Displacement\n\t%ld\t\t%d\t\t%f\n", jt, dir, *pdk);
    if (OPTFLAG == 2) {
		fprintf(IFP[1], "%ld,%d,%e\n", jt, dir, *pdk);
    }

    *pdkdof = *(pjcode+(jt-1)*7+dir-1); // Scan and load jcode

    // Check that displacement is not prescribed at fixed degree of freedom
    if (*pdkdof == 0) {
        fprintf(OFP[0], "\n*** ERROR *** Initially prescribed displacement at fixed");
        fprintf(OFP[0], " degree of freedom\n");
        return 1;
    }
    (*pdkdof)--;
    return 0;
}

int quad (double *pa, double *pb, double *pc, double *pdt, double *pdp, double *pddr,
          double *pddq, double *pdd, double *pdlpf, double *plpft)
{
    // Initialize function variables
    long i;
    double radical, root1, root2;
    double gamma1, gamma2;
    double *d1 = alloc_dbl(NEQ);
    if (d1 == NULL) {
        return 4;
    }
    double *d2 = alloc_dbl(NEQ);
    if (d2 == NULL) {
        if (d1 != NULL) {
            free (d1);
            d1 = NULL;
        }
        return 4;
    }

    // Check for negative quantity underneath the radical
    radical = (*pb) * (*pb) - 4 * (*pa) * (*pc);
    if (radical > 0) {
        /* Compute each root; then compute numerator of the equation of a dot product
           between displacement increments due to root and displacement increment from
           previous iteration */
        root1 = (- (*pb) + sqrt((*pb) * (*pb) - 4 * (*pa) * (*pc))) / (2 * (*pa));
        gamma1 = 0;
        for (i = 0; i < NEQ; ++i) {
            d1[i] = *(pdt+i) - *(pdp+i) + *(pddr+i) + root1 * (*(pddq+i));
            gamma1 += (*(pdt+i) - *(pdp+i)) * d1[i];
        }
        root2 = (- (*pb) - sqrt((*pb) * (*pb) - 4 * (*pa) * (*pc))) / (2 * (*pa));
        gamma2 = 0;
        for (i = 0; i < NEQ; ++i) {
            d2[i] = *(pdt+i) - *(pdp+i) + *(pddr+i) + root2 * (*(pddq+i));
            gamma2 += (*(pdt+i) - *(pdp+i)) * d2[i];
        }
        /* Choose root from largest positive gamma, i.e. smallest angle between
           displacement increment due to root and displacement increment from previous
           iteration */
        if (gamma1 > gamma2 && gamma1 > 0) {
            for (i = 0; i < NEQ; ++i) {
                *(pdt+i) = d1[i] + *(pdp+i);
                *(pdd+i) = *(pddr+i) + root1 * (*(pddq+i));
            }
            *pdlpf = root1;
            *plpft += *pdlpf;
        } else if (gamma2 > gamma1 && gamma2 > 0) {
            for (i = 0; i < NEQ; ++i) {
                *(pdt+i) = d2[i] + *(pdp+i);
                *(pdd+i) = *(pddr+i) + root2 * (*(pddq+i));
            }
            *pdlpf = root2;
            *plpft += *pdlpf;
        } else {
            if (d1 != NULL) {
                free (d1);
                d1 = NULL;
            }
            if (d2 != NULL) {
                free (d2);
                d2 = NULL;
            }
            return 3;
        }
    } else {
        if (d1 != NULL) {
            free (d1);
            d1 = NULL;
        }
        if (d2 != NULL) {
            free (d2);
            d2 = NULL;
        }
        return 2;
    }

    if (d1 != NULL) {
        free (d1);
        d1 = NULL;
    }
    if (d2 != NULL) {
        free (d2);
        d2 = NULL;
    }
    return 0;
}
