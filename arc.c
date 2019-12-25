//********************************************************************************
//**                                                                            **
//**  Pertains to CU-BEN ver 4.0                                                **
//**                                                                            **
//**  CU-BENs: a ship hull modeling finite element library                      **
//**  Copyright (c) 2019 C. J. Earls                                            **
//**  Developed by C. J. Earls, Cornell University                              **
//**  All rights reserved.                                                      **
//**                                                                            **
//**  Contributors:                                                             **
//**    Christopher Stull                                                       **
//**    Heather Reed                                                            **
//**    Justyna Kosianka                                                        **
//**    Wensi Wu                                                                **
//**                                                                            **
//**  This program is free software: you can redistribute it and/or modify it   **
//**  under the terms of the GNU General Public License as published by the     **
//**  Free Software Foundation, either version 3 of the License, or (at your    **
//**  option) any later version.                                                **
//**                                                                            **
//**  This program is distributed in the hope that it will be useful, but       **
//**  WITHOUT ANY WARRANTY; without even the implied warranty of                **
//**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General  **
//**  Public License for more details.                                          **
//**                                                                            **
//**  You should have received a copy of the GNU General Public License along   **
//**  with this program. If not, see <https://www.gnu.org/licenses/>.           **
//**                                                                            **
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
