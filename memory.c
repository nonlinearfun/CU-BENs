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
#include "prototypes.h"

extern FILE *OFP[8];

int * alloc_int (long arraylen)
{
    int *a;
    a = (int *) malloc(arraylen * sizeof(int));
    if (a == NULL) {
        fprintf(OFP[0], "\n***ERROR*** Unable to allocate memory\n");
        return NULL;
    }
    return a;
}

long * alloc_long (long arraylen)
{
    long *a;
    a = (long *) malloc(arraylen * sizeof(long));
    if (a == NULL) {
        fprintf(OFP[0], "\n***ERROR*** Unable to allocate memory\n");
        return NULL;
    }
    return a;
}

double * alloc_dbl (long arraylen)
{
    double *a;
    a = (double *) malloc(arraylen * sizeof(double));
    if (a == NULL) {
        fprintf(OFP[0], "\n***ERROR*** Unable to allocate memory\n");
        return NULL;
    }
    return a;
}

int free_all (int **pp2p2i, int ni, long **pp2p2l, int nl, double **pp2p2d, int nd,
    int flag)
{
    // Initialize function variables
    int i;

    // Free allocated memory for arrays of type int
    for (i = 0; i < ni; ++i) {
        if (*(pp2p2i+i) != NULL) {
            free (*(pp2p2i+i));
            *(pp2p2i+i) = NULL;
        }
    }

    // Free allocated memory for arrays of type long
    for (i = 0; i < nl; ++i) {
        if (*(pp2p2l+i) != NULL) {
            free (*(pp2p2l+i));
            *(pp2p2l+i) = NULL;
        }
    }

    // Free allocated memory for arrays of type double
    for (i = 0; i < nd; ++i) {
        if (*(pp2p2d+i) != NULL) {
            free (*(pp2p2d+i));
            *(pp2p2d+i) = NULL;
        }
    }

    return closeio(flag);
}
