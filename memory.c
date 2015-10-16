//********************************************************************************//
//**																			**//
//**  Pertains to CU-BEN ver 2.0												**//
//**																			**//
//**  Copyright (c) 2015 C. J. Earls                                            **//
//**  Developed by C. J. Earls, Cornell University                              **//
//**  All rights reserved.														**//
//**                                                                            **//
//**  Contributors:                                                             **//
//**    Christopher Stull                                                       **//
//**    Heather Reed                                                            **//
//**    Justyna Kosianka                                                        **//
//**                                                                            **//
//**  Redistribution and use in source and binary forms, with or without		**//
//**  modification, are permitted provided that the following conditions are	**//
//**  met:																		**//
//**																			**//
//**  - Redistributions of source code must retain the above copyright			**//
//**    notice, this list of conditions and the following disclaimer.			**//
//**																			**//
//**  - Redistributions in binary form must reproduce the above copyright		**//
//**    notice, this list of conditions and the following disclaimer listed		**//
//**    in this license in the documentation and/or other materials				**//
//**    provided with the distribution.											**//
//**																			**//
//**  - Neither the name of the copyright holders nor the name of Cornell		**//
//**    University may be used to endorse or promote products derived from		**//
//**    this software without specific prior written permission.				**//
//**																			**//
//**  Private, research, and institutional usage is without charge.				**//
//**  Distribution of modified versions of this soure code is admissible		**//
//**  UNDER THE CONDITION THAT THIS SOURCE CODE REMAINS UNDER COPYRIGHT OF		**//
//**  THE ORIGINAL DEVELOPERS, BOTH SOURCE AND OBJECT CODE ARE MADE FREELY      **//
//**  AVAILABLE WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE MODIFICATIONS.	**//
//**  Distribution of this code as part of a commercial system is permissible	**//
//**  ONLY BY DIRECT ARRANGEMENT WITH THE DEVELOPERS.							**//
//**																			**//
//**																			**//
//**  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS		**//
//**  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT			**//
//**  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR		**//
//**  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT		**//
//**  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,		**//
//**  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT			**//
//**  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,		**//
//**  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY		**//
//**  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT		**//
//**  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE		**//
//**  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.		**//
//**																			**//
//********************************************************************************//


#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"

extern FILE *OFP[7];

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
