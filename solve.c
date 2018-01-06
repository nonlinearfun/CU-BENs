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
#include <math.h>
#include <stdlib.h>
#include "prototypes.h"

// CLAPACK header files
#include <Accelerate/Accelerate.h>

extern long NJ, SNDOF, FNDOF, NEQ, NPDB, NTSTPS, NE_SBR, NE_FBR;
extern double dt, ttot;
extern int ANAFLAG, ALGFLAG, SLVFLAG, FSIFLAG, brFSI_FLAG, shFSI_FLAG, CHKPT, RFLAG;
extern FILE *IFP[4], *OFP[8];


int solve (long *pjcode, double *pss, double *pss_fsi, double *psm, double *psm_fsi, double *psd_fsi, double *pr, double *pdd, long *pmaxa, double *pssd, int *pdet,
           double *pum, double *pvm, double *pam, double *puc, double *pvc, double *pac, double *pqdyn, double *ptstps,
           double *pKeff, double *pReff, double *pMeff, double alpham, double alphaf, int *pipiv, int fact, double ddt, double *ppdisp, long *pkht, int *piter, int *pii, int *pij, int tstp)
{
    
    // Initialize function variables
    long i, j, k;
    int err, dum = 0;
    char trans = 'N';
    double time, sum = 0;
    double a0, a1, a2, a3, a4, a5, a6, a7;
    double alpha, delta;
    double dt_temp;// Variable for time stepping scheme
    
    // Initialize CLAPACK variables
    int n, lda, ldb, info, nrhs = 1;
    
    n = lda = ldb = NEQ;
    
    // Pass residual array to the incremental displacements array
    for (i = 0; i < NEQ; ++i) {
        *(pdd+i) = *(pr+i);
    }
    
    // Static analysis
    if (ALGFLAG < 4) {
        
        // Skyline solver for structural elements
        if (SLVFLAG == 0) {
            skyfact(pmaxa, pss, pssd, pdd, fact, pdet);
            err = skysolve (pmaxa, pss, pssd, pdd, fact, pdet);
        }
        
        // CLAPACK direct solver
        else if (SLVFLAG == 1) {
            
            // Call LAPACK routine for solving general matrices
            dgesv_(&n, &nrhs, pss, &lda, pipiv, pr, &ldb, &info);
            
            // Set displacements from CLAPACK to dd array
            for (i = 0; i < NEQ; ++i) {
                *(pdd+i) = *(pr+i);
            }
            
            // Pass control to output function
            output (pr, &dum, pdd, puc, 1);
        }
    }
    
    // Dynamic analysis
    else if (ALGFLAG > 3) {
        
        // Initialize function variables
        int m, n, lda, ldb, info, nrhs = 1;
        m = n = lda = ldb = NEQ;
        char trans = 'N';
        
        // Initialize effective stiffness matrix to zero
        if (SLVFLAG == 1) {
            for (i = 0; i < NEQ; ++i) {
                for (j = 0; j < NEQ; ++j) {
                    *(pKeff+i*NEQ+j) = 0;
                }
            }
        }
        else if (SLVFLAG == 0) {
            for (i = 0; i < *(pmaxa+NEQ)-1; ++i) {
                *(pKeff+i) = 0;
            }
        }
        
        alpha = pow(1-alpham+alphaf, 2)/4;
        delta = 0.5-alpham+alphaf;
        dt_temp = ddt * dt;
        
        // Calculate integration constants
        a0 = 1/(alpha*pow(dt_temp,2));
        a1 = delta/(alpha*dt_temp);
        a2 = 1/(alpha*dt_temp);
        a3 = 1/(2*alpha) - 1;
        a4 = delta/alpha - 1;
        a5 = (dt_temp)/2*(delta/alpha - 2);
        a6 = (dt_temp)*(1-delta);
        a7 = delta*(dt_temp);
        
        /* Calculate effective stiffness matrix */
        if (ANAFLAG == 4) { // FSI analysis, cannot use skyline
            for (i = 0; i < NEQ; ++i) {
                for (j = 0; j < NEQ; ++j) {
                    if (i != j) {
                        *(pKeff+j*NEQ+i) = (*(pss_fsi+i*NEQ+j))+a0*(1-alpham)*(*(psm_fsi+i*NEQ+j))/(1-alphaf);
                    }
                    else {
                        *(pKeff+j*NEQ+i) = (*(pss_fsi+i*NEQ+j))+a0*(1-alpham)*(*(psm_fsi+i*NEQ+j))/(1-alphaf)+a1*(*(psd_fsi+i));
                    }
                }
            }
        }
        else if (ANAFLAG != 4) { // Non-FSI analysis
            if (SLVFLAG == 0) { // using skyline function
                for (i = 0; i < *(pmaxa+NEQ)-1; ++i) {
                    *(pKeff+i) = *(pss+i); // Initialize Keff w K
                }
                for (i = 0; i <  NEQ; ++i) { // Add mass to diagonal elements
                    k = *(pmaxa+i);
                    *(pKeff+k-1) += a0*(1-alpham)*(*(psm+i))/(1-alphaf);
                }
            }
            else if (SLVFLAG == 1) { // using CLAPACK solver
                for (i = 0; i < NEQ; ++i) {
                    for (j = 0; j < NEQ; ++j) {
                        *(pKeff+j*NEQ+i) = (*(pss+i*NEQ+j))+a0*(1-alpham)*(*(psm+i*NEQ+j))/(1-alphaf);
                    }
                }
            }
        }

        if (NPDB == 0) {
            // Factorize Keff
            if (SLVFLAG == 0) {
                skyfact(pmaxa, pKeff, pssd, pdd, fact, pdet);
            }
            else if (SLVFLAG == 1) {
                dgetrf_(&m, &n, pKeff, &lda, pipiv, &info);
            }
        }
        
        if (ALGFLAG == 4){ // Dynamic: linear Newmark Intergration Method
            
            char CblasRowMajor, CblasNoTrans;
            // Initialize a copy of effective stiffness matrix for prescribed displacement matrix computation
            double *pKeffcp = alloc_dbl (*(pmaxa+NEQ)-1);
            if (pKeffcp == NULL) {
                // Pass control to closeio function
                return closeio(1);
            }
            // Initialize a dummy variable to temporarily store the product of matvec operations
            double *pddn = alloc_dbl (NEQ);
            if (pddn == NULL) {
                // Pass control to closeio function
                return closeio(1);
            }
            
            char file[20];
            
            if (RFLAG == 1) {
                // Open I/O for business!
                do {
                    IFP[3] = fopen("results8.txt", "r"); // Open input file for reading
                } while (IFP[3] == 0);
                
                fscanf(IFP[3], "%d\n", &tstp);
                
                // Read in displacements, velocities, and accelerations from checkpoint file
                for (i = 0; i < NEQ; ++i) {
                    fscanf(IFP[3],  "%le,%le,%le\n", &pum[i], &pvm[i], &pam[i]);
                }
                
                fscanf(IFP[3],  "%ld,%ld,%ld\n", &i, &j, &k);
                
                if (i == 0 && j == 0 && k == 0) {
                    printf ("Read in displacements, velocities, and accelerations complete\n");
                }
                
                if (tstp+1 == NTSTPS) {
                    // Pass control to output function
                    time = tstp*dt;
                    output (&time, &dum, pum, pum, 1);
                }
                tstp = tstp+1;
            }
            
            // Initialize current u, v, and a
            for (i = 0; i < NEQ; ++i) {
                *(puc+i) = *(pvc+i) = *(pac+i) = 0;
            }
            
            // Loop through each time step to solve for displacements
            for (k = tstp; k < NTSTPS; ++k) {
                
                // Initialize Meff and Reff for new time step to zero
                // Update displacement vector
                for (i = 0; i < NEQ; ++i) {
                    *(pMeff+i) = 0;	*(pReff+i) = 0;
                    *(pdd+i) = *(pum+i);
                }
                
                // Calculate effective mass matrix
                if (ANAFLAG == 4) { // FSI analysis, cannot use skyline
                    for (i = 0; i < NEQ; ++i) {
                        sum = 0;
                        for (j = 0; j < NEQ; ++j) {
                            sum += *(psm_fsi+i*NEQ+j)*((1-alpham)*((*(pum+j))*a0+(*(pvm+j))*a2+(*(pam+j))*a3)-alpham*(*(pam+j)))/(1-alphaf);
                        }
                        *(pMeff+i) = sum;
                    }
                }
                else if (ANAFLAG != 4) { // Non-FSI analysis
                    if (SLVFLAG == 0) { // using skyline function
                        for (i = 0; i < NEQ; ++i) {
                            *(pMeff+i) = *(psm+i)*((1-alpham)*((*(pum+i))*a0+(*(pvm+i))*a2+(*(pam+i))*a3)-alpham*(*(pam+i)))/(1-alphaf);
                        }
                    }
                    else if (SLVFLAG == 1) { // using CLAPACK solver
                        for (i = 0; i < NEQ; ++i) {
                            sum = 0;
                            for (j = 0; j < NEQ; ++j) {
                                sum += *(psm+i*NEQ+j)*((1-alpham)*((*(pum+j))*a0+(*(pvm+j))*a2+(*(pam+j))*a3)-alpham*(*(pam+j)))/(1-alphaf);
                            }
                            *(pMeff+i) = sum;
                        }
                    }
                }
                
                // Calculate effective load vector
                if (k == 0){ //First time step, cannot interpolate external force vectors
                    for (i = 0; i < NEQ; ++i) {
                        *(pReff+i) = (*(pqdyn+i*NTSTPS+k)) + *(pMeff+i);
                    }
                } else {
                    for (i = 0; i < NEQ; ++i) {
                        *(pReff+i) = (*(pqdyn+i*NTSTPS+k))+alphaf/(1-alphaf)*(*(pqdyn+i*NTSTPS+k-1))+ *(pMeff+i);
                    }
                }
                
                // Calculate effective damping vector (re-use Meff)
                // Add effective damping vector to effective load vector
                if (ANAFLAG == 4) {
                    for (i = 0; i < NEQ; ++i) {
                        *(pMeff+i) = *(psd_fsi+i)*((*(pum+i)*a1+(*(pvm+i))*a4+(*(pam+i))*a5)-alphaf/(1-alphaf)*(*(pum+i)));
                        *(pReff+i) += *(pMeff+i);
                    }
                }
                
                // Calculate static force vector if generalized-alpha method applied
                if (alphaf != 0){
                    for (i = 0; i < NEQ; ++i) {
                        *(pddn+i)=0;
                    }
                    
                    if(ANAFLAG != 4 && SLVFLAG == 0){// Non-FSI analysis, using skyline function
                        skymult (pmaxa, pss, pdd, pddn);
                    } else if (ANAFLAG == 4 || SLVFLAG == 1){// FSI analysis, cannot use skyline function
                        double beta, gamma;
                        int incx, incy;
                        incx = incy = 1;
                        beta = 1;
                        gamma = 0;
                        
                        cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, beta, pss, lda, pdd, incx, gamma, pddn, incy);
                    }
                    for (i = 0; i < NEQ; ++i) {
                        *(pReff+i) -= alphaf/(1-alphaf)*(*(pddn+i));
                    }
                }
                
                // Solve for displacements at current time step
                if (NPDB > 0) {
                    for (i = 0; i < *(pmaxa+NEQ)-1; ++i) {
                        *(pKeffcp+i) = *(pKeff+i);
                    }
                    
                    for (i = 0; i < NEQ; ++i) {
                        if (*(ppdisp+i*NTSTPS+k) != 0) {
                            *(pReff+i) = *(ppdisp+i*NTSTPS+k);
                        }
                    }
                    // Partition Keff matrix
                    matpart (pmaxa, pkht, pKeffcp, pReff, pum, pii, pij);
                    // Factorize Keff
                    if (SLVFLAG == 0) {
                        skyfact(pmaxa, pKeffcp, pssd, pdd, fact, pdet);
                    }
                    else if (SLVFLAG == 1) {
                        dgetrf_(&m, &n, pKeffcp, &lda, pipiv, &info);
                    }
                    
                    // Solve for displacements at current time step
                    if (SLVFLAG == 0) {
                        err = skysolve (pmaxa, pKeffcp, pssd, pReff, fact, pdet);
                    }
                    else if (SLVFLAG == 1) {
                        dgetrs_(&trans, &n, &nrhs, pKeffcp, &lda, pipiv, pReff, &ldb, &info);
                    }
                    
                    for (i = 0; i < *(pmaxa+NEQ)-1; ++i) {
                        *(pKeffcp+i) = *(pKeff+i);
                    }
                    
                } else {
                    if (SLVFLAG == 0) {
                        err = skysolve (pmaxa, pKeff, pssd, pReff, fact, pdet);
                    }
                    else if (SLVFLAG == 1) {
                        dgetrs_(&trans, &n, &nrhs, pKeff, &lda, pipiv, pReff, &ldb, &info);
                    }
                }
                
                if (SLVFLAG == 0 || SLVFLAG == 1) {
                    for (i = 0; i < NEQ; ++i) {
                        *(puc+i) = *(pReff+i);
                    }
                }
                
                // Calculate current velocities and accelerations
                for (i = 0; i < NEQ; ++i) {
                    *(pac+i) = a0*(*(puc+i) - *(pum+i)) - a2*(*(pvm+i)) - a3*(*(pam+i));
                    *(pvc+i) = *(pvm+i) + a6*(*(pam+i)) + a7*(*(pac+i));
                }
                
                time = k*dt;
                // Pass control to output function
                output (&time, &dum, puc, puc, 1);
                
                // Assign current u, v, a to be the previous values
                for (i = 0; i < NEQ; ++i) {
                    *(pum+i) = *(puc+i);
                    *(pvm+i) = *(pvc+i);
                    *(pam+i) = *(pac+i);
                }
                
                if (k % CHKPT == 0){
                    sprintf(file, "results8.txt");
                    do {
                        OFP[7] = fopen(file, "w"); // Open output file for writing
                    } while (OFP[7] == 0);
                    
                    // Print the time step of occurrence
                    fprintf(OFP[7], "%ld\n", k);
                    
                    // Print displacements, velocities, and accelerations
                    for (i = 0; i < NEQ; ++i) {
                        fprintf(OFP[7], "%le,%le,%le\n", *(pum+i), *(pvm+i), *(pam+i));
                    }
                    
                    // Signal the end of displacements, velocities, and accelerations
                    fprintf (OFP[7], "%d,%d,%d\n", 0, 0, 0);
                    
                    fclose(OFP[7]);
                }
            }
            
            if (pKeffcp != NULL) {
                free (pKeffcp);
                pKeffcp = NULL;
            }
            if (pddn != NULL) {
                free (pddn);
                pddn = NULL;
            }
            
        } else if (ALGFLAG == 5){ //Dynamic: nonlinear Newmark Intergration Method
            
            
            // Initialize Reff for new time step to zero
            for (i = 0; i < NEQ; ++i) {
                *(pReff+i) = 0;
            }
            
            //Compute the equivalent change in dynamic external force vector
            if (ANAFLAG != 4 && SLVFLAG == 0) { // Non-FSI analysis, using skyline funciton
                for (i = 0; i < NEQ; ++i){
                    if (*(ppdisp+i*NTSTPS+tstp) != 0 && *(piter) > 0) {
                        *(pReff+i) = 0;
                    }else if (*(ppdisp+i*NTSTPS+tstp) != 0 && *(piter) == 0){
                        *(pReff+i) = *(pum+i);
                    }else {
                        *(pReff+i) = *(pr+i) + *(psm+i)*((1-alpham)*((*(pvm+i))*a2+(*(pam+i))*a3)-alpham*(*(pam+i)))/(1-alphaf);
                    }
                }
            }else if (ANAFLAG ==4 || SLVFLAG == 1) { //using CLAPACK solver
                for (i = 0; i < NEQ; ++i) {
                    sum = 0;
                    for (j = 0; j < NEQ; ++j) {
                        sum += *(psm+i*NEQ+j);
                    }
                    if (*(ppdisp+i*NTSTPS+tstp) != 0 && *(piter) > 0) {
                        *(pReff+i) = 0;
                    }else if (*(ppdisp+i*NTSTPS+tstp) != 0 && *(piter) == 0){
                        *(pReff+i) = *(pum+i);
                    }else {
                        *(pReff+i) = *(pr+i)+sum*((1-alpham)*((*(pvm+i))*a2+(*(pam+i))*a3)-alpham*(*(pam+i)))/(1-alphaf);
                    }
                }
            }
            
            if (NPDB > 0) {
                // Pass residual array to the incremental displacements array
                for (i = 0; i < NEQ; ++i) {
                    *(pdd+i) = *(pr+i);
                }
                // Partition Keff matrix
                matpart (pmaxa, pkht, pKeff, pReff, pum, pii, pij);
                // Factorize Keff
                if (SLVFLAG == 0) {
                    skyfact(pmaxa, pKeff, pssd, pdd, fact, pdet);
                }
                else if (SLVFLAG == 1) {
                    dgetrf_(&m, &n, pKeff, &lda, pipiv, &info);
                }
            }
            
            /*Compute displacement*/
            // Solve for displacements at each iteration
            if (SLVFLAG == 0) {
                err = skysolve (pmaxa, pKeff, pssd, pReff, fact, pdet);
            }
            else if (SLVFLAG == 1) {
                dgetrs_(&trans, &n, &nrhs, pKeff, &lda, pipiv, pReff, &ldb, &info);
            }
            
            //Pass displacement to main for Newton-Raphson iteration
            for (i = 0; i < NEQ; ++i) {
                *(pdd+i) = *(pReff+i);
            }
        }
    }
    return 0;
}


int skyfact (long *pmaxa, double *pss_temp, double *pssd, double *pdd, int fact, int *pdet)
{
    
    // Initialize function variables
    long i, n, kn, kl, ku, kh, k, ic, klt, j, ki, nd, kk, l;
    double b, c;
    
    /* Initialize determinant sign flag to zero; zero indicates a positive definite
     stiffness matrix */
    *pdet = 0;
    // Perform LDL^t factorization of the stiffness matrix
    if (fact == 0) {
        for (n = 1; n <= NEQ; ++n) {
            kn = *(pmaxa+n-1);
            kl = kn + 1;
            ku = *(pmaxa+n) - 1;
            kh = ku - kl;
            
            if (kh < 0) {
                if (ALGFLAG != 3 && *(pss_temp+kn-1) <= 0) {
                    fprintf(OFP[0], "\n***ERROR*** Non-positive definite stiffness");
                    fprintf(OFP[0], " matrix\n");
                    return 1;
                    
                } else if (ALGFLAG == 3) {
                    *(pssd+n-1) = *(pss_temp+kn-1);
                    if (*(pss_temp+kn-1) == 0) {
                        fprintf(OFP[0], "\n***ERROR*** Singular stiffness matrix\n");
                        return 1;
                    } else if (*(pss_temp+kn-1) > 0 && *pdet != 1) {
                        *pdet = 0;
                    } else {
                        *pdet = 1;
                    }
                }
            } else if (kh > 0) {
                k = n - kh;
                ic = 0;
                klt = ku;
                for (j = 1; j <= kh; ++j) {
                    ic++;
                    klt--;
                    ki = *(pmaxa+k-1);
                    
                    nd = *(pmaxa+k) - ki - 1;
                    if (nd > 0) {
                        if (nd < ic) {
                            kk = nd;
                        } else {
                            kk = ic;
                        }
                        c = 0;
                        for (l = 1; l <= kk; ++l) {
                            c += *(pss_temp+ki-1+l) * (*(pss_temp+klt-1+l));
                        }
                        *(pss_temp+klt-1) -= c;
                    }
                    k++;
                }
                k = n;
                b = 0;
                for (kk = kl; kk <= ku; ++kk) {
                    k--;
                    ki = *(pmaxa+k-1);
                    c = *(pss_temp+kk-1) / *(pss_temp+ki-1);
                    b += c * (*(pss_temp+kk-1));
                    *(pss_temp+kk-1) = c;
                }
                *(pss_temp+kn-1) -= b;
                if (ALGFLAG != 3 && *(pss_temp+kn-1) <= 0) {
                    fprintf(OFP[0], "\n***ERROR*** Non-positive definite stiffness");
                    fprintf(OFP[0], " matrix\n");
                    return 1;
                } else if (ALGFLAG == 3) {
                    *(pssd+n-1) = *(pss_temp+kn-1);
                    if (*(pss_temp+kn-1) == 0) {
                        fprintf(OFP[0], "\n***ERROR*** Singular stiffness matrix\n");
                        return 1;
                    } else if (*(pss_temp+kn-1) > 0 && *pdet != 1) {
                        *pdet = 0;
                    } else {
                        *pdet = 1;
                    }
                }
            } else if (kh == 0) {
                k = n;
                b = 0;
                for (kk = kl; kk <= ku; ++kk) {
                    k--;
                    ki = *(pmaxa+k-1);
                    c = *(pss_temp+kk-1) / *(pss_temp+ki-1);
                    b += c * (*(pss_temp+kk-1));
                    *(pss_temp+kk-1) = c;
                }
                *(pss_temp+kn-1) -= b;
                if (ALGFLAG != 3 && *(pss_temp+kn-1) <= 0) {
                    fprintf(OFP[0], "\n***ERROR*** Non-positive definite stiffness");
                    fprintf(OFP[0], " matrix\n");
                    return 1;
                } else if (ALGFLAG == 3) {
                    *(pssd+n-1) = *(pss_temp+kn-1);
                    if (*(pss_temp+kn-1) == 0) {
                        fprintf(OFP[0], "\n***ERROR*** Singular stiffness matrix\n");
                        return 1;
                    } else if (*(pss_temp+kn-1) > 0 && *pdet != 1) {
                        *pdet = 0;
                    } else {
                        *pdet = 1;
                    }
                }
            }
        }
    }
    return 0;
}

int skysolve (long *pmaxa, double *pss_temp, double *pssd, double *pdd, int fact, int *pdet)
{
    
    // Initialize function variables
    long i, n, kl, ku, kh, k, kk;
    double c;
    
    // Reduce right-hand-side load vector
    for (n = 1; n <= NEQ; ++n) {
        kl = *(pmaxa+n-1) + 1;
        ku = *(pmaxa+n) - 1;
        kh = ku - kl;
        if (kh >= 0) {
            k = n;
            c = 0;
            for (kk = kl; kk <= ku; ++kk) {
                k--;
                c += *(pss_temp+kk-1) * (*(pdd+k-1));
            }
            *(pdd+n-1) -= c;
        }
    }
    
    // Back-substitute
    for (n = 0; n < NEQ; ++n) {
        k = *(pmaxa+n);
        *(pdd+n) /= (*(pss_temp+k-1));
    }
    n = NEQ;
    for (i = 2; i <= NEQ; ++i) {
        kl = *(pmaxa+n-1) + 1;
        ku = *(pmaxa+n) - 1;
        kh = ku - kl;
        if (kh >= 0) {
            k = n;
            for (kk = kl; kk <= ku; ++kk) {
                k--;
                *(pdd+k-1) -= *(pss_temp+kk-1) * (*(pdd+n-1));
            }
        }
        n--;
    }
    return 0;
}

int skymult (long *pmaxa, double *pss_temp, double *pdd, double *pddn)
{
    
    // Initialize function variables
    long i, n, kl, ku, kh, k, kk;
    double c;
    
    //Compute the product of the upper triangle in the stiffness matrix and the displacement vector
    n = NEQ;
    for (i = 2; i <= NEQ; ++i) {
        kl = *(pmaxa+n-1) + 1;
        ku = *(pmaxa+n) - 1;
        kh = ku - kl;
        if (kh >= 0) {
            k = n;
            for (kk = kl; kk <= ku; ++kk) {
                k--;
                *(pddn+k-1) += *(pss_temp+kk-1) * (*(pdd+n-1));
            }
        }
        n--;
    }
    
    // Compute product of the diagonal in the stiffness matrix and the displacement vector
    for (n = 0; n < NEQ; ++n) {
        k = *(pmaxa+n);
        *(pddn+n) += (*(pdd+n))*(*(pss_temp+k-1));
    }
    
    //Compute product of the lower triangle in the stiffness matrix and the displacement vector
    n = NEQ;
    for (i = 1; i <= NEQ; ++i) {
        kl = *(pmaxa+n-1) + 1;
        ku = *(pmaxa+n) - 1;
        kh = ku - kl;
        if (kh >= 0) {
            k = n;
            c = 0;
            for (kk = kl; kk <= ku; ++kk) {
                k--;
                c += *(pss_temp+kk-1) * (*(pdd+k-1));
            }
            *(pddn+n-1) += c;
        }
        n--;
    }
    
    for (i = 0; i < NEQ; ++i) {
        *(pdd+i) = *(pddn+i);
    }
    
    return 0;
}

int matpart (long *pmaxa, long *pkht, double *pss, double *pqtot, double *puc, int *pii, int *pij)
{
    int i, j, k;
    
    long kl, ku, kk;
    
    for (i = 0; i < NPDB; ++i) {
        // Rearrange the matrix such that the known terms in the lower triangle
        // of the matrix are on the rhs
        for (k = 0; k < NEQ-NPDB; ++k) {
            // Search for known terms in each row of the interior dofs
            if (*(pii+k)-(*(pkht+(*(pii+k)))) <= *(pij+i) && *(pii+k) > *(pij+i)) {
                kl = (*(pmaxa+(*(pii+k))+1));
                ku = (*(pij+i))-(*(pii+k))+(*(pkht+(*(pii+k))))+1;
                kk = kl-ku;
                *(pqtot + (*(pii+k))) -= (*(pss+kk-1))*(*(puc+(*(pij+i))));
                *(pss+kk-1) = 0;
            }
        }
        
        // Rearrange the matrix such that the known terms in the upper triangle
        // of the matrix are on the rhs
        for (j = 0; j <= *(pkht+*(pij+i)); ++j) {
            kk = *(pmaxa+(*(pij+i)))-1+j;
            if (j != 0) {
                for (k = 0; k < NPDB-1; ++k) {
                    if ((*(pij+i))-j != *(pij+k)){
                        *(pqtot+(*(pij+i))-j) -= (*(pss+kk))*(*(puc+(*(pij+i))));
                    }
                    *(pss+kk) = 0;
                }
            } else if (j == 0) {
                *(pss+kk) = 1;
            }
        }
    }
    
    return 0;
}


