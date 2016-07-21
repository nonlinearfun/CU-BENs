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
#include <math.h>
#include <stdlib.h>
#include "prototypes.h"

// CLAPACK header files
#include <Accelerate/Accelerate.h>

extern long NJ, SNDOF, FNDOF, NEQ, NTSTPS, NE_SBR, NE_FBR;
extern double dt, ttot;
extern int ANAFLAG, ALGFLAG, SLVFLAG, FSIFLAG, brFSI_FLAG, shFSI_FLAG;
extern FILE *IFP[2], *OFP[7];


int solve (long *pjcode, double *pss, double *pss_fsi, double *psm, double *psm_fsi, double *psd_fsi, double *pr, double *pdd, long *pmaxa, double *pssd, int *pdet,
           double *pum, double *pvm, double *pam, double *puc, double *pvc, double *pac, double *pqdyn, double *ptstps,
           double *pKeff, double *pReff, double *pMeff, double alpha, double delta, int *pipiv, int fact)
{
    
    // Initialize function variables
    long i, j, k;
    int err, dum = 0;
    char trans = 'N';
    double time, sum = 0;
    double a0, a1, a2, a3, a4, a5, a6, a7;
    
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
        
        // Calculate integration constants
        a0 = 1/(alpha*pow(dt,2));
        a1 = delta/(alpha*dt);
        a2 = 1/(alpha*dt);
        a3 = 1/(2*alpha) - 1;
        a4 = delta/alpha - 1;
        a5 = dt/2*(delta/alpha - 2);
        a6 = dt*(1-delta);
        a7 = delta*dt;
        
        /* Calculate effective stiffness matrix */
        if (ANAFLAG == 4) { // FSI analysis, cannot use skyline
            for (i = 0; i < NEQ; ++i) {
                for (j = 0; j < NEQ; ++j) {
                    if (i != j) {
                        *(pKeff+j*NEQ+i) = *(pss_fsi+i*NEQ+j) + a0*(*(psm_fsi+i*NEQ+j));
                    }
                    else {
                        *(pKeff+j*NEQ+i) = *(pss_fsi+i*NEQ+j) + a0*(*(psm_fsi+i*NEQ+j)) + a1*(*(psd_fsi+i));
                    }
                }
            }
        }
        else if (ANAFLAG != 4) { // Non-FSI analysis
            if (SLVFLAG == 0) { // using skyline function
                for (i = 0; i < *(pmaxa+NEQ)-1; ++i) {
                    *(pKeff+i) = *(pss_fsi+i); // Initialize Keff w K
                }
                for (i = 0; i <  NEQ; ++i) { // Add mass to diagonal elements
                    k = *(pmaxa+i);
                    *(pKeff+k-1) += *(psm_fsi+i)*a0;
                }
            }
            else if (SLVFLAG == 1) { // using CLAPACK solver
                for (i = 0; i < NEQ; ++i) {
                    for (j = 0; j < NEQ; ++j) {
                        *(pKeff+j*NEQ+i) = *(pss+i*NEQ+j) + a0*(*(psm+i*NEQ+j));
                    }
                }
            }
        }
        
        // Factorize Keff
        if (SLVFLAG == 0) {
            skyfact(pmaxa, pKeff, pssd, pdd, fact, pdet);
        }
        else if (SLVFLAG == 1) {
            dgetrf_(&m, &n, pKeff, &lda, pipiv, &info);
        }
        
        if (ALGFLAG == 4){ // Dynamic: linear Newmark Intergration Method
            
            // Initialize current u, v, and a
            for (i = 0; i < NEQ; ++i) {
                *(puc+i) = *(pvc+i) = *(pac+i) = 0;
            }
            
            // Loop through each time step to solve for displacements
            for (k = 0; k < NTSTPS; ++k) {
                
                // Initialize Meff and Reff for new time step to zero
                for (i = 0; i < NEQ; ++i) {
                    *(pMeff+i) = 0;	*(pReff+i) = 0;
                }
                
                // Calculate effective mass matrix
                if (ANAFLAG == 4) { // FSI analysis, cannot use skyline
                    for (i = 0; i < NEQ; ++i) {
                        sum = 0;
                        for (j = 0; j < NEQ; ++j) {
                            sum += *(psm_fsi+i*NEQ+j) * (*(pum+j)*a0 + *(pvm+j)*a2 + *(pam+j)*a3);
                        }
                        *(pMeff+i) = sum;
                    }
                }
                else if (ANAFLAG != 4) { // Non-FSI analysis
                    if (SLVFLAG == 0) { // using skyline function
                        for (i = 0; i < NEQ; ++i) {
                            *(pMeff+i) = *(psm+i) * (*(pum+i)*a0 + *(pvm+i)*a2 + *(pam+i)*a3);
                        }
                    }
                    else if (SLVFLAG == 1) { // using CLAPACK solver
                        for (i = 0; i < NEQ; ++i) {
                            sum = 0;
                            for (j = 0; j < NEQ; ++j) {
                                sum += *(psm+i*NEQ+j) * (*(pum+j)*a0 + *(pvm+j)*a2 + *(pam+j)*a3);
                            }
                            *(pMeff+i) = sum;
                        }
                    }
                }
                
                // Calculate effective load vector
                for (i = 0; i < NEQ; ++i) {
                    *(pReff+i) = *(pqdyn+i*NTSTPS+k) + *(pMeff+i);
                }
                
                // Calculate effective damping vector (re-use Meff)
                if (ANAFLAG == 4) {
                    for (i = 0; i < NEQ; ++i) {
                        *(pMeff+i) = *(psd_fsi+i) * (*(pum+i)*a1 + *(pvm+i)*a4 + *(pam+i)*a5);
                    }
                }
                
                // Add effective damping vector to effective load vector
                if (ANAFLAG == 4) {
                    for (i = 0; i < NEQ; ++i) {
                        *(pReff+i) += *(pMeff+i);
                    }
                }
                
                // Solve for displacements at current time step
                if (SLVFLAG == 0) {
                    err = skysolve (pmaxa, pKeff, pssd, pReff, fact, pdet);
                }
                else if (SLVFLAG == 1) {
                    dgetrs_(&trans, &n, &nrhs, pKeff, &lda, pipiv, pReff, &ldb, &info);
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
            }
        } else if (ALGFLAG == 5){ //Dynamic: nonlinear Newmark Intergration Method
            
            
            // Initialize Reff for new time step to zero
            for (i = 0; i < NEQ; ++i) {
                *(pReff+i) = 0;
            }
            
            //Compute the equivalent change of dynamic external force vector
            
            if (ANAFLAG == 4) { // FSI analysis, cannot use skyline
                for (i = 0; i<NEQ; ++i){
                    for (j = 0; j<NEQ; ++j){
                        *(pReff+i) = *(pr+i) + (a2 * (*(pvm+i))+ (*(pam+i))) * (*(psm_fsi+i*NEQ+j));
                    }
                }
            }
            else if (ANAFLAG != 4) { // Non-FSI analysis
                if (SLVFLAG == 0) { // using skyline funciton
                    for (i = 0; i < NEQ; ++i){
                        *(pReff+i) = *(pr+i) + (a2 * (*(pvm+i))+ (*(pam+i))) * (*(psm+i));
                    }
                }
                else if (SLVFLAG == 1) { //using CLAPACK solver
                    for (i = 0; i<NEQ; ++i){
                        for(j = 0; j<NEQ; ++j){
                            *(pReff+i) = *(pr+i) + (a2 * (*(pvm+i))+ (*(pam+i)))* (*(psm+i*NEQ+j));
                        }
                    }
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
            //    printf("n = %ld\t",n);
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


