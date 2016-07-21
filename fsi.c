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

#define phitol 1e-4 // Allowable +/- deviation from 1.0 of phi

extern long NJ, NE_TR, NE_FR, NE_SH, NE_SBR, NE_FBR, NEQ, SNDOF, FNDOF, NTSTPS, ntstpsinpt;
extern double dt, ttot;
extern int ANAFLAG, ALGFLAG, OPTFLAG, SLVFLAG, brFSI_FLAG, shFSI_FLAG;
extern FILE *IFP[2], *OFP[7];

void prop_fsi (double *px, double *pemod, double *pnu, double *pdens, double *pfdens, double *pbmod,
			  double *pfarea, 
			  double *pslength, double *pyield, long *pminc, long *pelface, long *pfsiinc, 
			  double *pnnorm, double *ptarea, double *pss, double *pss_fsi, double *psd_fsi, double *pabspt,
			  double *pnorpt, long *pmcode, long *pjcode, double *pL)
{
	
    // Initialize function variables
    long i, j, k, l, jt, ptr, ptr2, ptr3, NE_BR = NE_SBR + NE_FBR;
    
    int nnpfsif; // Num nodes per FSI face; = 3 if solids are shells; = 4 if solid are bricks
    int nnps; // Num nodes per solid
    int nfps; // Num faces per solid
    int nsolids;
    int ndofpe; // Number of dofs per solid 
    int ndofspn; // Number of dofs per node
    long dof;
	double Afact; // Factor to divide face area to get trib area for node
    
    if (shFSI_FLAG == 1) {
        nnpfsif = 3;
        nnps = 3;
        nfps = 1;
        nsolids = NE_SH;
        ndofpe = 18;
        ndofspn = 6;
		Afact = .333333;
    }
    
    if (brFSI_FLAG == 1) {
        nnpfsif = 4;
        nnps = 8;
        nfps = 6;
        nsolids = NE_SBR;
        ndofpe = 24;
        ndofspn = 3;
		Afact = 0.25;
    }
    
    NE_BR = NE_FBR+nsolids;
    
	long nfaces; // number of f-s faces each solid brick has
	
	double coords[4][3];
	double vec[4][3];
	double norm[4][3];
	
	double normals[nsolids*nfps*nnpfsif*3]; // local f-s face normals
	int nfpj[NJ]; // array for counting number of f-s faces per joint
	
	// Variables used to determine orientation of normal vectors
	double dotpos[4];
	
	// variables used to calculate global normal vectors
	int tot;
	double x, y, z, length;
	
	// Variables used to determine face areas
	double narea[3], farea[nsolids*nfps*nnpfsif];
	int fctyp[NJ];
	
	ptr = NE_TR + NE_FR + NE_SH;
    ptr2 = NE_TR * 2 + NE_FR * 2;
    ptr3 = NE_TR + NE_FR * 3;
        
    // Read in which joints are absorbing, if FSI analysis
    double absarea;
    if (ANAFLAG == 4) {
        for (i = 0; i < NJ; ++i) {
            fscanf(IFP[0], "%lf\n", &absarea);
            *(pabspt+i) = absarea/sqrt(*pbmod/ *pfdens);
        }
    }
        
    // Initialize system damping array to zero
    for (i = 0; i < NEQ; ++i) { *(psd_fsi+i) = 0;}
    
    // Calculate system damping array
    for (i = 0; i < NE_FBR; ++i) {
        for (j = 0; j < 8; ++j) {
            jt = *(pminc+ nsolids*nnps+i*8+j) - 1;
            if (*(pabspt+jt-1) > 0) {// absorbing point
                dof = *(pjcode+jt*7+6); //fdof
                if (dof != 0) {
                    *(psd_fsi+dof-1) = *(pabspt+jt);
                }
            }
        }
    }
	
	// Read in joint orientation points from input file
    for (i = 0; i < NJ; ++i) {
        fscanf(IFP[0], "%lf,%lf,%lf\n", pnorpt+i*3,pnorpt+i*3+1, pnorpt+i*3+2);
    }
	
	// Initialize number of f-s faces per joint to zero
	for (i = 0; i < NJ; ++i){
		*(nfpj+i) = 0;
		fctyp[i] = 0;
	}
	
	/* Determine the normal vectors and the tributary areas of the nodes on the
	 f-s interface*/
	for (i = 0; i < nsolids; ++i) {
        
		nfaces = *(pelface + i); // How many f-s faces does element i have?
		
		for (j = 0; j < nfaces; ++j) { // Loop through the number of f-s faces for current element
			
			for (l = 0; l < nnpfsif; ++l) { // Loop through num nodes per f-s face
				
				jt = *(pfsiinc + i*nfps*nnpfsif + j*nnpfsif + l);  // Global joint for current face
                
				// Get coordinates for use in evaluating local f-s face normals
				for (k = 0; k < 3; ++k) {
					coords[l][k] = *(px+(jt-1)*3+k);
				}
			}            
            
			// Evaluate vectors connecting each f-s node on the f-s face to be used to eval normals
			for (k = 0; k < 3; ++k) {
				
                if (brFSI_FLAG == 1) {
                    vec[0][k] = coords[2][k] - coords[0][k];
                    vec[1][k] = coords[3][k] - coords[2][k];
                    vec[2][k] = coords[1][k] - coords[3][k];
                    vec[3][k] = coords[0][k] - coords[1][k];
                }
                
                if (shFSI_FLAG == 1) {
                    vec[0][k] = coords[1][k] - coords[0][k];
                    vec[1][k] = coords[2][k] - coords[1][k];
                    vec[2][k] = coords[0][k] - coords[2][k];
                }                
			}
			
			// Calculate the normal vector for each node on f-s face by taking the cross product
			if (brFSI_FLAG == 1) {
                cross(&vec[0][0], &vec[3][0], &norm[0][0], 1);
                cross(&vec[3][0], &vec[2][0], &norm[1][0], 1);
                cross(&vec[1][0], &vec[0][0], &norm[2][0], 1);
                cross(&vec[2][0], &vec[1][0], &norm[3][0], 1);
            }
            
            if (shFSI_FLAG == 1) {
                cross(&vec[0][0], &vec[2][0], &norm[0][0], 1);
                cross(&vec[0][0], &vec[1][0], &norm[1][0], 1);
                cross(&vec[2][0], &vec[1][0], &norm[2][0], 1);
            }
			
			// Calculate face areas
            if (brFSI_FLAG == 1) {
                cross(&vec[0][0], &vec[3][0], narea, 0);
            }
            if (shFSI_FLAG == 1) {
                cross(&vec[0][0], &vec[1][0], narea, 0);
            }
            
			farea[i*nfps*nnpfsif+j*nnpfsif+0]=sqrt(pow(narea[0],2)+pow(narea[1],2)+pow(narea[2],2));
			farea[i*nfps*nnpfsif+j*nnpfsif+2]=farea[i*nfps*nnpfsif+j*nnpfsif+1]=farea[i*nfps*nnpfsif+j*nnpfsif+0];
            
            if (brFSI_FLAG == 1) {
                farea[i*nfps*nnpfsif+j*nnpfsif+3]=farea[i*nfps*nnpfsif+j*nnpfsif+0];
            }
			
			for (l = 0; l < nnpfsif; ++l) { // Loop through each joint on the f-s face
				
				jt = *(pfsiinc + i*nfps*nnpfsif + j*nnpfsif + l);  // Global joint for current face
				
				/* Evaluate dot products of local f-s face normal and orientation 
				 point to determine whether or not the f-s face normal is 
				 oriented correctly */
				dotpos[l] = dot((pnorpt+(jt-1)*3),&norm[l][0],3);
				
				if (dotpos[l] > -dotpos[l]) { // The normal is oriented outward from the structure
					
					for (k = 0; k < 3; ++k) {
						*(normals+i*nfps*nnpfsif*3+j*nnpfsif*3+l*3+k) = norm[l][k];
					}
				}
				else { // The normal is oriented inward to the structure
					for (k = 0; k < 3; ++k) {
						*(normals+i*nfps*nnpfsif*3+j*nnpfsif*3+l*3+k) = -norm[l][k];
					}
				}
			}			
		} 
		
		// Determine how many f-s interfaces each node is on
		for (k = 0; k < nnps; ++k) {
			
			jt = *(pminc + i*nnps + k); 
			
			for (j = 0; j < nfaces; ++j) {
				for (l = 0; l < nnpfsif; ++l) {
					
					// Increase nfpj array everytime the jt is found in fsiinc
					if (*(pfsiinc+i*nfps*nnpfsif + j*nnpfsif + l) == jt) {
						
						*(nfpj+jt-1) += 1;
					}
				}
			}
		}
	}
	
    // Calculate joint normals and tributary areas
	for (jt = 1; jt <= NJ; ++jt) {
		
		tot = 0;
		x = y = z = 0;
		length = 0;
		
		while (tot < *(nfpj+jt-1)) {
			
			for (i = 0; i < nsolids; ++i) {
				for (j = 0; j < nfps; ++j) {
					for (k = 0; k < nnpfsif; ++k) {
						
						if (*(pfsiinc+i*nfps*nnpfsif+j*nnpfsif+k) == jt) {
							
							/* Average the local f-s face normals for each joint 
							 to get global normal vectors */
							x = x + *(normals+i*nfps*nnpfsif*3+j*nnpfsif*3+k*3+0);
							y = y + *(normals+i*nfps*nnpfsif*3+j*nnpfsif*3+k*3+1);
							z = z + *(normals+i*nfps*nnpfsif*3+j*nnpfsif*3+k*3+2);
							
							tot += 1;
						}
					}
				}
			}
		}
		
		if (nfpj[jt-1] >= 1) {
			
			x/= nfpj[jt-1];
			y/= nfpj[jt-1];
			z/= nfpj[jt-1];
			
			length = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
			
			// Normalize normal vectors
			*(pnnorm+(jt-1)*3+0) = x/length;
			*(pnnorm+(jt-1)*3+1) = y/length;
			*(pnnorm+(jt-1)*3+2) = z/length;
		}
		
		else {
			*(pnnorm+(jt-1)*3+0) = 0;
			*(pnnorm+(jt-1)*3+1) = 0;
			*(pnnorm+(jt-1)*3+2) = 0;
		}
		
		/* Compute contribution to tributary area from each individual element's tributary area 
		 Check to see if the multiple faces on an f-s node are on different elements or the same
		 element*/
		for (i = 0; i < nsolids; ++i) {
			for (j = 0; j < nfps; ++j) {
				for (k = 0; k < nnpfsif; ++k) {
					
					if (*(pfsiinc+i*nfps*nnpfsif+j*nnpfsif+k) == jt) {
						
						// f-s faces are on different elements
						if (*(pnnorm+(jt-1)*3+0) == *(normals+i*nfps*nnpfsif*3+j*nnpfsif*3+k*3+0) && 
							*(pnnorm+(jt-1)*3+1) == *(normals+i*nfps*nnpfsif*3+j*nnpfsif*3+k*3+1) &&
							*(pnnorm+(jt-1)*3+2) == *(normals+i*nfps*nnpfsif*3+j*nnpfsif*3+k*3+2)) {
							
							*(ptarea+jt-1) += Afact*farea[i*nfps*nnpfsif+j*nnpfsif];
							fctyp[jt] = 1;
						}
						
						// f-s faces are on the same element
						else {
							*(ptarea+jt-1) += pow((Afact*farea[i*nfps*nnpfsif+j*nnpfsif]),2);
							fctyp[jt] = 0;
						}
					}
				}
			}
		}
        
		
		if (fctyp[jt] == 0) {
			*(ptarea+jt-1) = sqrt(*(ptarea+jt-1));
		}
	}
}

void stiff_fsi (long *pminc, long *pmcode, long *pjcode, double *pnnorm, double *ptarea, double *pfarea, double *pthick,
				double *pdeffarea, double *pslength, double *pdefslen, double *pL, double *pA, double *pss, double *pss_fsi, 
				double *px, double *pxlocal, double *pemod, double *pnu, double *pJinv, double *pjac, double *pyield, 
				double *pc1, double *pc2, double *pc3, double *pef, double *pd, double *pchi, double *pefN, double *pefM, long *pmaxa)
{

	long i, j;
	
	/* Pass control to the stiff_br function to build the partitioned
	 matrix within the system stiffness matrix*/

    stiff_br (pss, px, pemod, pnu, pminc, pmcode, pjcode, pJinv, pjac);	

	if (shFSI_FLAG == 1) {
		stiff_sh (pss, pemod, pnu, px, pxlocal, pthick, pfarea, pdeffarea, pslength,
				  pdefslen, pyield, pc1, pc2, pc3, pef, pd, pchi, pefN, pefM, pmaxa, pminc, pmcode);
	}
	
	// Initialize the system "stiffness" matrix to zero
	for (i = 0; i < NEQ; ++i) {
		for (j = 0; j < NEQ; ++j) {
			*(pss_fsi+i*NEQ+j) = 0;
		}
	}
	
	// Assemble K 
	for (i = 0; i < SNDOF; ++i) {
		for (j = 0; j < SNDOF; ++j) {
			if (*(pss+i*NEQ+j) != 0) {
				*(pss_fsi+i*NEQ+j) = *(pss+i*NEQ+j);
			}
		}
	}
	
	// Assemble H
	for (i = SNDOF; i < NEQ; ++i) {
		for (j = SNDOF; j < NEQ; ++j) {
			if (*(pss+i*NEQ+j) != 0) {
				*(pss_fsi+i*NEQ+j) = *(pss+i*NEQ+j);
			}
		}
	}
	
	// Assemble L
	for (i = 0; i < SNDOF; ++i) {
		for (j = SNDOF; j < NEQ; ++j) {
			if (*(pL+i*FNDOF+j-SNDOF) != 0) {
				*(pss_fsi+i*NEQ+j) = *(pL+i*FNDOF+j-SNDOF);
			}
		}
	}

}

void mass_fsi (long *pminc, long *pmcode, long *pjcode, double *pnnorm, double *ptarea, double *pcarea, double *pfarea, 
               double *pthick, double *pslength, double *pL, double *pLT, double *psm, double *psm_fsi, double *px, 
               double *pdens, double *pfdens, double *pJinv, double *pjac)
{
	
	long i, j;
    
	for (i = 0; i < FNDOF; ++i) {
		for (j = 0; j < SNDOF; ++j) {
			*(pLT+i*SNDOF+j) = 0;
		}
	}
	
	mass_br (psm, pdens, px, pminc, pmcode, pjac);	
    
    if (shFSI_FLAG == 1) {
		mass_sh (psm, pcarea, pdens, pthick, pfarea, pslength, px, pminc, pmcode, pjac);
	}
    
	// Initialize the system "mass" matrix to zero
	for (i = 0; i < NEQ; ++i) {
		for (j = 0; j < NEQ; ++j) {
			*(psm_fsi+i*NEQ+j) = 0;
		}
	}
	
	// Assemble M 
	for (i = 0; i < SNDOF; ++i) {
		for (j = 0; j < SNDOF; ++j) {
			if (*(psm+i*NEQ+j) != 0) {
				*(psm_fsi+i*NEQ+j) = *(psm+i*NEQ+j);
			}
		}
	}
	
	// Assemble Q
	for (i = SNDOF; i < NEQ; ++i) {
		for (j = SNDOF; j < NEQ; ++j) {
			if (*(psm+i*NEQ+j) != 0) {
				*(psm_fsi+i*NEQ+j) = *(psm+i*NEQ+j);
			}
		}
	}
	
	// Transpose(L)
	// Assemble L
	for (i = 0; i < SNDOF; ++i) {
		for (j = 0; j < FNDOF; ++j) {
			*(pLT+j*SNDOF+i) = *(pL+i*FNDOF+j);
		}
	}
	
	// Assemble lower left corner of system "mass" matrix
	for (i = SNDOF; i < NEQ; ++i) {
		for (j = 0; j < SNDOF; ++j) {
			*(psm_fsi+i*NEQ+j) = -1*(*pfdens) * (*(pLT+(i-SNDOF)*SNDOF+j));
		}
	}	
}

void L_br (long *pminc, long *pmcode, long *pjcode, long *pjcode_fsi, double *pnnorm, double *ptarea,
		   double *pL, double *pA, double *pG)
{
	long NE_BR = NE_SBR + NE_FBR;
	long i, j, m, jt;
    
    int nnpfsif; // Num nodes per FSI face; = 3 if solids are shells; = 4 if solid are bricks
    int nnps; // Num nodes per solid
    int nfps; // Num faces per solid
    int nsolids;
    int ndofpe; // Number of dofs per solid 
    int ndofspn; // Number of dofs per node
    
    if (shFSI_FLAG == 1) {
        nnpfsif = 3;
        nnps = 3;
        nfps = 1;
        nsolids = NE_SH;
        ndofpe = 18;
        ndofspn = 6;
    }
    
    if (brFSI_FLAG == 1) {
        nnpfsif = 4;
        nnps = 8;
        nfps = 6;
        nsolids = NE_SBR;
        ndofpe = 24;
        ndofspn = 3;
    }
    
    NE_BR = NE_FBR+nsolids;
	
	// Variables for building G, A, and L matrices
	int fdof, sdof[3];
	int found = 0;
	double sum = 0.0;
	
	// Initialize G and L matrices
	for (i = 0; i < SNDOF; ++i) {
		for (j = 0; j < FNDOF; ++j) {
     		*(pG+i*FNDOF+j) = 0;
			*(pL+i*FNDOF+j) = 0;
		}
	}
	
	// Initialize A matrix
	for (i = 0; i < FNDOF; ++i) {
		for (j = 0; j < FNDOF; ++j) {
			*(pA+i*FNDOF+j) = 0;
		}
	}

	/* Populate G and A matrices by looping through solid elements to find
	 f-s faces*/
	for (i = 0; i < nsolids; ++i) { // Solid elements
		
		for (j = 0; j < nnps; ++j) { //Solid element dofs
			
			jt = *(pminc+i*nnps+j) - 1; // Gloal joint
			
			// Check whether the joint has a pressure DOF
            if (*(pjcode+jt*7+6) != 0) { 
                
				fdof = *(pjcode+jt*7+6) - SNDOF; // Fluid DOF
				
				*(pA+(fdof-1)*FNDOF+fdof-1) = *(ptarea+jt);
				for (m = 0; m < 3; ++m) {
					
                    sdof[m] = *(pjcode+jt*7+m);
					*(pG+(sdof[m]-1)*FNDOF+fdof-1) = *(pnnorm+jt*3+m);
				}
			}
		}
	}

	for (i = 0; i < FNDOF; ++i) {
		if (*(pA+i*FNDOF+i) != 0) {
			for (j = 0; j < SNDOF; ++j) {
                sum = *(pG+j*FNDOF+i) * (*(pA+i*FNDOF+i));
				*(pL+j*FNDOF+i) = sum;
			}
		}
	}
}


void load_fsi (long *pjcode, double *ptinpt, double *ppinpt, double *ppresinpt, double *paccinpt, double *pfdens,
			   double *pum, double *pvm, double *pam)

{
	
	long i, j, k, kps, kpf, ks, kf, jt;
	int dir;
	double load, p, a;
	
	ntstpsinpt = ttot/dt + 1;
	
	// Initialize input loads, and fluid pressures and accelerations to zero
	for (i = 0; i < SNDOF; ++i) {
		for (j = 0; j < ntstpsinpt; ++j) {
			*(ppinpt+i*ntstpsinpt+j) = 0;
		}
	}
	
	for (i = 0; i < FNDOF; ++i) {
		for (j = 0; j < ntstpsinpt; ++j) {
			*(ppresinpt+i*ntstpsinpt+j) = 0;
			*(paccinpt+i*ntstpsinpt+j) = 0;
		}
	}

	// Initialize first element in time array to be zero
	*(ptinpt) = 0;
	
	// Assemble time array 
	for (i = 1; i < ntstpsinpt; ++i) {
		*(ptinpt+i) = *(ptinpt+i-1) + dt;
	}	
	
	/* Scan and load user input applied forces, fluid pressures and fluid accelerations.
	 Store only the applied forces corresponding to active solid DOFs and store only fluid
	 pressures and accelerations corresponding to active fluid DOFs */
	kps = 0;
	kpf = 0;
	i = 0;
	fscanf(IFP[0], "%ld,%d,%lf,%lf,%lf\n", &jt, &dir, &load, &p, &a);
    if (jt != 0) { // Check for joint loading
		
        fprintf(OFP[0], "\nJoint Loads:\n\tGlobal Joint\tDirection\tForce\t\tFInc Pressure\tNodal Acc\n");
        do {
			
            fprintf(OFP[0], "\t%ld\t\t%d\t\t%lf\t%lf\t%lf\n", jt, dir, load,p,a);
			
			// Use jcode_fsi if FSI analysis (because DOFs are renumbered)
			if (dir == 4) {
				ks = 0;
			}
			else {
                ks = *(pjcode+(jt-1)*7+dir-1); // Scan and load solid jcode
			}
			kf = *(pjcode+(jt-1)*7+6);
            
			if (ks != kps || kf != kpf) {
				i = 0;
			} 
            // Store only joint loads corresponding to active solid DOFs
            switch (ks) {
                case (0):
                    break; // Do not store loads at supports
                default:
					*(ppinpt+(ks-1)*ntstpsinpt+i) = load;
                    break;
            }
			kps = ks; 
			
			if (kf != kpf) {
				//i = 0;
			}
            
			switch (kf) {
                case (0):
                    break; // Do not store loads at supports
                default:
					*(ppresinpt+(kf-SNDOF-1)*ntstpsinpt+i) = p;
					*(paccinpt+(kf-SNDOF-1)*ntstpsinpt+i) = a;
                    break;
            }
			kpf = kf;
			
			++i;
			
            fscanf(IFP[0], "%ld,%d,%lf,%lf,%lf\n", &jt, &dir, &load, &p, &a);
            
        } while (jt != 0); // Check for last joint load
	}
    
	double d, v;
	
	/* Scan and load initial displ, vel and acc.  If dir = 4 u, v, and a refer to pressure dof
	 the its derivatives */
	fscanf(IFP[0], "%ld,%d,%lf,%lf,%lf\n", &jt, &dir, &d, &v, &a);
	fprintf(OFP[0], "\nInitial Displ, Vel and Acc:\n\tGlobal Joint\tDirection\tDisplacement\tVelocity\tAcceleration\n");
	if (jt != 0) { // Check for joint loading
		do {
			
			fprintf(OFP[0], "\t%ld\t\t%d\t\t%lf\t%lf\t%lf\n", jt, dir, d, v, a);
			
			// Use jcode_fsi if FSI analysis (because DOFs are renumbered)
            k = *(pjcode+(jt-1)*7+dir-1); // Scan and load jcode
			
			// Store only joint loads corresponding to active global DOFs
			switch (k) {
				case (0):
					break; // Do not store loads at supports
				default:
					//Are the if's necessary?
					if (d != 0) {
						*(pum+k-1) = d;
					}
					if (v != 0) {
						*(pvm+k-1) = v;
					}
					if (a != 0) {
						*(pam+k-1) = a;
					}
					break;
			}
			fscanf(IFP[0], "%ld,%d,%lf,%lf,%lf\n", &jt, &dir, &d, &v, &a);
		} while (jt != 0); // Check for last joint load		
	}
	
	/* Evaluate expression for actual dt. If actual dt < input dt, then linearlly
	 interpolate between the input loads to get load, pressure and fluid acceleration
	 values at each dt */
	double dtmax;
	dtmax = 1*dt;
	
	if (dtmax < dt) {
		dt = dtmax;
	}
	
	NTSTPS = ttot/dt + 1;
    
}


void q_fsi (long *pjcode, double *pqdyn, double *ptstps, double *papload, double *ppres, double *pacc, 
			double *pL, double *pA, double *pLp, double *pAu, double *pfdens, double *ptinpt, double *ppinpt, double *ppresinpt, double *paccinpt)

{
	// Initialize function variables
	long i, j, k;
	double sum;
	double load, p, a;
	
	double slopef, slopep, slopea; // Slopes for linear interpolation
	
	// Assemble time array based on actual time step
	*(ptstps) = *(ptinpt);
	for (i = 1; i < NTSTPS-1; ++i) {
		*(ptstps+i) = *(ptstps+i-1) + dt;
	}
	
	// Assign the last time value equal to the last input time value
	*(ptstps+NTSTPS-1) = *(ptinpt+ntstpsinpt-1);
	
	// Initialize applied load array
	for (i = 0; i < SNDOF; ++i) {
		for (j = 0; j < NTSTPS; ++j) {
			*(papload+i*NTSTPS+j) = 0;
		}
	}
	
	// Initialize fluid pressure and acceleration arrays
	for (i = 0; i < FNDOF; ++i) {
		for (j = 0; j < NTSTPS; ++j) {
			*(ppres+i*NTSTPS+j) = 0;
			*(pacc+i*NTSTPS+j) = 0;
		}
	}
	
	k = 0;
	i = 0;
	do {
		// Interpolate applied loads
		for (j = 0; j < SNDOF; ++j) {			
			slopef = (*(ppinpt+j*ntstpsinpt+k+1) - *(ppinpt+j*ntstpsinpt+k)) / (*(ptinpt+k+1) - *(ptinpt+k));			
			*(papload+j*NTSTPS+i) = slopef*(*(ptstps+i) - *(ptinpt+k)) + *(ppinpt+j*ntstpsinpt+k);
		}
		
		// Interpolate fluid  pressure and fluid incident accelerations
		for (j = 0; j < FNDOF; ++j) {
			
			slopep = (*(ppresinpt+j*ntstpsinpt+k+1) - *(ppresinpt+j*ntstpsinpt+k)) / (*(ptinpt+k+1) - *(ptinpt+k));
			slopea = (*(paccinpt+j*ntstpsinpt+k+1) - *(paccinpt+j*ntstpsinpt+k)) / (*(ptinpt+k+1) - *(ptinpt+k));
			
			*(ppres+j*NTSTPS+i) = slopep*(*(ptstps+i)-*(ptinpt+k)) + *(ppresinpt+j*ntstpsinpt+k);
			*(pacc+j*NTSTPS+i) = slopea*(*(ptstps+i)-*(ptinpt+k)) + *(paccinpt+j*ntstpsinpt+k);
		}
		++i;
		
		if (*(ptstps+i) >= *(ptinpt+k+1)) {
			++k;
		}
	}while(i < NTSTPS-1);
	
	// Assign the last load value to be equal to the last input load value
	for (i = 0; i < SNDOF; ++i) {
		*(papload+i*(NTSTPS)+NTSTPS-1) = *(ppinpt+i*(ntstpsinpt)+(ntstpsinpt-1));
	}
	
	// Assign the last fluid pressure and acceleration values to be equal to the last input values
	for (i = 0; i < FNDOF; ++i) {
		*(ppres+i*NTSTPS+NTSTPS-1) = *(ppresinpt+i*ntstpsinpt+ntstpsinpt-1);
		*(pacc+i*NTSTPS+NTSTPS-1) = *(paccinpt+i*ntstpsinpt+ntstpsinpt-1);
	}	
    
	
	// Initialize Lp matrix to zero
	for (i = 0; i < SNDOF; ++i) {
		for (j = 0; j < NTSTPS; ++j) {
			*(pLp+i*NTSTPS+j) = 0;
		}
	}
	
	// Initialize Au matrix to zero
	for (i = 0; i < FNDOF; ++i) {
		for (j = 0; j < NTSTPS; ++j) {
			*(pAu+i*NTSTPS+j) = 0;
		}
	}
	
	// Initialize the dynamic applied load array
	for (i = 0; i < NEQ; ++i) {
		for (j = 0; j < NTSTPS; ++j) {
			*(pqdyn+i*NTSTPS+j) = 0;
		}
	}
	
	// Evaluate L=L*p
	for (i = 0; i < NTSTPS; ++i) {
		for (j = 0; j < SNDOF; ++j) {
			sum = 0;
			for (k = 0; k < FNDOF; ++k) {
				sum += *(pL+j*FNDOF+k) * (*(ppres+k*NTSTPS+i));
			}
			*(pLp+j*NTSTPS+i) = sum;
		}
	}
	
	// Evaluate Au=A*u
	for (i = 0; i < NTSTPS; ++i) {
		for (j = 0; j < FNDOF; ++j) {
			sum = 0;
			for (k = 0; k < FNDOF; ++k) {
				sum += *(pA+j*FNDOF+k) * (*(pacc+k*NTSTPS+i));
			}
			*(pAu+j*NTSTPS+i) = sum;
		}
	}	
	
	// Calculate the solid DOFs of the load array
	for (i = 0; i < SNDOF; ++i) {
		for (j = 0; j < NTSTPS; ++j) {
			*(pqdyn+i*NTSTPS+j) = *(papload+i*NTSTPS+j) - *(pLp+i*NTSTPS+j);
		}
	}
	
	// Evaluate the fluid DOFs of the load array
	for (i = SNDOF; i < NEQ; ++i) {
		for (j = 0; j < NTSTPS; ++j) {
			*(pqdyn+i*NTSTPS+j) = -1*(*(pfdens)) * (*(pAu+(i-SNDOF)*NTSTPS+j));
		}
	}    
}

