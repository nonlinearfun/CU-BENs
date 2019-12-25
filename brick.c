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

#define phitol 1e-4 // Allowable +/- deviation from 1.0 of phi

extern long NJ, NE_TR, NE_FR, NE_SH, NE_SBR, NE_FBR, NEQ, SNDOF, FNDOF;
extern int ANAFLAG, ALGFLAG, OPTFLAG, SLVFLAG, brFSI_FLAG, shFSI_FLAG;
extern FILE *IFP[4], *OFP[8];

void prop_br (double *pemod, double *pnu, double *pyield, double *pdens, double *pfdens, double *pbmod,
			   double *pfarea)
{
	// Initialize function variables
    long i, ptr, ptr2, ptr3;
	
	ptr = NE_TR + NE_FR + NE_SH;
    ptr2 = NE_TR * 2 + NE_FR * 2;
    ptr3 = NE_TR + NE_FR * 3;
	
    // Read in element properties from input file
	if (NE_SBR > 0) {		
		for (i = 0; i < NE_SBR; ++i) {			
			fscanf(IFP[0], "%lf,%lf,%lf,%lf\n", pemod+ptr+i, pnu+ptr+i, pdens+ptr+i, pyield+ptr+i);
		}
	}
	
    // Scan fluid properties from input file
	if (NE_FBR > 0) {	
		fscanf(IFP[0], "%lf,%lf\n", pfdens, pbmod);		
		double wvsp = sqrt(*pbmod/ *pfdens); // Wave speed in the fluid
		
		/* Re-assign the fluid density to be equal to Ge/c^2 
		 Ge is arbitrarily chosen as unity (Everstine, 1997)*/
		for (i = 0; i < NE_FBR; ++i) {
			*(pdens+ptr+NE_SBR+i) = 1/(pow(wvsp,2));
		}		
		for (i = 0; i < NE_FBR; ++i) {			
			// Assign values for fluid properties
			*(pemod+ptr+NE_SBR+i) = 1*pow(10,20);
			*(pnu+ptr+NE_SBR+i) = .5*pow(10,20);
		}
	}
}


void stiff_br (double *pss, double *px, double *pemod, double *pnu, long *pminc, 
			   long *pmcode, long *pjcode, double *pJinv, double *pjac)
{
	long NE_BR = NE_SBR + NE_FBR;
	long ptr = NE_TR + NE_FR + NE_SH;
	long ptr3 = NE_TR * 6 + NE_FR * 14 + NE_SH * 18;
	long i, j, k, l, n, ie, je;
    int  r, s, t;
	double sum = 0;
	
	// Set up flags to determine whether the solid elements are shells or bricks
	if (brFSI_FLAG == 1) {
		NE_BR = NE_SBR + NE_FBR;
		ptr = NE_TR + NE_FR + NE_SH;
	}
	if (shFSI_FLAG == 1) {
		NE_BR =  NE_FBR;
		ptr = NE_TR + NE_FR + NE_SH;
	}
	
	double k_br[24][24]; // General element stiffness matrix in global coordinate system
	double temp[6][24], temp2[24][24]; // Temp matrices used for evaluating element stiffness
	
	double B[6][24]; // Strain-displacement matrix
	
	double C[6][6];  // Material matrix
	double e1, e2, e3;
	
	// Gauss-quadrature integration points
	double R[2], S[2], T[2];
	R[0] = S[0] = T[0] =  1.0/sqrt(3);
	R[1] = S[1] = T[1] = -1.0/sqrt(3);
	
	// Derivatives of shape functions w.r.t. global coords
	double dhdx[8], dhdy[8], dhdz[8];
	
	// Jacobian, and determinant of Jacobian
	double detJ;
	for (i = 0; i < NE_BR; ++i) {
		
		// Initialize element stiffness array and temp2 to zero
		for (j = 0; j < 24; ++j) {
			for (k = 0; k < 24; k++) {
				k_br[j][k] = 0;
			}
		}
		
		// Calculate elements of material matrix C
		e1 = *(pemod+ptr+i)*(*(pnu+ptr+i))/((1+*(pnu+ptr+i))*(1-2*(*(pnu+ptr+i)))); // Ev/((1+v)(1-2v))
		e2 = .5*(*(pemod+ptr+i)/((1+*(pnu+ptr+i)))); // E/2(1+v)
		e3 = *(pemod+ptr+i)*(1-*(pnu+ptr+i))/((1+*(pnu+ptr+i))*(1-2*(*(pnu+ptr+i)))); // E(1-v)/((1+v)(1-2v))
        
		// Initialize C to zero
		for (j = 0; j < 6; ++j) {
			for (k = 0; k < 6; k++) {
				C[j][k]= 0;
			}
		}
		
		// Evaluate material matrix C
		C[0][0] = C[1][1] = C[2][2] = e3;
		C[0][1] = C[0][2] = C[1][0] = C[1][2] = C[2][0] = C[2][1] = e1;
		C[3][3] = C[4][4] = C[5][5] = e2;
		
		// Loop over integration points in all three dimensions
		for (r = 0; r < 2; ++r) {
			for (s = 0; s < 2; ++s) {
				for (t = 0; t < 2; ++t) {
					
					// Initialize J to zero for each integration point
					for (j = 0; j < 3; ++j) {
						for (k = 0; k < 3; ++k) {
							*(pJinv+j*3+k) = 0;
						}
					}
                    
					// Initialize B and temp to zero for each integration point
					for (j = 0; j < 6; ++j) {
						for (k = 0; k < 24; k++) {
							B[j][k] = 0;
							temp[j][k]= 0;
						}
					}
					// Initialize temp2 to zero
					for (j = 0; j < 24; ++j) {
						for (k = 0; k < 24; k++) {
							temp2[j][k]= 0;
						}
					}
					
					// Initialize J to zero for each integration point
					for (j = 0; j < 3; ++j) {
						for (k = 0; k < 3; ++k) {
							*(pjac+j*3+k) = 0;
						}
					}
					
					/* Pass control to jacob in order to calculate the Jacobian for the current
					 integration point*/
					jacob (px, pminc, &i, &r, &s, &t, pjac);
					
					// Initialize Jinv array with Jacobian values 
					for (j = 0; j < 3; ++j) {
						for (k = 0; k < 3; ++k) {
							*(pJinv+j*3+k) = *(pjac+j*3+k);
						}
					}
					
					// Calculate the inverse of the Jacobian
					inverse (pJinv, 3);
					
					// Evaluate det(Jacobian)
					detJ = *(pjac+0*3+0)*(*(pjac+1*3+1))*(*(pjac+2*3+2)) - *(pjac+0*3+0)*(*(pjac+1*3+2))*(*(pjac+2*3+1)) - 
					*(pjac+0*3+1)*(*(pjac+1*3+0))*(*(pjac+2*3+2)) + *(pjac+0*3+1)*(*(pjac+1*3+2))*(*(pjac+2*3+0)) + 
					*(pjac+0*3+2)*(*(pjac+1*3+0))*(*(pjac+2*3+1)) - *(pjac+0*3+2)*(*(pjac+1*3+1))*(*(pjac+2*3+0));
                    
					// Loop through each node of the current element
					for (n = 0; n < 8; ++n) {
						
						// Calculate the derivative of the shape functions w.r.t. global coords
						switch (n) {
							case 0:
								
								dhdx[0] = *(pJinv+0*3+0)*((S[s] + 1.0)*(T[t] + 1.0))/8.0 + 
								*(pJinv+0*3+1)*(R[r]/8.0 + 1.0/8.0)*(T[t] + 1.0) +
								*(pJinv+0*3+2)*(R[r]/8.0 + 1.0/8.0)*(S[s] + 1.0);
								
								dhdy[0] = *(pJinv+1*3+0)*((S[s] + 1.0)*(T[t] + 1.0))/8.0 + 
								*(pJinv+1*3+1)*(R[r]/8.0 + 1.0/8.0)*(T[t] + 1.0) +
								*(pJinv+1*3+2)*(R[r]/8.0 + 1.0/8.0)*(S[s] + 1.0);
								
								dhdz[0] = *(pJinv+2*3+0)*((S[s] + 1.0)*(T[t] + 1.0))/8.0 + 
								*(pJinv+2*3+1)*(R[r]/8.0 + 1.0/8.0)*(T[t] + 1.0) +
								*(pJinv+2*3+2)*(R[r]/8.0 + 1.0/8.0)*(S[s] + 1.0);
								
								break;
							case 1:
								
								dhdx[1] = *(pJinv+0*3+0)*-((S[s] + 1.0)*(T[t] + 1.0))/8.0 + 
								*(pJinv+0*3+1)*-(R[r]/8.0 - 1.0/8.0)*(T[t] + 1.0) +
								*(pJinv+0*3+2)*-(R[r]/8.0 - 1.0/8.0)*(S[s] + 1.0);
								
								dhdy[1] = *(pJinv+1*3+0)*-((S[s] + 1.0)*(T[t] + 1.0))/8.0 + 
								*(pJinv+1*3+1)*-(R[r]/8.0 - 1.0/8.0)*(T[t] + 1.0) +
								*(pJinv+1*3+2)*-(R[r]/8.0 - 1.0/8.0)*(S[s] + 1.0);
								
								dhdz[1] = *(pJinv+2*3+0)*-((S[s] + 1.0)*(T[t] + 1.0))/8.0 + 
								*(pJinv+2*3+1)*-(R[r]/8.0 - 1.0/8.0)*(T[t] + 1.0) +
								*(pJinv+2*3+2)*-(R[r]/8.0 - 1.0/8.0)*(S[s] + 1.0);
								
								break;	
							case 2:
								
								dhdx[2] = *(pJinv+0*3+0)*((S[s] -1.0)*(T[t] +1.0))/8.0 + 
								*(pJinv+0*3+1)*(R[r]/8.0 -1.0/8.0)*(T[t] +1.0) +
								*(pJinv+0*3+2)*(R[r]/8.0 -1.0/8.0)*(S[s] -1.0);
								
								dhdy[2] = *(pJinv+1*3+0)*((S[s] -1.0)*(T[t] +1.0))/8.0 + 
								*(pJinv+1*3+1)*(R[r]/8.0 -1.0/8.0)*(T[t] +1.0) +
								*(pJinv+1*3+2)*(R[r]/8.0 -1.0/8.0)*(S[s] -1.0);
								
								dhdz[2] = *(pJinv+2*3+0)*((S[s] -1.0)*(T[t] +1.0))/8.0 + 
								*(pJinv+2*3+1)*(R[r]/8.0 -1.0/8.0)*(T[t] +1.0) +
								*(pJinv+2*3+2)*(R[r]/8.0 -1.0/8.0)*(S[s] -1.0);
								
								break;
							case 3:
								dhdx[3] = *(pJinv+0*3+0)*-((S[s] -1.0)*(T[t] +1.0))/8.0 + 
								*(pJinv+0*3+1)*-(R[r]/8.0 +1.0/8.0)*(T[t] +1.0) +
								*(pJinv+0*3+2)*-(R[r]/8.0 +1.0/8.0)*(S[s] -1.0);
								
								dhdy[3] = *(pJinv+1*3+0)*-((S[s] -1.0)*(T[t] +1.0))/8.0 + 
								*(pJinv+1*3+1)*-(R[r]/8.0 +1.0/8.0)*(T[t] +1.0) +
								*(pJinv+1*3+2)*-(R[r]/8.0 +1.0/8.0)*(S[s] -1.0);
								
								dhdz[3] = *(pJinv+2*3+0)*-((S[s] -1.0)*(T[t] +1.0))/8.0 + 
								*(pJinv+2*3+1)*-(R[r]/8.0 +1.0/8.0)*(T[t] +1.0) +
								*(pJinv+2*3+2)*-(R[r]/8.0 +1.0/8.0)*(S[s] -1.0);
								
								break;
							case 4:
								dhdx[4] = *(pJinv+0*3+0)*-((S[s] +1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+0*3+1)*-(R[r]/8.0 +1.0/8.0)*(T[t] -1.0) +
								*(pJinv+0*3+2)*-(R[r]/8.0 +1.0/8.0)*(S[s] +1.0);
								
								dhdy[4] = *(pJinv+1*3+0)*-((S[s] +1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+1*3+1)*-(R[r]/8.0 +1.0/8.0)*(T[t] -1.0) +
								*(pJinv+1*3+2)*-(R[r]/8.0 +1.0/8.0)*(S[s] +1.0);
								
								dhdz[4] = *(pJinv+2*3+0)*-((S[s] +1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+2*3+1)*-(R[r]/8.0 +1.0/8.0)*(T[t] -1.0) +
								*(pJinv+2*3+2)*-(R[r]/8.0 +1.0/8.0)*(S[s] +1.0);
								
								break;
							case 5:
								dhdx[5] = *(pJinv+0*3+0)*((S[s] +1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+0*3+1)*(R[r]/8.0 -1.0/8.0)*(T[t] -1.0) +
								*(pJinv+0*3+2)*(R[r]/8.0 -1.0/8.0)*(S[s] +1.0);
								
								dhdy[5] = *(pJinv+1*3+0)*((S[s] +1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+1*3+1)*(R[r]/8.0 -1.0/8.0)*(T[t] -1.0) +
								*(pJinv+1*3+2)*(R[r]/8.0 -1.0/8.0)*(S[s] +1.0);
								
								dhdz[5] = *(pJinv+2*3+0)*((S[s] +1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+2*3+1)*(R[r]/8.0 -1.0/8.0)*(T[t] -1.0) +
								*(pJinv+2*3+2)*(R[r]/8.0 -1.0/8.0)*(S[s] +1.0);
								
								break;
							case 6:
								dhdx[6] = *(pJinv+0*3+0)*-((S[s] -1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+0*3+1)*-(R[r]/8.0 -1.0/8.0)*(T[t] -1.0) +
								*(pJinv+0*3+2)*-(R[r]/8.0 -1.0/8.0)*(S[s] -1.0);
								
								dhdy[6] = *(pJinv+1*3+0)*-((S[s] -1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+1*3+1)*-(R[r]/8.0 -1.0/8.0)*(T[t] -1.0) +
								*(pJinv+1*3+2)*-(R[r]/8.0 -1.0/8.0)*(S[s] -1.0);
								
								dhdz[6] = *(pJinv+2*3+0)*-((S[s] -1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+2*3+1)*-(R[r]/8.0 -1.0/8.0)*(T[t] -1.0) +
								*(pJinv+2*3+2)*-(R[r]/8.0 -1.0/8.0)*(S[s] -1.0);
								
								break;
							case 7:
								dhdx[7] = *(pJinv+0*3+0)*((S[s] -1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+0*3+1)*(R[r]/8.0 +1.0/8.0)*(T[t] -1.0) +
								*(pJinv+0*3+2)*(R[r]/8.0 +1.0/8.0)*(S[s] -1.0);
								
								dhdy[7] = *(pJinv+1*3+0)*((S[s] -1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+1*3+1)*(R[r]/8.0 +1.0/8.0)*(T[t] -1.0) +
								*(pJinv+1*3+2)*(R[r]/8.0 +1.0/8.0)*(S[s] -1.0);
								
								dhdz[7] = *(pJinv+2*3+0)*((S[s] -1.0)*(T[t] -1.0))/8.0 + 
								*(pJinv+2*3+1)*(R[r]/8.0 +1.0/8.0)*(T[t] -1.0) +
								*(pJinv+2*3+2)*(R[r]/8.0 +1.0/8.0)*(S[s] -1.0);
								
								break;
						}
					}
					
					for (j = 0; j < 6; ++j) {
						for (k = 0; k < 24; k++) {
							B[j][k] = 0;
							temp[j][k]= 0;
						}
					}
					
					// Build B from derivative of shape functions w.r.t. global coords						
					B[0][0] = dhdx[0]; B[0][3] = dhdx[1]; B[0][6] = dhdx[2]; B[0][9] = dhdx[3]; 
					B[0][12] = dhdx[4]; B[0][15] = dhdx[5]; B[0][18] = dhdx[6]; B[0][21] = dhdx[7];
					
					B[1][1] = dhdy[0]; B[1][4] = dhdy[1]; B[1][7] = dhdy[2]; B[1][10] = dhdy[3]; 
					B[1][13] = dhdy[4]; B[1][16] = dhdy[5]; B[1][19] = dhdy[6]; B[1][22] = dhdy[7];
					
					B[2][2] = dhdz[0]; B[2][5] = dhdz[1]; B[2][8] = dhdz[2]; B[2][11] = dhdz[3]; 
					B[2][14] = dhdz[4]; B[2][17] = dhdz[5]; B[2][20] = dhdz[6]; B[2][23] = dhdz[7];
					
					B[3][0] = dhdy[0]; B[3][1] = dhdx[0]; B[3][3] = dhdy[1]; B[3][4] = dhdx[1];
					B[3][6] = dhdy[2]; B[3][7] = dhdx[2]; B[3][9] = dhdy[3]; B[3][10] = dhdx[3];
					B[3][12] = dhdy[4]; B[3][13] = dhdx[4]; B[3][15] = dhdy[5]; B[3][16] = dhdx[5];
					B[3][18] = dhdy[6]; B[3][19] = dhdx[6]; B[3][21] = dhdy[7]; B[3][22] = dhdx[7];
					
					B[4][1] = dhdz[0]; B[4][2] = dhdy[0]; B[4][4] = dhdz[1]; B[4][5] = dhdy[1];
					B[4][7] = dhdz[2]; B[4][8] = dhdy[2]; B[4][10] = dhdz[3]; B[4][11] = dhdy[3];
					B[4][13] = dhdz[4]; B[4][14] = dhdy[4]; B[4][16] = dhdz[5]; B[4][17] = dhdy[5];
					B[4][19] = dhdz[6]; B[4][20] = dhdy[6]; B[4][22] = dhdz[7]; B[4][23] = dhdy[7];
					
					B[5][0] = dhdz[0]; B[5][2] = dhdx[0]; B[5][3] = dhdz[1]; B[5][5] = dhdx[1];
					B[5][6] = dhdz[2]; B[5][8] = dhdx[2]; B[5][9] = dhdz[3]; B[5][11] = dhdx[3];
					B[5][12] = dhdz[4]; B[5][14] = dhdx[4]; B[5][15] = dhdz[5]; B[5][17] = dhdx[5];
					B[5][18] = dhdz[6]; B[5][20] = dhdx[6]; B[5][21] = dhdz[7]; B[5][23] = dhdx[7];
					
					// Evaluate C*B
					for (l = 0; l < 24; ++l) {
						for (j = 0; j < 6; ++j) {
							sum = 0;
							for (k = 0; k < 6; ++k) {
								sum += C[j][k]*B[k][l];
							}
							temp[j][l] = sum;
						}
					}
					
					// Evaluate B'*C*B*detJ
					for (l = 0; l < 24; ++l) {
						for (j = 0; j < 24; ++j) {
							sum = 0;
							for (k = 0; k < 6; ++k) {
								sum += B[k][j]*temp[k][l];
							}
							temp2[j][l] = sum*detJ;
						}
					}
					
					/* Add numerically integrated B'DB*detJ at each integration point
					 to the element stiffness matrix*/
					for (j = 0; j < 24; ++j) {
						for (k = 0; k < 24; ++k) {
							k_br[j][k] += temp2[j][k];
						}
					}
				}
			}
		}
		
		// Assemble system stiffness matrix - full order [NEQ][NEQ]
		for (ie = 0; ie < 24; ++ie) {
			for (je = 0; je < 24; ++je) {
				
                
				j = *(pmcode+ptr3+i*24+ie);
				k = *(pmcode+ptr3+i*24+je);
                
				if ((j != 0) && (k != 0)) {
					*(pss+(j-1)*NEQ+k-1) += k_br[je][ie];
				}
			}
		}
	}	
}

void mass_br (double *psm, double *pdens, double *px, long *pminc, long *pmcode, double *pjac)
{
	long NE_BR = NE_SBR + NE_FBR;
	long ptr = NE_TR + NE_FR + NE_SH;
	long ptr3 = NE_TR * 6 + NE_FR * 14 + NE_SH * 18;
    long i, j, k, l, n, ie, je;
    int  r, s, t;
	double sum = 0;
	
	// Set up flags to determine whether the solid elements are shells or bricks	
	if (brFSI_FLAG == 1) {
		NE_BR = NE_SBR + NE_FBR;
		ptr = NE_TR + NE_FR + NE_SH;
	}
	if (shFSI_FLAG == 1) {
		NE_BR =  NE_FBR;
		ptr = NE_TR + NE_FR + NE_SH;
	}
	
	double m_br[24][24]; // General element mass matrix 
	double h[8];  // The shape functions array
	double H[3][24]; // shape functions matrix
	double HT[24][3];
	double temp[24][24]; // Temp matrix used for evaluating element stiffness
	
	// Gauss-quadrature integration points
	double R[2], S[2], T[2];
	R[0] = S[0] = T[0] =  1.0/sqrt(3);
	R[1] = S[1] = T[1] = -1.0/sqrt(3);
	
	// Determinant of Jacobian
	double detJ;
	
	for (i = 0; i < NE_BR; ++i) {
		
		// Initialize element mass array to zero
		for (j = 0; j < 24; ++j) {
			for (k = 0; k < 24; k++) {
				m_br[j][k] = 0;
			}
		}
		
		// Loop over integration points in all three dimensions
		for (r = 0; r < 2; ++r) {
			for (s = 0; s < 2; ++s) {
				for (t = 0; t < 2; ++t) {
					
					// Initialize H to zero
					for (j = 0; j < 3; ++j) {
						for (k = 0; k < 24; ++k) {
							HT[k][j] = H[j][k] = 0;
						}
					}
					
					// Initialize Jacobian matrix to zero
					for (j = 0; j < 3; ++j) {
						for (k = 0; k < 3; ++k) {
							*(pjac+j*3+k) = 0;
						}
					}
					
					/* Pass control to jacob in order to calculate the Jacobian for the current
					 integration point*/
					jacob (px, pminc, &i, &r, &s, &t, pjac);
					
					// Evaluate det(Jacobian)
					detJ = *(pjac+0*3+0)*(*(pjac+1*3+1))*(*(pjac+2*3+2)) - *(pjac+0*3+0)*(*(pjac+1*3+2))*(*(pjac+2*3+1)) - 
					*(pjac+0*3+1)*(*(pjac+1*3+0))*(*(pjac+2*3+2)) + *(pjac+0*3+1)*(*(pjac+1*3+2))*(*(pjac+2*3+0)) + 
					*(pjac+0*3+2)*(*(pjac+1*3+0))*(*(pjac+2*3+1)) - *(pjac+0*3+2)*(*(pjac+1*3+1))*(*(pjac+2*3+0));
                    

					// Evaluate the shape functions H
					h[0] = (1.0+R[r])*(1.0+S[s])*(1.0+T[t])/8.0;
					h[1] = (1.0-R[r])*(1.0+S[s])*(1.0+T[t])/8.0;
					h[2] = (1.0-R[r])*(1.0-S[s])*(1.0+T[t])/8.0;
					h[3] = (1.0+R[r])*(1.0-S[s])*(1.0+T[t])/8.0;
					h[4] = (1.0+R[r])*(1.0+S[s])*(1.0-T[t])/8.0;
					h[5] = (1.0-R[r])*(1.0+S[s])*(1.0-T[t])/8.0;
					h[6] = (1.0-R[r])*(1.0-S[s])*(1.0-T[t])/8.0;
					h[7] = (1.0+R[r])*(1.0-S[s])*(1.0-T[t])/8.0;
					
					H[0][0] = H[1][1] = H[2][2] = h[0];
					H[0][3] = H[1][4] = H[2][5] = h[1];
					H[0][6] = H[1][7] = H[2][8] = h[2];
					H[0][9] = H[1][10] = H[2][11] = h[3];
					H[0][12] = H[1][13] = H[2][14] = h[4];
					H[0][15] = H[1][16] = H[2][17] = h[5];
					H[0][18] = H[1][19] = H[2][20] = h[6];
					H[0][21] = H[1][22] = H[2][23] = h[7];
					
					// transpose(H)
					for (j = 0; j < 3; ++j) {
						for (k = 0; k < 24; ++k) {
							HT[k][j] = H[j][k];
						}
					}
					
					// Initialize temp to zero for each integration point
					for (j = 0; j < 24; ++j) {
						for (k = 0; k < 24; k++) {
							temp[j][k]= 0;
						}
					}
					
					// Evaluate dens*H^T*H*detJ
					for (l = 0; l < 24; ++l) {
						for (j = 0; j < 24; ++j) {
							sum = 0;
							for (k = 0; k < 3; ++k) {
								sum += HT[j][k]*H[k][l];
							}
							temp[j][l] = *(pdens+ptr+i)*sum*detJ;
						}
					}
					
					/* Add numerically integrated H^T*H*detJ at each integration point
					 to the element mass matrix*/
					for (j = 0; j < 24; ++j) {
						for (k = 0; k < 24; ++k) {
							m_br[j][k] += temp[j][k];
						}
					}
				}
			}
		}
		
		// Assemble system mass matrix - FULL ORDER [NEQ][NEQ]
		for (ie = 0; ie < 24; ++ie) {
			for (je = 0; je < 24; ++je) {
				j = *(pmcode+ptr3+i*24+ie);
				k = *(pmcode+ptr3+i*24+je);
				
				if ((j != 0) && (k != 0)) {
					*(psm+(j-1)*NEQ+k-1) += m_br[je][ie];
				}
			}
		}
	}
}



void jacob (double *px, long *pminc, long *el, int *rval, int *sval, int *tval, double *pjac)

{
    
    // Jacobian function
    
	long ptr2 = NE_TR * 2 + NE_FR * 2 + NE_SH * 3;
	long i, n, r, s, t, jt;
	
	i = *el;
	r = *rval;
	s = *sval;
	t = *tval;
	
	// Gauss-quadrature integration points
	double R[2], S[2], T[2];
	R[0] = S[0] = T[0] =  1.0/sqrt(3);
	R[1] = S[1] = T[1] = -1.0/sqrt(3);
	
	// Derivatives of X, Y, Z w.r.t. natural coordinates
	double dxdr[8], dydr[8], dzdr[8];
	double dxds[8], dyds[8], dzds[8];
	double dxdt[8], dydt[8], dzdt[8];
	
	for (n = 0; n < 8; ++n) {
		dxdr[n] = 0; dxds[n] = 0; dxdt[n] = 0;
		dydr[n] = 0; dyds[n] = 0; dydt[n] = 0;
		dzdr[n] = 0; dzds[n] = 0; dzdt[n] = 0;
	}
	
	// Loop over the 8 nodes per brick element
	for (n = 0; n < 8; ++n) {
		
		jt = *(pminc+ptr2+i*8+n) - 1; // Determine global joint

		// Calculate dx/dr, etc. for each integration point at each node in the element
		switch (n) {
			case 0:
				dxdr[0] = ((S[s] + 1.0)*(T[t] + 1.0))/8.0*(*(px+jt*3+0));
				dxds[0] = (R[r]/8.0 + 1.0/8.0)*(T[t] + 1.0)*(*(px+jt*3+0));
				dxdt[0] = (R[r]/8.0 + 1.0/8.0)*(S[s] + 1.0)*(*(px+jt*3+0));
				
				dydr[0] = ((S[s] + 1.0)*(T[t] + 1.0))/8.0*(*(px+jt*3+1));
				dyds[0] = (R[r]/8.0 + 1.0/8.0)*(T[t] + 1.0)*(*(px+jt*3+1));
				dydt[0] = (R[r]/8.0 + 1.0/8.0)*(S[s] + 1.0)*(*(px+jt*3+1));
				
				dzdr[0] = ((S[s] + 1.0)*(T[t] + 1.0))/8.0*(*(px+jt*3+2));
				dzds[0] = (R[r]/8.0 + 1.0/8.0)*(T[t] + 1.0)*(*(px+jt*3+2));
				dzdt[0] = (R[r]/8.0 + 1.0/8.0)*(S[s] + 1.0)*(*(px+jt*3+2));
				
				break;
			case 1:
				dxdr[1] = -((S[s] + 1.0)*(T[t] + 1.0))/8.0*(*(px+jt*3+0));
				dxds[1] = -(R[r]/8.0 - 1.0/8.0)*(T[t] + 1.0)*(*(px+jt*3+0));
				dxdt[1] = -(R[r]/8.0 - 1.0/8.0)*(S[s] + 1.0)*(*(px+jt*3+0));
				
				dydr[1] = -((S[s] + 1.0)*(T[t] + 1.0))/8.0*(*(px+jt*3+1));
				dyds[1] = -(R[r]/8.0 - 1.0/8.0)*(T[t] + 1.0)*(*(px+jt*3+1));
				dydt[1] = -(R[r]/8.0 - 1.0/8.0)*(S[s] + 1.0)*(*(px+jt*3+1));
				
				dzdr[1] = -((S[s] + 1.0)*(T[t] + 1.0))/8.0*(*(px+jt*3+2));
				dzds[1] = -(R[r]/8.0 - 1.0/8.0)*(T[t] + 1.0)*(*(px+jt*3+2));
				dzdt[1] = -(R[r]/8.0 - 1.0/8.0)*(S[s] + 1.0)*(*(px+jt*3+2));
				
				break;
			case 2:
				dxdr[2] = ((S[s] -1.0)*(T[t] +1.0))/8.0*(*(px+jt*3+0));
				dxds[2] = (R[r]/8.0 -1.0/8.0)*(T[t] +1.0)*(*(px+jt*3+0));
				dxdt[2] = (R[r]/8.0 -1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+0));
				
				dydr[2] = ((S[s] -1.0)*(T[t] +1.0))/8.0*(*(px+jt*3+1));
				dyds[2] = (R[r]/8.0 -1.0/8.0)*(T[t] +1.0)*(*(px+jt*3+1));
				dydt[2] = (R[r]/8.0 -1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+1));
				
				dzdr[2] = ((S[s] -1.0)*(T[t] +1.0))/8.0*(*(px+jt*3+2));
				dzds[2] = (R[r]/8.0 -1.0/8.0)*(T[t] +1.0)*(*(px+jt*3+2));
				dzdt[2] = (R[r]/8.0 -1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+2));
				
				
				break;
			case 3:
				dxdr[3] = -((S[s] -1.0)*(T[t] +1.0))/8.0*(*(px+jt*3+0));
				dxds[3] = -(R[r]/8.0 +1.0/8.0)*(T[t] +1.0)*(*(px+jt*3+0));
				dxdt[3] = -(R[r]/8.0 +1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+0));
				
				dydr[3] = -((S[s] -1.0)*(T[t] +1.0))/8.0*(*(px+jt*3+1));
				dyds[3] = -(R[r]/8.0 +1.0/8.0)*(T[t] +1.0)*(*(px+jt*3+1));
				dydt[3] = -(R[r]/8.0 +1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+1));
				
				dzdr[3] = -((S[s] -1.0)*(T[t] +1.0))/8.0*(*(px+jt*3+2));
				dzds[3] = -(R[r]/8.0 +1.0/8.0)*(T[t] +1.0)*(*(px+jt*3+2));
				dzdt[3] = -(R[r]/8.0 +1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+2));
				
				break;
			case 4:
				dxdr[4] = -((S[s] +1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+0));
				dxds[4] = -(R[r]/8.0 +1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+0));
				dxdt[4] = -(R[r]/8.0 +1.0/8.0)*(S[s] +1.0)*(*(px+jt*3+0));
				
				dydr[4] = -((S[s] +1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+1));
				dyds[4] = -(R[r]/8.0 +1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+1));
				dydt[4] = -(R[r]/8.0 +1.0/8.0)*(S[s] +1.0)*(*(px+jt*3+1));
				
				dzdr[4] = -((S[s] +1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+2));
				dzds[4] = -(R[r]/8.0 +1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+2));
				dzdt[4] = -(R[r]/8.0 +1.0/8.0)*(S[s] +1.0)*(*(px+jt*3+2));
				
				break;
			case 5:
				dxdr[5] = ((S[s] +1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+0));
				dxds[5] = (R[r]/8.0 -1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+0));
				dxdt[5] = (R[r]/8.0 -1.0/8.0)*(S[s] +1.0)*(*(px+jt*3+0));
				
				dydr[5] = ((S[s] +1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+1));
				dyds[5] = (R[r]/8.0 -1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+1));
				dydt[5] = (R[r]/8.0 -1.0/8.0)*(S[s] +1.0)*(*(px+jt*3+1));
				
				dzdr[5] = ((S[s] +1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+2));
				dzds[5] = (R[r]/8.0 -1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+2));
				dzdt[5] = (R[r]/8.0 -1.0/8.0)*(S[s] +1.0)*(*(px+jt*3+2));
				
				break;
			case 6:
				dxdr[6] = -((S[s] -1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+0));
				dxds[6] = -(R[r]/8.0 -1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+0));
				dxdt[6] = -(R[r]/8.0 -1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+0));
				
				dydr[6] = -((S[s] -1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+1));
				dyds[6] = -(R[r]/8.0 -1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+1));
				dydt[6] = -(R[r]/8.0 -1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+1));
				
				dzdr[6] = -((S[s] -1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+2));
				dzds[6] = -(R[r]/8.0 -1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+2));
				dzdt[6] = -(R[r]/8.0 -1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+2));
				
				break;
			case 7:
				dxdr[7] = ((S[s] -1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+0));
				dxds[7] = (R[r]/8.0 +1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+0));
				dxdt[7] = (R[r]/8.0 +1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+0));
				
				dydr[7] = ((S[s] -1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+1));
				dyds[7] = (R[r]/8.0 +1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+1));
				dydt[7] = (R[r]/8.0 +1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+1));
				
				dzdr[7] = ((S[s] -1.0)*(T[t] -1.0))/8.0*(*(px+jt*3+2));
				dzds[7] = (R[r]/8.0 +1.0/8.0)*(T[t] -1.0)*(*(px+jt*3+2));
				dzdt[7] = (R[r]/8.0 +1.0/8.0)*(S[s] -1.0)*(*(px+jt*3+2));
				
				break;
		}
		
		*(pjac+0*3+0) += dxdr[n]; *(pjac+0*3+1) += dydr[n]; *(pjac+0*3+2) += dzdr[n];
		*(pjac+1*3+0) += dxds[n]; *(pjac+1*3+1) += dyds[n]; *(pjac+1*3+2) += dzds[n];
		*(pjac+2*3+0) += dxdt[n]; *(pjac+2*3+1) += dydt[n]; *(pjac+2*3+2) += dzdt[n];
		
	}
	
}




