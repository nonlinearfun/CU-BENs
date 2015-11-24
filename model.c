
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
#include <math.h>
#include "prototypes.h"

extern long NJ, NE_TR, NE_FR, NE_SH, NE_SBR, NE_FBR, NEQ, SNDOF, FNDOF, NTSTPS, ntstpsinpt;
extern double dt, ttot;
extern int ANAFLAG, ALGFLAG, OPTFLAG, SLVFLAG, FSIFLAG, FSIINCFLAG, brFSI_FLAG, shFSI_FLAG;
extern FILE *IFP[2], *OFP[7];

int struc (long *pjcode, long *pminc, int *pwrpres, long *pjnt)
{
    long i, j, k, l, n, o, p, q, r, ptr; // Initialize function variables
	long NE_BR = NE_SBR + NE_FBR;
    int m, flag = 0, errchk;
    //int *jflag = alloc_int (NJ*3);
	int *jflag = alloc_int (NJ*4);
    if (jflag == NULL) {
        return 1;
    }
    long *jinc = alloc_long (NE_TR*2+NE_FR*2+NE_SH*3+NE_SBR*8+NE_FBR*8);
    if (jinc == NULL) {
        if (jflag != NULL) {
            free (jflag);
            jflag = NULL;
        }
        return 1;
    }
    long *jincloc = alloc_long (NJ*3+1);
    if (jincloc == NULL) {
        if (jflag != NULL) {
            free (jflag);
            jflag = NULL;
        }
        if (jinc != NULL) {
            free (jinc);
            jinc = NULL;
        }
        return 1;
    }
    long *xadj = alloc_long (NJ+1);
    if (xadj == NULL) {
        if (jflag != NULL) {
            free (jflag);
            jflag = NULL;
        }
        if (jinc != NULL) {
            free (jinc);
            jinc = NULL;
        }
        if (jincloc != NULL) {
            free (jincloc);
            jincloc = NULL;
        }
        return 1;
    }

    /* Joint-element connectivity flags track element-types connected to each joint;
       initialized to zeros, which assumes "floating" joint */
    /* Warping restraint flags track warping restraint at a joint (array 1) and number of
       frame elements connected to a joint (array 2); initialized to zeros, which assumes
       member ends are fixed against warping */
    for (i = 0; i < NJ*3; ++i) {
        jflag[i] = *(pwrpres+i) = 0;
    }

    // Establish truss member incidences
    for (i = 0; i < NE_TR; ++i) {
        // Read in Ends 1 and 2 from input file
        fscanf(IFP[0], "%ld,%ld\n", pminc+i*2, pminc+i*2+1);
        jflag[(*(pminc+i*2)-1)*3]++; // Truss element is connected to Joint j
        jflag[(*(pminc+i*2+1)-1)*3]++; // Truss element is connected to Joint k
    }

    // Establish frame member incidences
    ptr = NE_TR * 2;
    for (i = 0; i < NE_FR; ++i) {
        // Read in Ends 1 and 2 from input file
        fscanf(IFP[0], "%ld,%ld\n", pminc+ptr+i*2, pminc+ptr+i*2+1);
        jflag[(*(pminc+ptr+i*2)-1)*3+1]++; // Frame element is connected to Joint j
        jflag[(*(pminc+ptr+i*2+1)-1)*3+1]++; // Frame element is connected to Joint k
        /* Increment warping restraint flags on joints to reflect number of frame
           members framing into joint */
        (*(pwrpres+(*(pminc+ptr+i*2)-1)*3+1))++;
        (*(pwrpres+(*(pminc+ptr+i*2+1)-1)*3+1))++;
    }

    // Establish shell member incidences
    ptr += NE_FR * 2;
    for (i = 0; i < NE_SH; ++i) {
        // Read in Vertices 1, 2, and 3 from input file
        fscanf(IFP[0], "%ld,%ld,%ld\n", pminc+ptr+i*3, pminc+ptr+i*3+1, pminc+ptr+i*3+2);
        jflag[(*(pminc+ptr+i*3)-1)*3+2]++; // Shell element is connected to Joint j
        jflag[(*(pminc+ptr+i*3+1)-1)*3+2]++; // Shell element is connected to Joint k
        jflag[(*(pminc+ptr+i*3+2)-1)*3+2]++; // Shell element is connected to Joint l
    }
	
	// Establish brick member incidences
    ptr += NE_SH * 3;
    for (i = 0; i < NE_BR; ++i) {
        // Read in nodes 1-8 from input file
        //fscanf(IFP[0], "%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n", pminc+ptr+i*3, pminc+ptr+i*3+1, pminc+ptr+i*3+2, pminc+ptr+i*3+3, pminc+ptr+i*3+4, pminc+ptr+i*3+5, pminc+ptr+i*3+6, pminc+ptr+i*3+7);
        fscanf(IFP[0], "%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n", pminc+ptr+i*8, pminc+ptr+i*8+1, pminc+ptr+i*8+2, pminc+ptr+i*8+3, pminc+ptr+i*8+4, pminc+ptr+i*8+5, pminc+ptr+i*8+6, pminc+ptr+i*8+7);
		jflag[(*(pminc+ptr+i*8)-1)*3+3]++; // Brick element is connected to Joint j
        jflag[(*(pminc+ptr+i*8+1)-1)*3+3]++; // Brick element is connected to Joint k
        jflag[(*(pminc+ptr+i*8+2)-1)*3+3]++; // Brick element is connected to Joint l
		jflag[(*(pminc+ptr+i*8+3)-1)*3+3]++; // Brick element is connected to Joint j
        jflag[(*(pminc+ptr+i*8+4)-1)*3+3]++; // Brick element is connected to Joint k
        jflag[(*(pminc+ptr+i*8+5)-1)*3+3]++; // Brick element is connected to Joint l
		jflag[(*(pminc+ptr+i*8+6)-1)*3+3]++; // Brick element is connected to Joint j
        jflag[(*(pminc+ptr+i*8+7)-1)*3+3]++; // Brick element is connected to Joint k
    }

    // Build joint incidence and adjacency "locator arrays"
    jincloc[0] = xadj[0] = 0;

	for (i = 0; i < NJ; ++i) {
        *(pjnt+i) = i;
        jincloc[i*3+1] = jincloc[i*3] + jflag[i*3];
        jincloc[i*3+2] = jincloc[i*3+1] + jflag[i*3+1];
        jincloc[i*3+3] = jincloc[i*3+2] + jflag[i*3+2];
        xadj[i+1] = xadj[i] + jflag[i*3] + jflag[i*3+1] + 2 * jflag[i*3+2];
		jflag[i*3] = jflag[i*3+1] = jflag[i*3+2] = 0;
    }

    // Establish truss joint incidences
    for (i = 0; i < NE_TR; ++i) {
        j = *(pminc+i*2) - 1;
        k = *(pminc+i*2+1) - 1;
        jinc[jincloc[j*3]+jflag[j*3]] = i + 1;
        jinc[jincloc[k*3]+jflag[k*3]] = i + 1;
        jflag[j*3]++;
        jflag[k*3]++;
    }

    // Establish frame joint incidences
    ptr = NE_TR * 2;
    for (i = 0; i < NE_FR; ++i) {
        j = *(pminc+ptr+i*2) - 1;
        k = *(pminc+ptr+i*2+1) - 1;
        jinc[jincloc[j*3+1]+jflag[j*3+1]] = i + 1;
        jinc[jincloc[k*3+1]+jflag[k*3+1]] = i + 1;
        jflag[j*3+1]++;
        jflag[k*3+1]++;
    }

    // Establish shell joint incidences
    ptr += NE_FR * 2;
    for (i = 0; i < NE_SH; ++i) {
        j = *(pminc+ptr+i*3) - 1;
        k = *(pminc+ptr+i*3+1) - 1;
        l = *(pminc+ptr+i*3+2) - 1;
        jinc[jincloc[j*3+2]+jflag[j*3+2]] = i + 1;
        jinc[jincloc[k*3+2]+jflag[k*3+2]] = i + 1;
        jinc[jincloc[l*3+2]+jflag[l*3+2]] = i + 1;
        jflag[j*3+2]++;
        jflag[k*3+2]++;
        jflag[l*3+2]++;
    }
	
	// Establish brick joint incidences
    ptr += NE_SH * 3;
	
	for (i = 0; i < NE_BR; ++i) {
        j = *(pminc+ptr+i*8+0) - 1;
        k = *(pminc+ptr+i*8+1) - 1;
        l = *(pminc+ptr+i*8+2) - 1;
		n = *(pminc+ptr+i*8+3) - 1;
        o = *(pminc+ptr+i*8+4) - 1;
		p = *(pminc+ptr+i*8+5) - 1;
        q = *(pminc+ptr+i*8+6) - 1;
		r = *(pminc+ptr+i*8+7) - 1;
        jinc[jincloc[j*3+2]+jflag[j*3+2]] = i + 1;
        jinc[jincloc[k*3+2]+jflag[k*3+2]] = i + 1;
        jinc[jincloc[l*3+2]+jflag[l*3+2]] = i + 1;
		jinc[jincloc[n*3+2]+jflag[j*3+2]] = i + 1;
        jinc[jincloc[o*3+2]+jflag[k*3+2]] = i + 1;
        jinc[jincloc[p*3+2]+jflag[l*3+2]] = i + 1;
		jinc[jincloc[q*3+2]+jflag[j*3+2]] = i + 1;
        jinc[jincloc[r*3+2]+jflag[k*3+2]] = i + 1;
        jflag[j*3+2]++;
        jflag[k*3+2]++;
        jflag[l*3+2]++;
		jflag[n*3+2]++;
        jflag[o*3+2]++;
        jflag[p*3+2]++;
		jflag[q*3+2]++;
        jflag[r*3+2]++;
    }

    // Establish joint constraints; initialized to ones, which assumes free DOF */
    for (i = 0; i < NJ*7; ++i) {
        *(pjcode+i) = -1;
    }
		
    // Potentially all joints are constrained
    for (i = 0; i < NJ*7; ++i) {
        // Read in joint number and constraint direction from input file
        fscanf(IFP[0], "%ld,%ld\n", &j, &k);
        if (j != 0) {
            switch (k) {
                case 1:
                    fscanf(IFP[0], "\n");
                    *(pjcode+(j-1)*7+(k-1)) = 0;
                    break;
                case 2:
                    fscanf(IFP[0], "\n");
                    *(pjcode+(j-1)*7+(k-1)) = 0;
                    break;
                case 3:
                    fscanf(IFP[0], "\n");
                    *(pjcode+(j-1)*7+(k-1)) = 0;
                    break;
                case 4:
                    fscanf(IFP[0], "\n");
                    *(pjcode+(j-1)*7+(k-1)) = 0;
                    break;
                case 5:
                    fscanf(IFP[0], "\n");
                    *(pjcode+(j-1)*7+(k-1)) = 0;
                    break;
                case 6:
                    fscanf(IFP[0], "\n");
                    *(pjcode+(j-1)*7+(k-1)) = 0;
                    break;
                case 7:
                    fscanf(IFP[0], ",%d\n", pwrpres+(j-1)*3);
                    *(pjcode+(j-1)*7+(k-1)) = 0;
                    break;
                default:
                    fprintf(OFP[0], "\n***ERROR*** Joint constraint input not");
                    fprintf(OFP[0], "  recognized\n");
                    return 1;
                    break;
            }
        } else {
            break;
        }
    }

    /* Determine if joint DOF(s) is free due to no joint connection with truss
       (DOFs 1 thru 3), shell (DOFs 4 thru 6), or frame (DOF 7) elements */
     
    for (i = 0; i < NJ; ++i) {
        if (jflag[i*3+1] == 0) {
            if (jflag[i*3+2] == 0) {
                if (jflag[i*3] == 0) {
                    for (j = 0; j < 7; ++j) {
                        *(pjcode+i*7+j) = 0;
                    }
                } else {
                    for (j = 3; j < 7; ++j) {
                        *(pjcode+i*7+j) = 0;
                    }
                }
            } else {
               if (ANAFLAG != 4){
                    *(pjcode+i*7+6) = 0;
               }
            }
        }
    }

    
	if (OPTFLAG == 2) {
		// Pass control to graph function
		errchk = graph (pjnt, xadj, pjcode, pwrpres, pminc, jinc, jincloc);

		// Terminate program if errors encountered
		if (errchk == 1) {
			if (jflag != NULL) {
				free (jflag);
				jflag = NULL;
			}
			if (jinc != NULL) {
				free (jinc);
				jinc = NULL;
			}
			if (jincloc != NULL) {
				free (jincloc);
				jincloc = NULL;
			}
			if (xadj != NULL) {
				free (xadj);
				xadj = NULL;
			}
			return 1;
		}
	}

    if (NE_TR > 0) {
        // Print truss member incidences
        fprintf(OFP[0], "\nTruss Member Incidences:\n\tMember\t\tEnd-1\t\tEnd-2\n");
		for (i = 0; i < NE_TR; ++i) {
			fprintf(OFP[0], "\t%ld\t\t%ld\t\t%ld\n", i + 1, *(pminc+i*2),
                    *(pminc+i*2+1));
		}
        if (OPTFLAG == 2) {
			for (i = 0; i < NE_TR; ++i) {
				fprintf(IFP[1], "%ld,%ld\n", *(pminc+i*2), *(pminc+i*2+1));
			}
        }
    }

    if (NE_FR > 0) {
        // Print frame member incidences
        fprintf(OFP[0], "\nFrame Member Incidences:\n\tMember\t\tEnd-1\t\tEnd-2\n");
        ptr = NE_TR * 2;
		for (i = 0; i < NE_FR; ++i) {
			fprintf(OFP[0], "\t%ld\t\t%ld\t\t%ld\n", i + 1, *(pminc+ptr+i*2),
				*(pminc+ptr+i*2+1));
		}
        if (OPTFLAG == 2) {
			for (i = 0; i < NE_FR; ++i) {
				fprintf(IFP[1], "%ld,%ld\n", *(pminc+ptr+i*2), *(pminc+ptr+i*2+1));
			}
        }
    }

    if (NE_SH > 0) {
        // Print shell member incidences
        fprintf(OFP[0], "\nShell Member Incidences:\n\tElement\t\tVertex-1\tVertex-2\t");
        fprintf(OFP[0], "Vertex-3\n");
        ptr = NE_TR * 2 + NE_FR * 2;
		for (i = 0; i < NE_SH; ++i) {
			fprintf(OFP[0], "\t%ld\t\t%ld\t\t%ld\t\t%ld\n", i + 1, *(pminc+ptr+i*3),
				*(pminc+ptr+i*3+1), *(pminc+ptr+i*3+2));
		}
        if (OPTFLAG == 2) {
 			for (i = 0; i < NE_SH; ++i) {
				fprintf(IFP[1], "%ld,%ld,%ld\n", *(pminc+ptr+i*3), *(pminc+ptr+i*3+1),
					*(pminc+ptr+i*3+2));
			}
        }
    }
	
    if (NE_SBR > 0) {
        // Print brick member incidences
        fprintf(OFP[0], "\nSBrick Member Incidences:\n\tElement\t\tVertex-1\tVertex-2\tVertex-3\tVertex-4\tVertex-5\tVertex-6\tVertex-7\tVertex-8\n");
        ptr = NE_TR * 2 + NE_FR * 2 + NE_SH * 3;
		for (i = 0; i < NE_SBR; ++i) {
			fprintf(OFP[0], "\t%ld\t\t%ld\t\t%ld\t\t%ld\t\t%ld\t\t%ld\t\t%ld\t\t%ld\t\t\%ld\n", i + 1, *(pminc+ptr+i*8),
					*(pminc+ptr+i*8+1), *(pminc+ptr+i*8+2), *(pminc+ptr+i*8+3), *(pminc+ptr+i*8+4), *(pminc+ptr+i*8+5), *(pminc+ptr+i*8+6), *(pminc+ptr+i*8+7));
		}
	}
    
    if (NE_FBR > 0) {
        // Print brick member incidences
        fprintf(OFP[0], "\nFBrick Member Incidences:\n\tElement\t\tVertex-1\tVertex-2\tVertex-3\tVertex-4\tVertex-5\tVertex-6\tVertex-7\tVertex-8\n");
        ptr = NE_TR * 2 + NE_FR * 2 + NE_SH * 3 + NE_SBR * 8;
		for (i = 0; i < NE_FBR; ++i) {
			fprintf(OFP[0], "\t%ld\t\t%ld\t\t%ld\t\t%ld\t\t%ld\t\t%ld\t\t%ld\t\t%ld\t\t\%ld\n", i + 1, *(pminc+ptr+i*8),
					*(pminc+ptr+i*8+1), *(pminc+ptr+i*8+2), *(pminc+ptr+i*8+3), *(pminc+ptr+i*8+4), *(pminc+ptr+i*8+5), *(pminc+ptr+i*8+6), *(pminc+ptr+i*8+7));
		}
	}

    // Print user-prescribed joint constraints
    fprintf(OFP[0], "\nJoint Constraints:\n\tJoint\t\tDirection\tWarping (optional)\n");
    // Potentially all joints are constrained
	for (i = 0; i < NJ; ++i) {
		for (j = 0; j < 6; ++j) {
			if (*(pjcode+i*7+j) == 0) {
				fprintf(OFP[0], "\t%ld\t\t%ld\n", i + 1, j + 1);
			}
		}
		if (*(pjcode+i*7+6) == 0) {
			fprintf(OFP[0], "\t%ld\t\t%ld\t\t%d\n", i + 1, j + 1, *(pwrpres+i*3));
		}
	}
	if (OPTFLAG == 2) {
		for (i = 0; i < NJ; ++i) {
			for (j = 0; j < 6; ++j) {
				if (*(pjcode+i*7+j) == 0) {
					fprintf(IFP[1], "%ld,%ld\n", i + 1, j + 1);
				}
			}
			if (*(pjcode+i*7+6) == 0) {
				fprintf(IFP[1], "%ld,%ld,%d\n", i + 1, j + 1, *(pwrpres+i*3));
			}
		}
		fprintf(IFP[1], "0,0\n");
	}

    if (jflag != NULL) {
        free (jflag);
        jflag = NULL;
    }
    if (jinc != NULL) {
        free (jinc);
        jinc = NULL;
    }
    if (jincloc != NULL) {
        free (jincloc);
        jincloc = NULL;
    }
    if (xadj != NULL) {
        free (xadj);
        xadj = NULL;
    }
    return 0;
}

void fsi (long *pmcode, long *pjcode, long *pminc, long *pelface, long *pfsiinc) {
    
	// Initialize function variables
	long i, j, k, l, m, n, NE_BR;
    
    int nnpfsif; // Num nodes per FSI face; = 3 if solids are shells; = 4 if solid are bricks
    int nnps; // Num nodes per solid
    int nfps; // Num faces per solid
    int nsolids;
    
    if (shFSI_FLAG == 1) {
        nnpfsif = 3;
        nnps = 3;
        nfps = 1;
        nsolids = NE_SH;
    }
    
    if (brFSI_FLAG == 1) {
        nnpfsif = 4;
        nnps = 8;
        nfps = 6;
        nsolids = NE_SBR;
    }

    NE_BR = NE_FBR+nsolids;
    long ptr = nsolids*nnps;
    
	int face[4]; //global joints associated with an element's particular face
	int lcjts[4], lcjtf[4]; //local fluid and solid joints for use in evaluating f-s faces
	
	// Initialize elface and fsi incidence array
	for (i = 0; i < nsolids; ++i) {
		*(pelface+i) = 0;
		
		for (j = 0; j < nfps; ++j) {
			for (k = 0; k < nnpfsif; ++k) {
				*(pfsiinc+i*nfps*nnpfsif+j*nnpfsif+k) = 0;
			}
		}
	}
    
    if (FSIINCFLAG == 0) { // read in fsi incidence array from input file
        FILE *fsifile;
        long j1, j2, j3;
        
        // Read in solid shell coordinates
        do {
            fsifile = fopen("fsiinc.txt", "r"); // Open input file for reading
        } while (fsifile == 0);
        
        // Load dented shell coordinates
        
        for (i = 0; i < nsolids; ++i) {
            fscanf(fsifile, "%ld,%ld,%ld\n",&j1,&j2,&j3);
            if (j1 != 0 || j2 != 0 || j3 != 0) { *(pelface+i) = 1;}
            *(pfsiinc+i*3+0) = j1; *(pfsiinc+i*3+1) = j2; *(pfsiinc+i*3+2) = j3;
        }
        fclose(fsifile); // Close input file
    }
    else if (FSIINCFLAG == 1) { // auto find interface
        /* *** READING IN FSIINC*** */
        /* Determine which element faces are f-s faces by finding the faces that have a
         fluid element and a solid element in common */
        if (brFSI_FLAG == 1) {
            // Loop through fluid elements
            for (i = nsolids; i < NE_BR; ++i) {
                
                l=0; // l = current solid element
                m=0;
                
                // Loop through solid elements until a matching face for each fluid element is found
                while ((m < nnpfsif) && l < nsolids) {
                    
                    m = 0; // m counts how many nodes on a potential matching face have ben found
                    
                    for (j = 0; j < 8; ++j) { // Loop through nodes of fluid element i
                        for (k = 0; k < nnps; ++k) { // Loop through nodes of solid element l
                            
                            // If the global joints of fluid and solid element are a match
                            if (*(pminc+ptr+(i-nsolids)*8+j) == *(pminc+l*nnps+k)) {
                                
                                face[m] = *(pminc+l*nnps+k); // Assign global joint to local face array
                                ++m;
                                
                                /* If a full face has been found, increase the
                                 number of fsi faces for current solid element */
                                if (m == nnpfsif) {
                                    *(pelface+l) += 1;
                                }
                            }
                        }
                    }
                    ++l; // Move on to next solid element if matching face not found
                }
                
                // If a full face has been found, assign the joints to the fsi incidence array
                if (*(pelface+l-1) != 0 && m == nnpfsif){
                    
                    for (n = 0; n < nnpfsif; ++n) {
                        
                        // Store the global joint in the fsi incidence array
                        *(pfsiinc + (l-1)*nfps*nnpfsif + (*(pelface+l-1)-1)*nnpfsif + n) = face[n];
                    }
                }
            }
        }
        
        if (shFSI_FLAG == 1) {
            // Loop through solid elements
            for (l = 0; l < nsolids; ++l) {
                
                i = nsolids; // i = current fluid element
                m = 0;
                
                // Loop through fluid elements until a matching face for each solid element is found
                while ((m < nnpfsif) && i < NE_BR) {
                    
                    m = 0; // m counts how many nodes on a potential matching face have ben found
                    
                    for (k = 0; k < nnps; ++k) { // Loop through nodes of solid element l
                        for (j = 0; j < 8; ++j) { // Loop through nodes of fluid element i
                            
                            // If the global joints of fluid and solid element are a match
                            if (*(pminc+ptr+(i-nsolids)*8+j) == *(pminc+l*nnps+k)) {
                                
                                face[m] = *(pminc+l*nnps+k); // Assign global joint to local face array
                                ++m;
                                
                                /* If a full face has been found, increase the
                                 number of fsi faces for current solid element */
                                if (m == nnpfsif) {
                                    *(pelface+l) += 1;
                                }
                            }
                        }
                    }
                    ++i; // Move on to next fluid element if matching face not found
                }
                
                // If a full face has been found, assign the joints to the fsi incidence array
                if (*(pelface+l-1) != 0 && m == nnpfsif){
                    
                    for (n = 0; n < nnpfsif; ++n) {
                        
                        // Store the global joint in the fsi incidence array
                        *(pfsiinc + l*nfps*nnpfsif + (*(pelface+l-1)-1)*nnpfsif + n) = face[n];
                    }
                }
            }
        }
    }
}

void renumfsi (long *pminc, long *pmcode, long *pjcode, long *pelface, long *pfsiinc) {
	
	// Initialize function variables
	long i, j, k, l, m, n, o, NE_BR;
	long jt, DOF;
    
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
    
	long jcode_temp[NJ][7];
    
    NE_BR = NE_FBR+nsolids;
    
	long NEQ_repl, NEQ_prev; //used for renumbering mcode
	double fsitemp[nsolids*nfps*nnpfsif];

	// Initialize fsitemp with fsi_inc values
	for (i = 0; i < nsolids; ++i) {
		for (j = 0; j < nfps; ++j) {
			for (k = 0; k < nnpfsif; ++k) {
				fsitemp[i*nnpfsif*nfps+j*nnpfsif+k] = *(pfsiinc+i*nnpfsif*nfps+j*nnpfsif+k);
			}
		}
	}
	
	for (i = 0; i < NJ; ++i) {
		for (j = 0; j < 7; ++j) {
			jcode_temp[i][j] = 0;
		}
	}
	
    for (i = 0; i < NJ; ++i) {
        fprintf(OFP[4],"%ld\t\t\t",i+1);
        for (j = 0; j < 7; ++j) {
            fprintf(OFP[4],"%ld\t",*(pjcode+i*7+j));
        }
        fprintf(OFP[4],"\n");
    }
    
    for (i = 0; i < NE_BR; ++i) {
        fprintf(OFP[6],"%ld\t\t\t",i+1);
        for (j = 0; j < 24; ++j) {
            fprintf(OFP[6],"%ld\t",*(pmcode+i*24+j));
        }
        fprintf(OFP[6],"\n");
    }
}		
						
int graph (long *pjnt, long *pxadj, long *pjcode, int *pwrpres, long *pminc, long *pjinc,
    long *pjincloc)
{
    // Initialize function variables
    long i, j, k, l, m, n, o, ptr, ptr2, max, band, prof;
    int errchk, flag = 0;
    long *adjncy = alloc_long (*(pxadj+NJ));
    if (adjncy == NULL) {
        return 1;
    }

    // Generate DOF-weights and adjacency structure for each joint
    ptr = NE_TR * 2;
    ptr2 = NE_TR * 2 + NE_FR * 2;
    band = prof = 0;
    for (i = 0; i < NJ; ++i) {
        o = 0;
        // Generate adjacency structure related to truss members
        for (j = *(pjincloc+i*3); j < *(pjincloc+i*3+1); ++j) {
            k = *(pjinc+j) - 1;
            l = *(pminc+k*2) - 1;
            m = *(pminc+k*2+1) - 1;
            if (l == i) {
                adjncy[*(pxadj+i)+o] = m;
                o++;
            } else {
                adjncy[*(pxadj+i)+o] = l;
                o++;
            }
        }
        // Generate adjacency structure related to frame members
        for (j = *(pjincloc+i*3+1); j < *(pjincloc+i*3+2); ++j) {
            k = *(pjinc+j) - 1;
            l = *(pminc+ptr+k*2) - 1;
            m = *(pminc+ptr+k*2+1) - 1;
            if (l == i) {
                adjncy[*(pxadj+i)+o] = m;
                o++;
            } else {
                adjncy[*(pxadj+i)+o] = l;
                o++;
            }
        }
        // Generate adjacency structure related to shell members
        for (j = *(pjincloc+i*3+2); j < *(pjincloc+i*3+3); ++j) {
            k = *(pjinc+j) - 1;
            l = *(pminc+ptr2+k*3) - 1;
            m = *(pminc+ptr2+k*3+1) - 1;
            n = *(pminc+ptr2+k*3+2) - 1;
            if (l == i) {
                adjncy[*(pxadj+i)+o] = m;
                adjncy[*(pxadj+i)+o+1] = n;
                o += 2;
            } else if (m == i) {
                adjncy[*(pxadj+i)+o] = l;
                adjncy[*(pxadj+i)+o+1] = n;
                o += 2;
            } else {
                adjncy[*(pxadj+i)+o] = l;
                adjncy[*(pxadj+i)+o+1] = m;
                o += 2;
            }
        }

        /* If duplicate joints exist in adjacency structure for a given joint, eliminate
           duplicates and reduce adjacency "locator array" value by one */
        for (j = 0; j < *(pxadj+i+1) - *(pxadj+i); ++j) {
			for (k = j + 1; k < *(pxadj+i+1) - *(pxadj+i); ++k) {
				while (adjncy[*(pxadj+i)+k] == adjncy[*(pxadj+i)+j] &&
					*(pxadj+i)+k != *(pxadj+i+1)) {
					for (l = k; l < *(pxadj+i+1) - *(pxadj+i) - 1; ++l) {
						adjncy[*(pxadj+i)+l] = adjncy[*(pxadj+i)+l+1];
	                }
	                for (l = i + 1; l < NJ + 1; ++l) {
                        (*(pxadj+l))--;
					}
	            }
			}
        }

        /* Determine maximum "bandwidth" associated with joint, add contribution to
           "profile," and compare to most recent maximum "bandwidth" associated with
           structure */
        max = 0;
        for (j = *(pxadj+i); j < *(pxadj+i+1); ++j) {
            if (labs(adjncy[j] - i) > max) {
                max = labs(adjncy[j] - i);
            }
        }
        prof += max;
        if (max > band) {
            band = max;
        }
    }

    // Pass control to optnum function
    errchk = optnum (pjnt, adjncy, pxadj, &band, &prof, &flag);

    // Terminate program if errors encountered
    if (errchk == 1) {
        if (adjncy != NULL) {
            free (adjncy);
            adjncy = NULL;
        }
        return 1;
    }

    if (flag == 1) {
        fprintf(OFP[0], "\n***WARNING*** Original node-numbering scheme modified\n");

        // Pass control to updatenum function
        errchk = updatenum (pjcode, pminc, pwrpres, pjnt);

        // Terminate program if errors encountered
        if (errchk == 1) {
            if (adjncy != NULL) {
                free (adjncy);
                adjncy = NULL;
            }
            return 1;
        }
    }

    if (adjncy != NULL) {
        free (adjncy);
        adjncy = NULL;
    }
    return 0;
}

int optnum (long *pjnt, long *padjncy, long *pxadj, long *pband, long *pprof,
    int *pflag)
{
    // Initialize function variables
    long ik, i, max, k, kmax, k4, jj, k5, nband, nprof, j;
    long *newjt = alloc_long (NJ);
    if (newjt == NULL) {
        return 1;
    }
    long *joint = alloc_long (NJ);
    if (joint == NULL) {
        if (newjt != NULL) {
            free (newjt);
            newjt = NULL;
        }
        return 1;
    }

	*pflag = 0; // Initialize flag for updating node numbering

    // Attempt a node-numbering scheme setting each old node as the starting node
    for (ik = 0; ik < NJ; ++ik) {
        // Initialize sub-routine variables for the ik'th starting node
        for (i = 0; i < NJ; ++i) {
            joint[i] = newjt[i] = -1;
        }
        max = 0;
        i = 0;
        newjt[0] = ik;
        joint[ik] = 0;
        k = kmax = 0;

        do {
            // Determine number of nodes connected to i'th new node
            if (newjt[i] != -1) {
				k4 = *(pxadj+newjt[i]+1) - *(pxadj+newjt[i]);

				// Consecutively number nodes connected to i'th new node
				for (jj = 0; jj < k4; ++jj) {
					k5 = *(padjncy+*(pxadj+newjt[i])+jj);
					if (joint[k5] == -1) {
						k++;
						kmax++;
						newjt[k] = k5;
						joint[k5] = k;

						/* Calculate "bandwidth" and compare to existing; terminate attempt
						   if "bandwidth" is greater than or equal to existing */
						nband = labs(i - k);
						if (nband <= *pband) {
							if (nband > max) {
								max = nband;
							}
						} else {
							k = NJ;
							break;
						}
					}
				}
            } else if (newjt[i] == -1 && k != 0) {
				k++;
				newjt[k] = newjt[k-1] - 1;
            } else {
            	k = NJ;
            	break;
            }
            i++;
        } while (k < NJ - 1);

        /* If all joints were assigned new joint numbers, compute profile of proposed
           joint-numbering scheme */
        if (k == NJ - 1) {
            nband = max;
            nprof = 0;
            for (i = 0; i < NJ; ++i) {
                // Determine maximum "bandwidth" associated with joint
                max = 0;
                for (j = *(pxadj+joint[i]); j < *(pxadj+joint[i]+1); ++j) {
                    if (labs(*(padjncy+j) - i) > max) {
                        max = labs(*(padjncy+j) - i);
                    }
                }
                nprof += max;
                if (nprof >= *pprof) {
                    break;
                }
            }

            /* If profile of proposed joint-numbering scheme is less than current
               profile, store joint-numbering scheme */
            if (i == NJ && nprof < *pprof) {
                *pflag = 1;
                *pband = nband;
                *pprof = nprof;
                for (i = 0; i < NJ; ++i) {
                	if (joint[i] != -1) {
						*(pjnt+i) = joint[i];
                	} else {
                		kmax++;
                		*(pjnt+i) = kmax;
                	}
                }
            }
        }
    }

    if (newjt != NULL) {
        free (newjt);
        newjt = NULL;
    }
    if (joint != NULL) {
        free (joint);
        joint = NULL;
    }
    return 0;
}

int updatenum (long *pjcode, long *pminc, int *pwrpres, long *pjnt)
{
    // Define function variables
    long i, j, ptr;
    long *jcode = alloc_long (NJ*7);
    if (jcode == NULL) {
        return 1;
    }
    int *wrpres = alloc_int (NJ*3);
    if (wrpres == NULL) {
        if (jcode != NULL) {
            free (jcode);
            jcode = NULL;
        }
        return 1;
    }

    // Transfer variables with new node-numbering scheme to temporary variables
    for (i = 0; i < NJ; ++i) {
        // Transfer jcode
        for (j = 0; j < 7; ++j) {
            jcode[*(pjnt+i)*7+j] = *(pjcode+i*7+j);
        }

        // Transfer wrpres
        for (j = 0; j < 3; ++j) {
            wrpres[*(pjnt+i)*3+j] = *(wrpres+i*3+j);
        }
    }

    // Update variables with new node-numbering scheme
    for (i = 0; i < NJ; ++i) {
        // Update jcode
        for (j = 0; j < 7; ++j) {
            *(pjcode+i*7+j) = jcode[i*7+j];
        }

        // Update wrpres
        for (j = 0; j < 3; ++j) {
            *(wrpres+i*3+j) = wrpres[i*3+j];
        }
    }

    // Update truss member incidences
    for (i = 0; i < NE_TR; ++i) {
        *(pminc+i*2) = *(pjnt+*(pminc+i*2)-1) + 1;
        *(pminc+i*2+1) = *(pjnt+*(pminc+i*2+1)-1) + 1;
    }

    // Update frame member incidences
    ptr = NE_TR * 2;
    for (i = 0; i < NE_FR; ++i) {
        *(pminc+ptr+i*2) = *(pjnt+*(pminc+ptr+i*2)-1) + 1;
        *(pminc+ptr+i*2+1) = *(pjnt+*(pminc+ptr+i*2+1)-1) + 1;
    }

    // Update shell member incidences
    ptr += NE_FR * 2;
    for (i = 0; i < NE_SH; ++i) {
        *(pminc+ptr+i*3) = *(pjnt+*(pminc+ptr+i*3)-1) + 1;
        *(pminc+ptr+i*3+1) = *(pjnt+*(pminc+ptr+i*3+1)-1) + 1;
        *(pminc+ptr+i*3+2) = *(pjnt+*(pminc+ptr+i*3+2)-1) + 1;
    }

    if (jcode != NULL) {
        free (jcode);
        jcode = NULL;
    }
    if (wrpres != NULL) {
        free (wrpres);
        wrpres = NULL;
    }
    return 0;
}

void codes (long *pmcode, long *pjcode, long *pminc, int *pwrpres)
{
    long i, j, k, l, m, n, o, p, q, r, ptr, ptr2; // Initialize function variables

    if (ANAFLAG != 4){
        // Generate jcode
        NEQ = 0;
        for (i = 0; i < NJ; ++i) {
            *(pwrpres+i*3+2) = *(pwrpres+i*3+1);
            for (j = 0; j < 6; ++j) {
                if (*(pjcode+i*7+j) != 0) {
                    NEQ++;
                    *(pjcode+i*7+j) = NEQ;
                }
            }
            if (*(pjcode+i*7+6) != 0) {
                NEQ++;
                *(pjcode+i*7+6) = NEQ;
            } else if (*(pjcode+i*7+6) == 0 && *(pwrpres+i*3) == 1) {
                NEQ++;
                *(pjcode+i*7+6) = NEQ;
                NEQ += *(pwrpres+i*3+1) - 1;
            }
        }
    }

    /* Generate jcode for FSI; DOF 7 is pressure, not warping
     Number DOFs s.t. displacement dofs are first, then pressure dofs*/
    if (ANAFLAG == 4){
        NEQ = 0;
        for (i = 0; i < NJ; ++i) {
            for (j = 0; j < 6; ++j) {
                if (*(pjcode+i*7+j) != 0) {
                    NEQ++;
                    *(pjcode+i*7+j) = NEQ;
                }
            }
        }
    }

    SNDOF = NEQ;

    if (ANAFLAG == 4){
        // Label pressure dofs
        for (i = 0; i < NJ; ++i) {
            if (*(pjcode+i*7+6) != 0) {
                NEQ++;
                *(pjcode+i*7+6) = NEQ;
            }
        }
	}
    
    // Fluid dofs are total - structural
    FNDOF = NEQ-SNDOF;

    // Generate mcode for truss elements from jcode using member incidence matrix
    for (i = 0; i < NE_TR; ++i) {
        j = *(pminc+i*2) - 1; // Establish End-1 joint number
        k = *(pminc+i*2+1) - 1; // Establish End-2 joint number
        // Increment of DOFs for each end
        for (m = 0; m < 3; ++m) {
            *(pmcode+i*6+m) = *(pjcode+j*7+m); // Assign mcode entry using jcode
            *(pmcode+i*6+3+m) = *(pjcode+k*7+m); // As above, for DOFs 4 thru 6
        }
    }

    // Generate mcode for frame elements from jcode using member incidence matrix
    ptr = NE_TR * 2;
    ptr2 = NE_TR * 6;
    for (i = 0; i < NE_FR; ++i) {
        j = *(pminc+ptr+i*2) - 1; // Establish End-1 joint number
        k = *(pminc+ptr+i*2+1) - 1; // Establish End-2 joint number
        // Increment of DOFs for each end
        for (m = 0; m < 6; ++m) {
            // Assign mcode entry using jcode
            *(pmcode+ptr2+i*14+m) = *(pjcode+j*7+m);
            // As above, for DOFs 8 thru 14
            *(pmcode+ptr2+i*14+7+m) = *(pjcode+k*7+m);
        }
        // Assign mcode entry using jcode
        if (*(pwrpres+j*3) == 1) {
            *(pmcode+ptr2+i*14+6) = *(pjcode+j*7+6) + *(pwrpres+j*3+1) -
                *(pwrpres+j*3+2);
            (*(pwrpres+j*3+2))--;
        } else {
            *(pmcode+ptr2+i*14+6) = *(pjcode+j*7+6);
        }
        // As above, for DOFs 8 thru 14
        if (*(pwrpres+k*3) == 1) {
            *(pmcode+ptr2+i*14+13) = *(pjcode+k*7+6) + *(pwrpres+k*3+1) -
                *(pwrpres+k*3+2);
            (*(pwrpres+k*3+2))--;
        } else {
            *(pmcode+ptr2+i*14+13) = *(pjcode+k*7+6);
        }
    }

    // Generate mcode for shell elements from jcode using member incidence matrix
    ptr = NE_TR * 2 + NE_FR * 2;
    ptr2 = NE_TR * 6 + NE_FR * 14;
    for (i = 0; i < NE_SH; ++i) {
        j = *(pminc+ptr+i*3) - 1; // Establish Vertex-1 joint number
        k = *(pminc+ptr+i*3+1) - 1; // Establish Vertex-2 joint number
        l = *(pminc+ptr+i*3+2) - 1; // Establish Vertex-3 joint number
        // Increment of DOFs for each vertex
        for (m = 0; m < 6; ++m) {
            // Assign mcode entry using jcode
            *(pmcode+ptr2+i*18+m) = *(pjcode+j*7+m);
            // As above, for DOFs 7 thru 12
            *(pmcode+ptr2+i*18+6+m) = *(pjcode+k*7+m);
            // As above, for DOFs 13 thru 18
            *(pmcode+ptr2+i*18+12+m) = *(pjcode+l*7+m);
        }
    }

    
	// Generate mcode for solid brick elements from jcode using member incidence matrix
    ptr = NE_TR * 2 + NE_FR * 2 + NE_SH * 3;
    ptr2 = NE_TR * 6 + NE_FR * 14 + NE_SH * 18;
    for (i = 0; i < NE_SBR; ++i) {
        
        j = *(pminc+ptr+i*8+0) - 1;
        k = *(pminc+ptr+i*8+1) - 1;
        l = *(pminc+ptr+i*8+2) - 1;
		n = *(pminc+ptr+i*8+3) - 1;
        o = *(pminc+ptr+i*8+4) - 1;
		p = *(pminc+ptr+i*8+5) - 1;
        q = *(pminc+ptr+i*8+6) - 1;
		r = *(pminc+ptr+i*8+7) - 1;
        // Increment of DOFs for each vertex
        for (m = 0; m < 3; ++m) {
            // Assign mcode entry using jcode
            *(pmcode+ptr2+i*24+m) = *(pjcode+j*7+m);
            // As above, for DOFs 4 thru 6
            *(pmcode+ptr2+i*24+3+m) = *(pjcode+k*7+m);
            // As above, for DOFs 7 thru 9
            *(pmcode+ptr2+i*24+6+m) = *(pjcode+l*7+m);
            // As above, for DOFs 10 thru 12
            *(pmcode+ptr2+i*24+9+m) = *(pjcode+n*7+m);
            // As above, for DOFs 13 thru 15
            *(pmcode+ptr2+i*24+12+m) = *(pjcode+o*7+m);
            // As above, for DOFs 16 thru 18
            *(pmcode+ptr2+i*24+15+m) = *(pjcode+p*7+m);
            // As above, for DOFs 19 thru 21
            *(pmcode+ptr2+i*24+18+m) = *(pjcode+q*7+m);
            // As above, for DOFs 22 thru 24
            *(pmcode+ptr2+i*24+21+m) = *(pjcode+r*7+m);
        }
    }

    ptr = NE_TR * 2 + NE_FR * 2 + NE_SH * 3 + NE_SBR * 8;
    ptr2 = NE_TR * 6 + NE_FR * 14 + NE_SH * 18 + NE_SBR * 24;
    for (i = 0; i < NE_FBR; ++i) {
        
        j = *(pminc+ptr+i*8+0) - 1;
        k = *(pminc+ptr+i*8+1) - 1;
        l = *(pminc+ptr+i*8+2) - 1;
		n = *(pminc+ptr+i*8+3) - 1;
        o = *(pminc+ptr+i*8+4) - 1;
		p = *(pminc+ptr+i*8+5) - 1;
        q = *(pminc+ptr+i*8+6) - 1;
		r = *(pminc+ptr+i*8+7) - 1;

        // Increment of DOFs for each vertex
        for (m = 0; m < 3; ++m) {
            // x and y DOFs are fixed in fluid
            if (m < 2) {
                // Assign mcode entry using jcode
                *(pmcode+ptr2+i*24+m) = 0;
                // As above, for DOFs 4 thru 6
                *(pmcode+ptr2+i*24+3+m) = 0;
                // As above, for DOFs 7 thru 9
                *(pmcode+ptr2+i*24+6+m) = 0;
                // As above, for DOFs 10 thru 12
                *(pmcode+ptr2+i*24+9+m) = 0;
                // As above, for DOFs 13 thru 15
                *(pmcode+ptr2+i*24+12+m) = 0;
                // As above, for DOFs 16 thru 18
                *(pmcode+ptr2+i*24+15+m) = 0;
                // As above, for DOFs 19 thru 21
                *(pmcode+ptr2+i*24+18+m) = 0;
                // As above, for DOFs 22 thru 24
                *(pmcode+ptr2+i*24+21+m) = 0;
            }
            
            else { // m = 3, i.e. pressure dof is dof 7 in jcode
            // Assign mcode entry using jcode
            *(pmcode+ptr2+i*24+m) = *(pjcode+j*7+6);
            // As above, for DOFs 4 thru 6
            *(pmcode+ptr2+i*24+3+m) = *(pjcode+k*7+6);
            // As above, for DOFs 7 thru 9
            *(pmcode+ptr2+i*24+6+m) = *(pjcode+l*7+6);
            // As above, for DOFs 10 thru 12
            *(pmcode+ptr2+i*24+9+m) = *(pjcode+n*7+6);
            // As above, for DOFs 13 thru 15
            *(pmcode+ptr2+i*24+12+m) = *(pjcode+o*7+6);
            // As above, for DOFs 16 thru 18
            *(pmcode+ptr2+i*24+15+m) = *(pjcode+p*7+6);
            // As above, for DOFs 19 thru 21
            *(pmcode+ptr2+i*24+18+m) = *(pjcode+q*7+6);
            // As above, for DOFs 22 thru 24
            *(pmcode+ptr2+i*24+21+m) = *(pjcode+r*7+6);

            }
        }
    }
}

int skylin (long *pmaxa, long *pmcode, long *plss)
{
	
    long i, j, k, min, ptr; // Initialize function variables
    long *kht = alloc_long (NEQ);
    if (kht == NULL) {
        return 1;
    }

    // Zero-out kht
    for (i = 0; i < NEQ; ++i) {
        kht[i] = 0;
    }

    /* Define column height array, kht.  Each address in kht corresponds to a column in
       the stiffness matrix; the value of the address defines the skyline height above
       the diagonal entry. */
    /* Iterate over truss elements to span columns in mcode (LM array using Bathe's
       notation) */
    for (i = 0; i < NE_TR; ++i) {
        min = NEQ; // Guess at a reasonably large "minimum" to get things started
        for (j = 0; j < 6; ++j) {
            // Does mcode entry correspond to a global DOF? Smaller than min?
            if ((*(pmcode+i*6+j) > 0) && (*(pmcode+i*6+j) < min)) {
                min = *(pmcode+i*6+j); // New min...
            }
        }
        for (j = 0; j < 6; ++j) {
            k = *(pmcode+i*6+j);
            // Does the mcode entry correspond to a global DOF?
            if (k != 0) {
                /* Use the mcode to discern column height kht.  The maximum difference
                   between non-zero mcode entries, (k - min), corresponding to given
                   element, defines the column height in the stiffness matrix,
                   corresponding to the degree of freedom "k". */
                if ((k - min) > kht[k-1]) {
                    kht[k-1] = k - min;
                }
            }
        }
    }
    /* Iterate over frame elements to span columns in mcode (LM array using Bathe's
       notation) */
    ptr = NE_TR * 6;
    for (i = 0; i < NE_FR; ++i) {
        min = NEQ; // Guess at a reasonably large "minimum" to get things started
        for (j = 0; j < 14; ++j) {
            // Does mcode entry correspond to a global DOF? Smaller than min?
            if ((*(pmcode+ptr+i*14+j) > 0) && (*(pmcode+ptr+i*14+j) < min)) {
                min = *(pmcode+ptr+i*14+j); // New min...
            }
        }
        for (j = 0; j < 14; ++j) {
            k = *(pmcode+ptr+i*14+j);
            // Does the mcode entry correspond to a global DOF?
            if (k != 0) {
                /* Use the mcode to discern column height kht.  The maximum difference
                   between non-zero mcode entries, (k - min), corresponding to given
                   element, defines the column height in the stiffness matrix,
                   corresponding to the degree of freedom "k". */
                if ((k - min) > kht[k-1]) {
                    kht[k-1] = k - min;
                }
            }
        }
    }
    /* Iterate over shell elements to span columns in mcode (LM array using Bathe's
       notation) */
    ptr = NE_TR * 6 + NE_FR * 14;
    for (i = 0; i < NE_SH; ++i) {
        min = NEQ; // Guess at a reasonably large "minimum" to get things started
        for (j = 0; j < 18; ++j) {
            // Does mcode entry correspond to a global DOF? Smaller than min?
            if ((*(pmcode+ptr+i*18+j) > 0) && (*(pmcode+ptr+i*18+j) < min)) {
                min = *(pmcode+ptr+i*18+j); // New min...
            }
        }
        for (j = 0; j < 18; ++j) {
            k = *(pmcode+ptr+i*18+j);
            // Does the mcode entry correspond to a global DOF?
            if (k != 0) {
                /* Use the mcode to discern column height kht.  The maximum difference
                   between non-zero mcode entries, (k - min), corresponding to given
                   element, defines the column height in the stiffness matrix,
                   corresponding to the degree of freedom "k". */
                if ((k - min) > kht[k-1]) {
                    kht[k-1] = k - min;
                }
            }
        }
    }

    /* Generate maxa which provides the addresses in the stiffness array for the diagonal
       elements of original stiffness matrix */
    *pmaxa = 1;
    for (i = 0; i < NEQ; ++i) { // Set counter to number of diagonal elements
        /* Last diagonal term + column height + one equals new diagonal address
           Last array element is used to define length of stiffness array */
        *(pmaxa+i+1) = *(pmaxa+i) + kht[i] + 1;
    }			
    /* Length of stiffness matrix.  This is true since last element of maxa is the
       address of the final skyline element, plus 1 */
	if (SLVFLAG == 0) {
		*plss = *(pmaxa+NEQ) - 1;
	}
	else {
		*plss = NEQ*NEQ;
	}
    if (kht != NULL) {
        free(kht);
        kht = NULL;
    }
    return 0;
}

int load (double *pq, double *pefFE_ref, double *px, double *pllength, double *poffset,
    int *posflag, double *pc1, double *pc2, double *pc3, long *pjnt, long *pmcode, long *pjcode, 
    long *pminc, double *ptinpt, double *ppinpt, double *pdinpt, double *pum, double *pvm, double *pam)
{
    // Initialize function variables
    long i, j, k, l, jt, fr, ptr, DOF;
    int dir, flag;
    double mag, sum;
	//long ntstpsinpt;
    // Distributed load vector in global and local coordinate systems, respectively
    double W[3], w[3];
    /* Element fixed-end force vector in global coordinate system, w.r.t. local joints i
       and j */
    double QFij[14];
    double T[14][14]; // Coordinate transformation matrix
    double T_rl[14][14]; // Rigid link transformation matrix

    // Zero-out generalized load vector
    for (i = 0; i < NEQ; ++i) {
		*(pq+i) = 0;
    }
	
	for (i = 0; i < NEQ; ++i) {
		*(pum+i) = 0;
		*(pvm+i) = 0;
		*(pam+i) = 0;
	}

	if (ALGFLAG != 4) {
		fscanf(IFP[0], "%ld,%d,%lf\n", &jt, &dir, &mag);
		if (jt != 0) { // Check for joint loading
			flag = 0;
			fprintf(OFP[0], "\nJoint Loads:\n\tGlobal Joint\tDirection\tForce\n");
			do {
				jt = *(pjnt+jt-1) + 1;
				fprintf(OFP[0], "\t%ld\t\t%d\t\t%lf\n", jt, dir, mag);
				if (OPTFLAG == 2) {
					fprintf(IFP[1], "%ld,%d,%lf\n", jt, dir, mag);
				}
				
				k = *(pjcode+(jt-1)*7+dir-1); // Scan and load jcode
				
				// Store only joint loads corresponding to active global DOFs
				switch (k) {
					case (0):
						break; // Do not store loads at supports
					default:
						*(pq+k-1) = mag;
						break;
				}
				fscanf(IFP[0], "%ld,%d,%lf\n", &jt, &dir, &mag);
			} while (jt != 0); // Check for last joint load
		} else {
			flag = 1;
		}
		if (OPTFLAG == 2) {
			fprintf(IFP[1], "0,0,0\n");
		}
		if (NE_FR != 0) {
			// Initialize reference fixed-end force vector
			for (i = 0; i < NE_FR; ++i) {
				for (j = 0; j < 14; ++j) {
					*(pefFE_ref+i*14+j) = 0;
				}
			}
			fscanf(IFP[0], "%ld,%d,%lf\n", &fr, &dir, &mag);
			if (fr != 0) { // Check for frame element distributed load
				if (dir >= 1 && dir <= 3) {
					flag = 0;
					ptr = NE_TR * 2;
					fprintf(OFP[0], "\nDistributed Loads:\n\tFrame Element\tDirection\t");
					fprintf(OFP[0], "Magnitude\n");
					do {
						fprintf(OFP[0], "\t%ld\t\t%d\t\t%lf\n", fr, dir, mag);
						if (OPTFLAG == 2) {
							fprintf(IFP[1], "%ld,%d,%le\n", fr, dir, mag);
						}
						
						/* Assign magnitude of distributed load to appropriate location
						 within global distributed load vector */
						W[0] = W[1] = W[2] = 0;
						W[dir-1] = mag;
						
						// Initialize all elements to zero
						for (i = 0; i < 14; ++i) {
							for (j = 0; j < 14; ++j) {
								T[i][j] = T_rl[i][j] = 0;
							}
						}
						
						// Assign non-zero elements of coordinate transformation matrix
						fr--; // Decrement for convenience
						T[0][0] = T[3][3] = T[7][7] = T[10][10] = *(pc1+NE_TR+fr*3);
						T[0][1] = T[3][4] = T[7][8] = T[10][11] = *(pc1+NE_TR+fr*3+1);
						T[0][2] = T[3][5] = T[7][9] = T[10][12] = *(pc1+NE_TR+fr*3+2);
						T[1][0] = T[4][3] = T[8][7] = T[11][10] = *(pc2+NE_TR+fr*3);
						T[1][1] = T[4][4] = T[8][8] = T[11][11] = *(pc2+NE_TR+fr*3+1);
						T[1][2] = T[4][5] = T[8][9] = T[11][12] = *(pc2+NE_TR+fr*3+2);
						T[2][0] = T[5][3] = T[9][7] = T[12][10] = *(pc3+NE_TR+fr*3);
						T[2][1] = T[5][4] = T[9][8] = T[12][11] = *(pc3+NE_TR+fr*3+1);
						T[2][2] = T[5][5] = T[9][9] = T[12][12] = *(pc3+NE_TR+fr*3+2);
						T[6][6] = T[13][13] = 1;
						
						/* Transform global distributed load vector into local coordinate
						 system */
						for (i = 0; i < 3; ++i) {
							w[i] = 0;
							for (j = 0; j < 3; ++j) {
								w[i] += T[i][j] * W[j];
							}
						}
						
						/* Compute fixed-end axial and shear forces and fixed-end moments,
						 w.r.t. local coordinate system */
						*(pefFE_ref+fr*14) = *(pefFE_ref+fr*14+7) =
                        -w[0] * (*(pllength+NE_TR+fr)) / 2;
						*(pefFE_ref+fr*14+1) = *(pefFE_ref+fr*14+8) =
                        -w[1] * (*(pllength+NE_TR+fr)) / 2;
						*(pefFE_ref+fr*14+2) = *(pefFE_ref+fr*14+9) =
                        -w[2] * (*(pllength+NE_TR+fr)) / 2;
						*(pefFE_ref+fr*14+4) = w[2] * pow(*(pllength+NE_TR+fr),2) / 12;
						*(pefFE_ref+fr*14+11) = -w[2] * pow(*(pllength+NE_TR+fr),2) / 12;
						*(pefFE_ref+fr*14+5) = -w[1] * pow(*(pllength+NE_TR+fr),2) / 12;
						*(pefFE_ref+fr*14+12) = w[1] * pow(*(pllength+NE_TR+fr),2) / 12;
						*(pefFE_ref+fr*14+3) = *(pefFE_ref+fr*14+10) = *(pefFE_ref+fr*14+6) =
                        *(pefFE_ref+fr*14+13) = 0;
						
						if (*(posflag+fr) == 0) {
							for (i = 0; i < 2; ++i) {
								jt = *(pminc+ptr+fr*2+i) - 1;
								for (j = 0; j < 7; ++j) {
									k = *(pjcode+jt*7+j); // Scan and load jcode
									// Store only loads corresponding to active global DOFs
									switch (k) {
										case (0):
											break; // Do not store loads at supports
										default:
											/* Transform element fixed-end force vector from
											 local into global coordinate system */
											sum = 0;
											for (l = 0; l < 14; ++l) {
												sum += T[l][i*7+j] * (*(pefFE_ref+fr*14+l));
											}
											*(pq+k-1) -= sum;
											break;
									}
								}
							}
						} else {
							/* Transform element fixed-end force vector from local into
							 global coordinate system */
							for (i = 0; i < 14; ++i) {
								sum = 0;
								for (j = 0; j < 14; ++j) {
									sum += T[j][i] * (*(pefFE_ref+fr*14+j));
								}
								QFij[i] = sum;
							}
							
							// Assign non-zero elements of rigid link transformation matrix
							for (i = 0; i < 14; ++i) {
								T_rl[i][i] = 1;
							}
							T_rl[3][1] = -(*(poffset+fr*6+2));
							T_rl[3][2] = *(poffset+fr*6+1);
							T_rl[4][0] = *(poffset+fr*6+2);
							T_rl[4][2] = -(*(poffset+fr*6));
							T_rl[5][0] = -(*(poffset+fr*6+1));
							T_rl[5][1] = *(poffset+fr*6);
							T_rl[10][8] = -(*(poffset+fr*6+5));
							T_rl[10][9] = *(poffset+fr*6+4);
							T_rl[11][7] = *(poffset+fr*6+5);
							T_rl[11][9] = -(*(poffset+fr*6+3));
							T_rl[12][7] = -(*(poffset+fr*6+4));
							T_rl[12][8] = *(poffset+fr*6+3);
							
							for (i = 0; i < 2; ++i) {
								jt = *(pminc+ptr+fr*2+i) - 1;
								for (j = 0; j < 7; ++j) {
									k = *(pjcode+jt*7+j); // Scan and load jcode
									// Store only loads corresponding to active global DOFs
									switch (k) {
										case (0):
											break; // Do not store loads at supports
										default:
											/* Transform element fixed-end force vector to
											 global joints 1 and 2 */
											sum = 0;
											for (l = 0; l < 14; ++l) {
												sum += T_rl[i*7+j][l] * QFij[l];
											}
											*(pq+k-1) -= sum;
											break;
									}
								}
							}
						}
						fscanf(IFP[0], "%ld,%d,%lf\n", &fr, &dir, &mag);
					} while (fr != 0);
				} else {
					fprintf(OFP[0], "\n***ERROR*** Distributed moments and bi-moments not");
					fprintf(OFP[0], " yet implemented\n");
					return 1;
				}
			} else if (flag == 1) {
				fprintf(OFP[0], "\n***ERROR*** No forces present in input file\n");
				return 1;
			}
		} else if (flag == 1) {
			fprintf(OFP[0], "\n***ERROR*** No forces present in input file\n");
			return 1;
		}
	}
	else {  // Dynamic analysis
		
		//ntstpsinpt = ttot/dt + 1;
		long ks, kps;
		
		// Initialize input loads to zero
		for (i = 0; i < NEQ; ++i) {
			for (j = 0; j < ntstpsinpt; ++j) {
				*(ppinpt+i*ntstpsinpt+j) = 0;
			}
		}
		
		// Initialize first element in time array to be zero
		*(ptinpt) = 0;
		
		// Assemble time array 
		for (i = 1; i < ntstpsinpt; ++i) {
			*(ptinpt+i) = *(ptinpt+i-1) + dt;
		}	
		
		/* Scan and load user input applied forces.  Store only the applied forces corresponding 
		 to active DOFs */
		kps = 0;
		i = 0;
		fscanf(IFP[0], "%ld,%d,%lf\n", &jt, &dir, &mag);
		if (jt != 0) { // Check for joint loading
			
			do {
				ks = *(pjcode+(jt-1)*7+dir-1);
				
				if (ks != kps) {
					i = 0;
				} 
				
				// Store only joint loads corresponding to active solid DOFs
				switch (ks) {
					case (0):
						break; // Do not store loads at supports
					default:
						*(ppinpt+(ks-1)*ntstpsinpt+i) = mag;
						break;
				}
				kps = ks; 
				
				++i;
				
				fscanf(IFP[0], "%ld,%d,%lf\n", &jt, &dir, &mag);
			} while (jt != 0); // Check for last joint load
		}
		
		double d,v,a;
		/* Scan and load initial displ, vel and acc.  If dir = 4 u, v, and a refer to pressure dof
		 the its derivatives */
		fscanf(IFP[0], "%ld,%d,%lf,%lf,%lf\n", &jt, &dir, &d, &v, &a);
		if (jt != 0) { // Check for joint loading
			do {
                
				k = *(pjcode+(jt-1)*7+dir-1); // Scan and load jcode
				
				// Store only joint loads corresponding to active global DOFs
				switch (k) {
					case (0):
						break; // Do not store loads at supports
					default:
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
	}
	return 0;
}
