#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#define OUTPUT "master.pvd" // Map of path to input file

//********************************************************************************//
//**																			**//
//**  BenPost Version 1.2                                                       **//
//**																			**//
//**  Copyright (c) 2016 C. J. Earls                                            **//
//**  Developed by C. J. Earls, Cornell University                              **//
//**  All rights reserved.														**//
//**                                                                            **//
//**  Contributors:                                                             **//
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

// Define file pointers.
FILE *ifp[5]; // Input file
FILE *ofp; // Output file
FILE *ofp2; // Output file

// Miscellaneous functions:
int closeio (int flag);

// Memory management functions:
long * alloc_long (long arraylen);
double * alloc_dbl (long arraylen);
int free_all (long **pp2p2l, int nl, double **pp2p2d, int nd);

int main (void)
{
    // Initialize function variables
    long nj, ne_tr, ne_fr, ne_sh, neq;	
	long ne_sbr = 0;
	long ne_fbr = 0;
    long i, j, k, ptr;
    double lpf;
    char junk_char[20], file[25];
    int fsi_flag;
    
    for (i = 0; i < 5; ++i) {
        sprintf(file, "results%ld.txt", i + 1);
        // Force input file to open
        do {
            ifp[i] = fopen(file, "r");
        } while (ifp[i] == 0);
    }
    // Force output file to open
    do {
        ofp = fopen(OUTPUT, "w");
    } while (ofp == 0);
    
	
	// Begin skipping through data in input file
    for (i = 0; i < 2; ++i) {
        fscanf(ifp[0], "%s", junk_char);
    }
    fscanf(ifp[0], "%s", junk_char);
    
    // Determine whether FSI analysis
    if (strcmp(junk_char, "Fluid") == 0) {fsi_flag = 1;}
    else {fsi_flag = 0;}
    printf("fsi_flag = %d\n",fsi_flag);
    fsi_flag = 1;

    // Resume skipping through data in input file
    for (i = 0; i < 4; ++i) {
        fscanf(ifp[0], "%s", junk_char);
    }
    fscanf(ifp[0], "%s", junk_char);
    
    if (strcmp(junk_char, "Direct") == 0 || strcmp(junk_char, "Newton") == 0 || strcmp(junk_char, "Dynamic") == 0) {
        for (i = 0; i < 6; ++i) {
            fscanf(ifp[0], "%s", junk_char);
        }
    } else if (strcmp(junk_char, "Modified") == 0) {
        fscanf(ifp[0], "%s", junk_char);
        if (strcmp(junk_char, "Newton") == 0) {
            for (i = 0; i < 6; ++i) {
                fscanf(ifp[0], "%s", junk_char);
            }
        } else if (strcmp(junk_char, "Spherical") == 0) {
            for (i = 0; i < 7; ++i) {
                fscanf(ifp[0], "%s", junk_char);
            }
        }
    }
    
    // Read number of joints and elements from input file
    fscanf(ifp[0], "%ld", &nj);
    for (i = 0; i < 4; ++i) {
        fscanf(ifp[0], "%s", junk_char);
    }
    fscanf(ifp[0], "%ld", &ne_tr);
    for (i = 0; i < 4; ++i) {
        fscanf(ifp[0], "%s", junk_char);
    }
    fscanf(ifp[0], "%ld", &ne_fr);
    for (i = 0; i < 4; ++i) {
        fscanf(ifp[0], "%s", junk_char);
    }
    fscanf(ifp[0], "%ld", &ne_sh);
	for (i = 0; i < 5; ++i) {
		fscanf(ifp[0], "%s", junk_char);
	}
	fscanf(ifp[0], "%ld", &ne_sbr);
	for (i = 0; i < 5; ++i) {
		fscanf(ifp[0], "%s", junk_char);
	}
	fscanf(ifp[0], "%ld", &ne_fbr);

    // Memory management variables
    /* Pointer-to-pointer-to-double array (2 arrays of type double are defined during
     program execution) */
    double *p2p2d[5];
    // Counter to track number of arrays of type double for which memory is allocated
    int nd = 0;
    /* Pointer-to-pointer-to-long array (2 arrays of type long are defined during program
     execution) */
    long *p2p2l[2];
    // Counter to track number of arrays of type long for which memory is allocated
    int nl = 0;
    
    // Joint coordinates
    double *x = alloc_dbl (nj*3);
    if (x == NULL) {
        // Pass control to closeio function
        getchar();
        return closeio(0);
    }
    p2p2d[nd] = x;
    nd++;
    // Joint constraint code
    long *jcode = alloc_long (nj*7);
    if (jcode == NULL) {
        // Pass control to free_all function
        return free_all (p2p2l, nl, p2p2d, nd);
    }
    p2p2l[nl] = jcode;
    nl++;
    // Member incidence
    long *minc = alloc_long (ne_tr*2+ne_fr*2+ne_sh*3+ne_sbr*8+ne_fbr*8);
    if (minc == NULL) {
        // Pass control to free_all function
        return free_all (p2p2l, nl, p2p2d, nd);
    }
    p2p2l[nl] = minc;
    nl++;
    // Truss element forces
    double *ef_tr = alloc_dbl (ne_tr);
    if (ef_tr == NULL) {
        // Pass control to free_all function
        return free_all (p2p2l, nl, p2p2d, nd);
    }
    p2p2d[nd] = ef_tr;
    nd++;
    // Frame element forces
    double *ef_fr = alloc_dbl (ne_fr * 7);
    if (ef_fr == NULL) {
        // Pass control to free_all function
        return free_all (p2p2l, nl, p2p2d, nd);
    }
    p2p2d[nd] = ef_fr;
    nd++;
    // Shell element forces
    double *ef_sh = alloc_dbl (ne_sh * 6);
    if (ef_sh == NULL) {
        // Pass control to free_all function
        return free_all (p2p2l, nl, p2p2d, nd);
    }
    p2p2d[nd] = ef_sh;
    nd++;
    
    // Resume skipping through data in input file
    fscanf(ifp[0], "%s", junk_char);
    if (strcmp(junk_char, "***WARNING***") == 0) {
        for (i = 0; i < 5; ++i) {
            fscanf(ifp[0], "%s", junk_char);
        }
	}
    
    // Read truss member incidence data from input file
    if (ne_tr > 0) {
        for (i = 0; i < 5; ++i) {
            fscanf(ifp[0], "%s", junk_char);
        }
        for (i = 0; i < ne_tr; ++i) {
            fscanf(ifp[0], "%s", junk_char);
            fscanf(ifp[0], "%ld\t%ld", &minc[i*2], &minc[i*2+1]);
        }
    }
    
    // Read frame member incidence data from input file
    if (ne_fr > 0 && ne_tr > 0) {
        ptr = ne_tr * 2;
        for (i = 0; i < 6; ++i) {
            fscanf(ifp[0], "%s", junk_char);
        }
        for (i = 0; i < ne_fr; ++i) {
            fscanf(ifp[0], "%s", junk_char);
            fscanf(ifp[0], "%ld\t%ld", &minc[ptr+i*2], &minc[ptr+i*2+1]);
        }
    } else if (ne_fr > 0) {
        ptr = ne_tr * 2;
        for (i = 0; i < 5; ++i) {
            fscanf(ifp[0], "%s", junk_char);
        }
        for (i = 0; i < ne_fr; ++i) {
            fscanf(ifp[0], "%s", junk_char);
            fscanf(ifp[0], "%ld\t%ld", &minc[ptr+i*2], &minc[ptr+i*2+1]);
        }
	}
    
    // Read shell member incidence data from input file
    if (ne_sh > 0 && (ne_tr > 0 || ne_fr > 0)) {
        ptr = ne_tr * 2 + ne_fr * 2;
        for (i = 0; i < 7; ++i) {
            fscanf(ifp[0], "%s", junk_char);
        }
        for (i = 0; i < ne_sh; ++i) {
            fscanf(ifp[0], "%s", junk_char);
            fscanf(ifp[0], "%ld\t%ld\t%ld", &minc[ptr+i*3], &minc[ptr+i*3+1],
				   &minc[ptr+i*3+2]);
        }
    } else if (ne_sh > 0) {
        ptr = ne_tr * 2 + ne_fr * 2;
        for (i = 0; i < 6; ++i) {
            fscanf(ifp[0], "%s", junk_char);
        }
        for (i = 0; i < ne_sh; ++i) {
            fscanf(ifp[0], "%s", junk_char);
            fscanf(ifp[0], "%ld\t%ld\t%ld", &minc[ptr+i*3], &minc[ptr+i*3+1],
				   &minc[ptr+i*3+2]);
        }
    }
	
	// Read solid brick member incidence data from input file
    if (ne_sbr > 0 && (ne_tr > 0 || ne_fr > 0 || ne_sh > 0)) {
        ptr = ne_tr * 2 + ne_fr * 2 + ne_sh * 3;
        for (i = 0; i < 12; ++i) {
            fscanf(ifp[0], "%s", junk_char);
        }
        for (i = 0; i < ne_sbr; ++i) {
            fscanf(ifp[0], "%s", junk_char);
            fscanf(ifp[0], "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld", &minc[ptr+i*8], &minc[ptr+i*8+1],
				   &minc[ptr+i*8+2], &minc[ptr+i*8+3], &minc[ptr+i*8+4], 
				   &minc[ptr+i*8+5], &minc[ptr+i*8+6], &minc[ptr+i*8+7]);			
        }
    } else if (ne_sbr > 0) {
        ptr = ne_tr * 2 + ne_fr * 2 + ne_sh * 3;
        for (i = 0; i < 11; ++i) {
            fscanf(ifp[0], "%s", junk_char);
			printf("%s\n", junk_char);
        }
        for (i = 0; i < ne_sbr; ++i) {
            fscanf(ifp[0], "%s", junk_char);
            fscanf(ifp[0], "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld", &minc[ptr+i*8], &minc[ptr+i*8+1],
				   &minc[ptr+i*8+2], &minc[ptr+i*8+3], &minc[ptr+i*8+4], 
				   &minc[ptr+i*8+5], &minc[ptr+i*8+6], &minc[ptr+i*8+7]);
        }
    }
    
	// Read fluid brick member incidence data from input file
    if (ne_fbr > 0 && (ne_tr > 0 || ne_fr > 0 || ne_sh > 0 || ne_sbr > 0)) {
        ptr = ne_tr * 2 + ne_fr * 2 + ne_sh * 3 + ne_sbr * 8;
        for (i = 0; i < 12; ++i) {
            fscanf(ifp[0], "%s", junk_char);
        }
        for (i = 0; i < ne_fbr; ++i) {
            fscanf(ifp[0], "%s", junk_char);
            fscanf(ifp[0], "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld", &minc[ptr+i*8], &minc[ptr+i*8+1],
				   &minc[ptr+i*8+2], &minc[ptr+i*8+3], &minc[ptr+i*8+4], 
				   &minc[ptr+i*8+5], &minc[ptr+i*8+6], &minc[ptr+i*8+7]);
            
        }
    } else if (ne_fbr > 0) {
        ptr = ne_tr * 2 + ne_fr * 2 + ne_sh * 3 + ne_sbr * 8;
        for (i = 0; i < 11; ++i) {
            fscanf(ifp[0], "%s", junk_char);
        }
        for (i = 0; i < ne_fbr; ++i) {
            fscanf(ifp[0], "%s", junk_char);
            fscanf(ifp[0], "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld", &minc[ptr+i*8], &minc[ptr+i*8+1],
				   &minc[ptr+i*8+2], &minc[ptr+i*8+3], &minc[ptr+i*8+4], 
				   &minc[ptr+i*8+5], &minc[ptr+i*8+6], &minc[ptr+i*8+7]);
        }
    }
    
    // Resume skipping through data in input file
    do {
        fscanf(ifp[0], "%s", junk_char);
    } while (strcmp(junk_char, "Number") != 0);
    
    for (i = 0; i < 4; ++i) {
        fscanf(ifp[0], "%s", junk_char);
    }
    fscanf(ifp[0], "%ld", &neq);
    for (i = 0; i < 11; ++i) {
        fscanf(ifp[0], "%s", junk_char);
    }
    
    // Total generalized nodal displacement
    double *d = alloc_dbl (neq);
    if (d == NULL) {
        // Pass control to free_all function
        return free_all (p2p2l, nl, p2p2d, nd);
    }
    p2p2d[nd] = d;
    nd++;
    
    // Read initial joint coordinates from input file
    for (i = 0; i < nj; ++i) {
        fscanf(ifp[0], "%s", junk_char);
        fscanf(ifp[0], "%lf\t%lf\t%lf", &x[i*3], &x[i*3+1], &x[i*3+2]);
    }
    
    // Resume skipping through data in input file
    for (i = 0; i < 12; ++i) {
        fscanf(ifp[0], "%s", junk_char);
    }
    
    // Read jcode from input file
    for (i = 0; i < nj; ++i) {
        fscanf(ifp[0], "%s", junk_char);
        fscanf(ifp[0], "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld", &jcode[i*7], &jcode[i*7+1],
               &jcode[i*7+2], &jcode[i*7+3], &jcode[i*7+4], &jcode[i*7+5], &jcode[i*7+6]);
    }
    
    // Write front matter for *.pvd file
	printf("**Begin writing to master.pvd**\n");
    fprintf(ofp, "<VTKFile type=\"Collection\" version=\"0.1\"");
    fprintf(ofp, " byte_order=\"LittleEndian\">\n");
	printf("**Finish writing to master.pvd**\n");    fprintf(ofp, " <Collection>\n");
    
    //  Create directory for storage of vtu files
    if (mkdir("./vtu_files", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
        rmdir("./vtu_files", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir("./vtu_files", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
  
	printf("**Begin writing to slave_0.vtu**\n");
    sprintf(file, "vtu_files//slave_0.vtu");
    printf("hey\n");
    // Force output file to open
    do {
        ofp2 = fopen(file, "w");
    } while (ofp2 == 0);
    
    // Write front matter for *.vtu file
    fprintf(ofp2, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
    fprintf(ofp2, " byte_order=\"LittleEndian\">\n");
    fprintf(ofp2, " <UnstructuredGrid>\n");
    fprintf(ofp2, "  <Piece NumberOfPoints=\"%ld\"", nj);
    fprintf(ofp2, " NumberOfCells=\"%ld\">\n", ne_tr + ne_fr + ne_sh + ne_sbr + ne_fbr);
    fprintf(ofp2, "   <PointData Scalars=\"scalars\">\n");
    
    // Write initial load proportioniality factor data to output file
    fprintf(ofp2, "    <DataArray type=\"Float32\"");
    if (fsi_flag == 0) {fprintf(ofp2, " Name=\"Load Proportionality Factor\">\n");}
    else {fprintf(ofp2, " Name=\"Time\">\n");}
    for (i = 0; i < nj; ++i) {
        fprintf(ofp2, "     %0.15lf\n", 0.0);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
	printf("**Write initial x-translation data to output file**\n");
    // Write initial x-translation data to output file
    fprintf(ofp2, "    <DataArray type=\"Float32\"");
    fprintf(ofp2, " Name=\"X-Translation\">\n");
    for (i = 0; i < nj; ++i) {
        fprintf(ofp2, "     %0.15lf\n", 0.0);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
    printf("**Write initial y-translation data to output file**\n");
	// Write initial y-translation data to output file
    fprintf(ofp2, "    <DataArray type=\"Float32\"");
    fprintf(ofp2, " Name=\"Y-Translation\">\n");
    for (i = 0; i < nj; ++i) {
        fprintf(ofp2, "     %0.15lf\n", 0.0);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
    printf("**Write initial z-translation data to output file**\n");
	// Write initial z-translation data to output file
    fprintf(ofp2, "    <DataArray type=\"Float32\"");
    fprintf(ofp2, " Name=\"Z-Translation\">\n");
    for (i = 0; i < nj; ++i) {
        fprintf(ofp2, "     %0.15lf\n", 0.0);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
    printf("**Write initial x-rotation data to output file**\n");
	// Write initial x-rotation data to output file
    fprintf(ofp2, "    <DataArray type=\"Float32\"");
    fprintf(ofp2, " Name=\"X-Rotation\">\n");
    for (i = 0; i < nj; ++i) {
        fprintf(ofp2, "     %0.15lf\n", 0.0);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
	printf("**Write initial y-rotation data to output file**\n");
	// Write initial y-rotation data to output file
    fprintf(ofp2, "    <DataArray type=\"Float32\"");
    fprintf(ofp2, " Name=\"Y-Rotation\">\n");
    for (i = 0; i < nj; ++i) {
        fprintf(ofp2, "     %0.15lf\n", 0.0);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
	printf("**Write initial z-rotation data to output file**\n");
	// Write initial z-rotation data to output file
    fprintf(ofp2, "    <DataArray type=\"Float32\"");
    fprintf(ofp2, " Name=\"Z-Rotation\">\n");
    for (i = 0; i < nj; ++i) {
        fprintf(ofp2, "     %0.15lf\n", 0.0);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
    // Write initial warping data to output file
    fprintf(ofp2, "    <DataArray type=\"Float32\"");
    if (fsi_flag == 0) {fprintf(ofp2, " Name=\"Warping\">\n");}
    else {fprintf(ofp2, " Name=\"Pressure\">\n");}
    for (i = 0; i < nj; ++i) {
        fprintf(ofp2, "     %0.15lf\n", 0.0);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
    fprintf(ofp2, "   </PointData>\n");
    
    // Write element force data if non fsi analysis
    if (fsi_flag == 0) {
        
        fprintf(ofp2, "   <CellData>\n");
        
        // Write initial x-force data to output file
        fprintf(ofp2, "    <DataArray type=\"Float32\"");
        fprintf(ofp2, " Name=\"X-Force\">\n");
        
        for (i = 0; i < ne_tr + ne_fr + ne_sh + ne_sbr + ne_fbr; ++i) {
            for (i = 0; i < ne_tr + ne_fr + ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write initial y-force data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"Y-Force\">\n");
            for (i = 0; i < ne_tr + ne_fr + ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write initial z-force data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"Z-Force\">\n");
            for (i = 0; i < ne_tr + ne_fr + ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write initial x-moment data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"X-Moment\">\n");
            for (i = 0; i < ne_tr + ne_fr + ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write initial y-moment data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"Y-Moment\">\n");
            for (i = 0; i < ne_tr + ne_fr + ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write initial z-moment data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"Z-Moment\">\n");
            for (i = 0; i < ne_tr + ne_fr + ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write initial bi-moment data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"Bi-Moment\">\n");
            for (i = 0; i < ne_tr + ne_fr + ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            fprintf(ofp2, "   </CellData>\n");
            
        }
    }
    
    fprintf(ofp2, "   <Points>\n");
    
    // Write initial joint coordinates to output files
    fprintf(ofp2, "    <DataArray type=\"Float32\" Name=\"\"");
    fprintf(ofp2, " NumberOfComponents=\"3\">\n");
    for (i = 0; i < nj; ++i) {
        fprintf(ofp2, "     %0.15lf %0.15lf %0.15lf\n", x[i*3], x[i*3+1], x[i*3+2]);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
    fprintf(ofp2, "   </Points>\n");
    fprintf(ofp2, "   <Cells>\n");
    
    // Write truss member incidence data to output file
    fprintf(ofp2, "    <DataArray type=\"Int32\" Name=\"connectivity\">\n");
    for (i = 0; i < ne_tr; ++i) {
        fprintf(ofp2, "     %ld %ld\n", minc[i*2] - 1, minc[i*2+1] - 1);
    }
    
    // Write frame member incidence data to output file
    ptr = ne_tr * 2;
    for (i = 0; i < ne_fr; ++i) {
        fprintf(ofp2, "     %ld %ld\n", minc[ptr+i*2] - 1, minc[ptr+i*2+1] - 1);
    }
    
    // Write shell member incidence data to output file
    ptr += ne_fr * 2;
    for (i = 0; i < ne_sh; ++i) {
        fprintf(ofp2, "     %ld %ld %ld\n", minc[ptr+i*3] - 1, minc[ptr+i*3+1] - 1,
                minc[ptr+i*3+2] - 1);
    }
    
    // Write solid brick element incidence data to output file
    ptr += ne_sh * 3;
    for (i = 0; i < ne_sbr; ++i) {
        fprintf(ofp2, "     %ld %ld %ld %ld %ld %ld %ld %ld\n", minc[ptr+i*8] - 1, minc[ptr+i*8+1] - 1,
                minc[ptr+i*8+2] - 1, minc[ptr+i*8+3] - 1, minc[ptr+i*8+4] - 1, minc[ptr+i*8+5] - 1,
                minc[ptr+i*8+6] - 1, minc[ptr+i*8+7] - 1);
    }
    
    // Write fluid brick element incidence data to output file
    ptr += ne_sbr * 8;
    for (i = 0; i < ne_fbr; ++i) {
        fprintf(ofp2, "     %ld %ld %ld %ld %ld %ld %ld %ld\n", minc[ptr+i*8] - 1, minc[ptr+i*8+1] - 1,
                minc[ptr+i*8+2] - 1, minc[ptr+i*8+3] - 1, minc[ptr+i*8+4] - 1, minc[ptr+i*8+5] - 1,
                minc[ptr+i*8+6] - 1, minc[ptr+i*8+7] - 1);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
    // Write member incidence offsets to output file
    fprintf(ofp2, "    <DataArray type=\"Int32\" Name=\"offsets\">\n");
    for (i = 0; i < ne_tr; ++i) {
        fprintf(ofp2, "     %ld\n", i * 2 + 2);
    }
    ptr = ne_tr * 2;
    for (i = 0; i < ne_fr; ++i) {
        fprintf(ofp2, "     %ld\n", ptr + i * 2 + 2);
    }
    ptr += ne_fr * 2;
    for (i = 0; i < ne_sh; ++i) {
        fprintf(ofp2, "     %ld\n", ptr + i * 3 + 3);
    }
    ptr += ne_sh * 3;
    for (i = 0; i < ne_sbr; ++i) {
        fprintf(ofp2, "     %ld\n", ptr + i * 8 + 8);
    }
    ptr += ne_sbr * 8;
    for (i = 0; i < ne_fbr; ++i) {
        fprintf(ofp2, "     %ld\n", ptr + i * 8 + 8);
    }
    fprintf(ofp2, "    </DataArray>\n");
    
    // Write cell type data to output file
    fprintf(ofp2, "    <DataArray type=\"UInt8\" Name=\"types\">\n");
    for (i = 0; i < ne_tr; ++i) {
        fprintf(ofp2, "     3\n");
    }
    for (i = 0; i < ne_fr; ++i) {
        fprintf(ofp2, "     3\n");
    }
    for (i = 0; i < ne_sh; ++i) {
        fprintf(ofp2, "     5\n");
    }
    for (i = 0; i < ne_sbr; ++i) {
        fprintf(ofp2, "     12\n");
    }
    for (i = 0; i < ne_fbr; ++i) {
        fprintf(ofp2, "     12\n");
    }
    fprintf(ofp2, "    </DataArray>\n");
    
    // Write back matter for *.vtu file
    fprintf(ofp2, "   </Cells>\n");
    fprintf(ofp2, "  </Piece>\n");
    fprintf(ofp2, " </UnstructuredGrid>\n");
    fprintf(ofp2, "</VTKFile>\n");
    
    closeio(1);
    
    // Write *vtu file data to *.pvd file
    fprintf(ofp, "  <DataSet timestep=\"0\" part=\"0\"");
    fprintf(ofp, " file=\"vtu_files/slave_0.vtu\"/>\n");
    
    // Begin skipping through data in model displacements input file
    for (i = 0; i < neq * 2 + 4; ++i) {
        fscanf(ifp[1], "%s", junk_char);
    }
    
    // Begin skipping through data in truss element forces input file
    for (i = 0; i < ne_tr * 2 + 6; ++i) {
        fscanf(ifp[2], "%s", junk_char);
    }
    
    // Begin skipping through data in frame element forces input file
    for (i = 0; i < ne_fr * 9 + 6; ++i) {
        fscanf(ifp[3], "%s", junk_char);
    }
    
    // Begin skipping through data in shell element forces input file
    for (i = 0; i < ne_sh * 8 + 6; ++i) {
        fscanf(ifp[4], "%s", junk_char);
    }
    
    // Step through load increments until solution was terminated
    k = 1;
    do {
        sprintf(file, "vtu_files//slave_%ld.vtu", k);
        // Force output file to open
        do {
            ofp2 = fopen(file, "w");
        } while (ofp2 == NULL);
        
        // Write front matter for *.vtu file
        fprintf(ofp2, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
        fprintf(ofp2, " byte_order=\"LittleEndian\">\n");
        fprintf(ofp2, " <UnstructuredGrid>\n");
        fprintf(ofp2, "  <Piece NumberOfPoints=\"%ld\"", nj);
        fprintf(ofp2, " NumberOfCells=\"%ld\">\n", ne_tr + ne_fr + ne_sh + ne_sbr + ne_fbr);
        fprintf(ofp2, "   <PointData Scalars=\"scalars\">\n");
        
        // Read load proportionality factor from model displacements input file
        fscanf(ifp[1], "%lf", &lpf);
        
        // Write load proportioniality factor data to output file
        fprintf(ofp2, "    <DataArray type=\"Float32\"");
        if (fsi_flag == 0) {fprintf(ofp2, " Name=\"Load Proportionality Factor\">\n");}
        else {fprintf(ofp2, " Name=\"Time\">\n");}
        for (i = 0; i < nj; ++i) {
            fprintf(ofp2, "     %0.15lf\n", lpf);
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        // Read total generalized nodal displacements from model displacements input file
        fscanf(ifp[1], "%s", junk_char);
        for (i = 0; i < neq; ++i) {
            fscanf(ifp[1], "%lf", &d[i]);
        }
        
        // Write x-translation data to output file
        fprintf(ofp2, "    <DataArray type=\"Float32\"");
        fprintf(ofp2, " Name=\"X-Translation\">\n");
        for (i = 0; i < nj; ++i) {
            if (jcode[i*7] != 0) {
                fprintf(ofp2, "     %0.15lf\n", d[jcode[i*7]-1]);
            } else {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        // Write y-translation data to output file
        fprintf(ofp2, "    <DataArray type=\"Float32\"");
        fprintf(ofp2, " Name=\"Y-Translation\">\n");
        for (i = 0; i < nj; ++i) {
            if (jcode[i*7+1] != 0) {
                fprintf(ofp2, "     %0.15lf\n", d[jcode[i*7+1]-1]);
            } else {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        // Write z-translation data to output file
        fprintf(ofp2, "    <DataArray type=\"Float32\"");
        fprintf(ofp2, " Name=\"Z-Translation\">\n");
        for (i = 0; i < nj; ++i) {
            if (jcode[i*7+2] != 0) {
                fprintf(ofp2, "     %0.15lf\n", d[jcode[i*7+2]-1]);
            } else {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        // Write x-rotation data to output file
        fprintf(ofp2, "    <DataArray type=\"Float32\"");
        fprintf(ofp2, " Name=\"X-Rotation\">\n");
        for (i = 0; i < nj; ++i) {
            if (jcode[i*7+3] != 0) {
                fprintf(ofp2, "     %0.15lf\n", d[jcode[i*7+3]-1]);
            } else {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        // Write y-rotation data to output file
        fprintf(ofp2, "    <DataArray type=\"Float32\"");
        fprintf(ofp2, " Name=\"Y-Rotation\">\n");
        for (i = 0; i < nj; ++i) {
            if (jcode[i*7+4] != 0) {
                fprintf(ofp2, "     %0.15lf\n", d[jcode[i*7+4]-1]);
            } else {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        // Write z-rotation data to output file
        fprintf(ofp2, "    <DataArray type=\"Float32\"");
        fprintf(ofp2, " Name=\"Z-Rotation\">\n");
        for (i = 0; i < nj; ++i) {
            if (jcode[i*7+5] != 0) {
                fprintf(ofp2, "     %0.15lf\n", d[jcode[i*7+5]-1]);
            } else {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        // Write warping data to output file
        fprintf(ofp2, "    <DataArray type=\"Float32\"");
        if (fsi_flag == 0) {fprintf(ofp2, " Name=\"Warping\">\n");}
        else {fprintf(ofp2, " Name=\"Pressure\">\n");}
        for (i = 0; i < nj; ++i) {
            if (jcode[i*7+6] != 0) {
                fprintf(ofp2, "     %0.15lf\n", d[jcode[i*7+6]-1]);
            } else {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        fprintf(ofp2, "   </PointData>\n");
        
        if (fsi_flag == 0) {
            
            fprintf(ofp2, "   <CellData>\n");
            
            // Read truss element forces from truss element forces input file
            fscanf(ifp[2], "%s", junk_char);
            fscanf(ifp[2], "%s", junk_char);
            for (i = 0; i < ne_tr; ++i) {
                fscanf(ifp[2], "%lf", &ef_tr[i]);
            }
            
            // Read frame element forces from frame element forces input file
            fscanf(ifp[3], "%s", junk_char);
            fscanf(ifp[3], "%s", junk_char);
            for (i = 0; i < ne_fr; ++i) {
                for (j = 0; j < 7; ++j) {
                    fscanf(ifp[3], "%lf", &ef_fr[i*7+j]);
                }
            }
            
            // Read shell element forces from shell element forces input file
            fscanf(ifp[4], "%s", junk_char);
            fscanf(ifp[4], "%s", junk_char);
            for (i = 0; i < ne_sh; ++i) {
                for (j = 0; j < 6; ++j) {
                    fscanf(ifp[4], "%lf", &ef_sh[i*6+j]);
                }
            }
            
            // Write x-force data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"X-Force\">\n");
            for (i = 0; i < ne_tr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_tr[i]);
            }
            for (i = 0; i < ne_fr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_fr[i*7]);
            }
            for (i = 0; i < ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_sh[i*6]);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write y-force data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"Y-Force\">\n");
            for (i = 0; i < ne_tr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            for (i = 0; i < ne_fr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_fr[i*7+1]);
            }
            for (i = 0; i < ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_sh[i*6+1]);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write z-force data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"Z-Force\">\n");
            for (i = 0; i < ne_tr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            for (i = 0; i < ne_fr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_fr[i*7+2]);
            }
            for (i = 0; i < ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_sh[i*6+2]);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write x-moment data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"X-Moment\">\n");
            for (i = 0; i < ne_tr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            for (i = 0; i < ne_fr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_fr[i*7+3]);
            }
            for (i = 0; i < ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_sh[i*6+3]);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write y-moment data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"Y-Moment\">\n");
            for (i = 0; i < ne_tr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            for (i = 0; i < ne_fr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_fr[i*7+4]);
            }
            for (i = 0; i < ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_sh[i*6+4]);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write z-moment data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"Z-Moment\">\n");
            for (i = 0; i < ne_tr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            for (i = 0; i < ne_fr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_fr[i*7+5]);
            }
            for (i = 0; i < ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_sh[i*6+5]);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            // Write bi-moment data to output file
            fprintf(ofp2, "    <DataArray type=\"Float32\"");
            fprintf(ofp2, " Name=\"Bi-Moment\">\n");
            for (i = 0; i < ne_tr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            for (i = 0; i < ne_fr; ++i) {
                fprintf(ofp2, "     %0.15lf\n", ef_fr[i*7+6]);
            }
            for (i = 0; i < ne_sh; ++i) {
                fprintf(ofp2, "     %0.15lf\n", 0.0);
            }
            fprintf(ofp2, "    </DataArray>\n");
            
            fprintf(ofp2, "   </CellData>\n");
            
        }
        
        fprintf(ofp2, "   <Points>\n");
        
        // Write current joint coordinates to output files
        fprintf(ofp2, "    <DataArray type=\"Float32\" Name=\"\"");
        fprintf(ofp2, " NumberOfComponents=\"3\">\n");
        for (i = 0; i < nj; ++i) {
            if (jcode[i*7] != 0) {
                fprintf(ofp2, "     %0.15lf", x[i*3] + d[jcode[i*7]-1]);
            } else {
                fprintf(ofp2, "     %0.15lf", x[i*3]);
            }
            if (jcode[i*7+1] != 0) {
                fprintf(ofp2, " %0.15lf", x[i*3+1] + d[jcode[i*7+1]-1]);
            } else {
                fprintf(ofp2, " %0.15lf", x[i*3+1]);
            }
            if (jcode[i*7+2] != 0) {
                fprintf(ofp2, " %0.15lf\n", x[i*3+2] + d[jcode[i*7+2]-1]);
            } else {
                fprintf(ofp2, " %0.15lf\n", x[i*3+2]);
            }
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        fprintf(ofp2, "   </Points>\n");
        fprintf(ofp2, "   <Cells>\n");
        
        // Write truss member incidence data to output file
        fprintf(ofp2, "    <DataArray type=\"Int32\" Name=\"connectivity\">\n");
        for (i = 0; i < ne_tr; ++i) {
            fprintf(ofp2, "     %ld %ld\n", minc[i*2] - 1, minc[i*2+1] - 1);
        }
        
        // Write frame member incidence data to output file
        ptr = ne_tr * 2;
        for (i = 0; i < ne_fr; ++i) {
            fprintf(ofp2, "     %ld %ld\n", minc[ptr+i*2] - 1,
                    minc[ptr+i*2+1] - 1);
        }
        
        // Write shell member incidence data to output file
        ptr += ne_fr * 2;
        for (i = 0; i < ne_sh; ++i) {
            fprintf(ofp2, "     %ld %ld %ld\n", minc[ptr+i*3] - 1,
                    minc[ptr+i*3+1] - 1, minc[ptr+i*3+2] - 1);
        }
        
        // Write solid brick member incidence data to output file
        ptr += ne_sh * 3;
        for (i = 0; i < ne_sbr; ++i) {
            fprintf(ofp2, "     %ld %ld %ld %ld %ld %ld %ld %ld\n", minc[ptr+i*8] - 1,
                    minc[ptr+i*8+1] - 1, minc[ptr+i*8+2] - 1, minc[ptr+i*8+3] - 1, minc[ptr+i*8+4] - 1,
                    minc[ptr+i*8+5] - 1, minc[ptr+i*8+6] - 1, minc[ptr+i*8+7] - 1);
        }
        
        // Write fluid brick member incidence data to output file
        ptr += ne_sbr * 8;
        for (i = 0; i < ne_fbr; ++i) {
            fprintf(ofp2, "     %ld %ld %ld %ld %ld %ld %ld %ld\n", minc[ptr+i*8] - 1,
                    minc[ptr+i*8+1] - 1, minc[ptr+i*8+2] - 1, minc[ptr+i*8+3] - 1, minc[ptr+i*8+4] - 1,
                    minc[ptr+i*8+5] - 1, minc[ptr+i*8+6] - 1, minc[ptr+i*8+7] - 1);
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        // Write member incidence offset data to output file
        fprintf(ofp2, "    <DataArray type=\"Int32\" Name=\"offsets\">\n");
        for (i = 0; i < ne_tr; ++i) {
            fprintf(ofp2, "     %ld\n", i * 2 + 2);
        }
        ptr = ne_tr * 2;
        for (i = 0; i < ne_fr; ++i) {
            fprintf(ofp2, "     %ld\n", ptr + i * 2 + 2);
        }
        ptr += ne_fr * 2;
        for (i = 0; i < ne_sh; ++i) {
            fprintf(ofp2, "     %ld\n", ptr + i * 3 + 3);
        }
        ptr += ne_sh * 3;
        for (i = 0; i < ne_sbr; ++i) {
            fprintf(ofp2, "     %ld\n", ptr + i * 8 + 8);
        }
        ptr += ne_sbr * 8;
        for (i = 0; i < ne_fbr; ++i) {
            fprintf(ofp2, "     %ld\n", ptr + i * 8 + 8);
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        // Write cell type data to output file
        fprintf(ofp2, "    <DataArray type=\"UInt32\" Name=\"types\">\n");
        for (i = 0; i < ne_tr; ++i) {
            fprintf(ofp2, "     3\n");
        }
        for (i = 0; i < ne_fr; ++i) {
            fprintf(ofp2, "     3\n");
        }
        for (i = 0; i < ne_sh; ++i) {
            fprintf(ofp2, "     5\n");
        }
        for (i = 0; i < ne_sbr; ++i) {
            fprintf(ofp2, "     12\n");
        }
        for (i = 0; i < ne_fbr; ++i) {
            fprintf(ofp2, "     12\n");
        }
        fprintf(ofp2, "    </DataArray>\n");
        
        // Write back matter for *.vtu file
        fprintf(ofp2, "   </Cells>\n");
        fprintf(ofp2, "  </Piece>\n");
        fprintf(ofp2, " </UnstructuredGrid>\n");
        fprintf(ofp2, "</VTKFile>\n");
        
        k++; // Increment step counter
        
        closeio(1);
        
        // Write *vtu file data to *.pvd file
        fprintf(ofp, "  <DataSet timestep=\"%ld\" part=\"0\"", k - 1);
        fprintf(ofp, " file=\"vtu_files/slave_%ld.vtu\"/>\n", k - 1);
    } while (getc(ifp[1]) != EOF);
    
    // Write back matter for *.pvd file
    fprintf(ofp, " </Collection>\n");
    fprintf(ofp, "</VTKFile>");
    
    return free_all (p2p2l, nl, p2p2d, nd);
}

// This function closes the input and output files
int closeio (int flag)
{
    // Initialize function variables
    int i;
    
    if (flag == 0) {
        // Close the input files
        for (i = 0; i < 5; ++i) {
            if (fclose(ifp[i]) != 0) {
                printf("***ERROR*** Unable to close the input file\n");
            }
        }
        // Close the output file
        if (fclose(ofp) != 0) {
            printf("***ERROR*** Unable to close the input file\n");
        }
    } else {
        // Close the output file
        if (fclose(ofp2) != 0) {
            printf("***ERROR*** Unable to close the output file\n");
        }
        ofp2 = NULL;
    }
    return 0;
}

// This function allocates memory for an array of type long
long * alloc_long (long arraylen)
{
    long *a;
    a = (long *) malloc(arraylen * sizeof(long));
    if (a == NULL) {
        printf("***ERROR*** Unable to allocate memory");
        return NULL;
    }
    return a;
}

// This function allocates memory for an array of type double
double * alloc_dbl (long arraylen)
{
    double *a;
    a = (double *) malloc(arraylen * sizeof(double));
    if (a == NULL) {
        printf("***ERROR*** Unable to allocate memory");
        return NULL;
    }
    return a;
}

// This function frees all allocated memory
int free_all (long **pp2p2l, int nl, double **pp2p2d, int nd)
{
    // Initialize function variables
    int i;
    
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
    
    return closeio(0);
}
