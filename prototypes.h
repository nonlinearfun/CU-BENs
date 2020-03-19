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


/* Define function prototypes - convention: "p" precedes the name of the variable in
   main being pointed to in the function */

/*
model.c
*/

/* This function reads in member incidences and joint constraints, checks for errors in
   the joint constraint input and employs METIS functions to decompose the domain by
   DOF-weighted joints */
int struc (long *pjcode, long *pminc, int *pwrpres, long *pjnt);

/* This function handles all acoustic fluid structure details for the monolithic coupling with the structural / solid mesh */
void fsi (long *pmcode,  long *pjcode, long *pminc, long *pelface, long *pfsiinc);

void renumfsi (long *pminc, long *pmcode, long *pjcode, long *pelface, long *pfsiinc);

/* This function generates the adjacency structure of the graph associated with the
   inputted finite element mesh. */
int graph (long *pjnt, long *pxadj, long *pjcode, int *pwrpres, long *pminc, long *pjinc,
    long *pjincloc);

/* This function provides the optimum node-numbering scheme, with respect to minimization
   of the bandwidth of the stiffness array; it follows the procedure outlined by Collins
   (1973) */
int optnum (long *pjnt, long *padjncy, long *pxadj, int *pflag);

/* This function updates all previously defined joint-dependent variables so that they
   reflect the optimum node-numbering scheme determined from optnum() */
int updatenum (long *pjcode, long *pminc, int *pwrpres,  long *pjnt);

/* This function generates joint code, jcode, by assigning integers in sequence, by
   columns, to all nonzero elements of jcode from 1 to neq and generates the member code,
   mcode, by transferring, via minc, columns of jcode into columns of mcode; it follows
   procedure in Chapter 6 of Bathe and Wilson as well as Chapter 6 of Holzer */
void codes (long *pmcode, long *pjcode, long *pminc, int *pwrpres);

/* This function reads in: (1) the magnitude (mag) and the direction (dir) of the
   concentrated force, and the joint (jt) on which it is applied AND / OR (2) the
   magnitude (mag) and the direction (dir) of the distributed force, and the frame
   element (fr) on which it is applied */
int load (double *pq, double *pefFE_ref, double *px, double *pllength, double *poffset,
		  int *posflag, double *pc1, double *pc2, double *pc3, long *pjnt, long *pmcode, long *pjcode, 
		  long *pminc, double *ptinpt, double *ppinpt, double *pdinpt, double *ppdisp, double *pum, double *pvm, double *pam);

// This function determines kht using mcode, and determines maxa from kht
int skylin (long *pmaxa, long *pmcode, long *plss, long *pjcode, long *pkht, long *ppmot);

/*
truss.c
*/

/* This function reads and echoes truss element properties, as well as computing element
   lengths and direction cosines */
void prop_tr (double *px, double *pemod, double *pcarea, double *pdens, double *pllength,
              double *pyield, double *pc1, double *pc2, double *pc3, long *pminc);

/* This function computes the contribution to the generalized stiffness matrix from the
   truss elements and stores it as an array */
void stiff_tr (double *pss, double *pemod, double *pcarea, double *plength,
    double *pdefllen_ip, double *pyield, double *pc1_ip, double *pc2_ip, double *pc3_ip,
    double *pef_ip, long *pmaxa, long *pmcode);

// This function adds contribution of material nonlinearity to element stiffness matrix
void stiffm_tr (double *pk_tr, double *pef_ip, double *pPy, long n);

/* This function computes the contribution to the generalized internal force vector from
   the truss elements and truss element internal force vectors */
void forces_tr (double *pf_temp, double *pef_i, double *pd, double *pemod,
    double *pcarea, double *pllength, double *pdefllen_i, double *pyield, double *pc1_i,
    double *pc2_i, double *pc3_i, long *pmcode);

/* This function computes the contribution to the generalized mass matrix from the
   truss elements and stores it as an array */
void mass_tr (double *psm, double *pcarea, double *pllength, double *pdens, double *px, 
         long *pminc, long *pmcode, double *pjac);

/*
frame.c
*/

/* This function reads and echoes frame element properties, as well as computing element
   lengths and direction cosines */
void prop_fr (double *px, double *pxfr, double *pemod, double *pgmod, double *pdens, double *poffset,
              int *posflag, double *pauxpt, double *pcarea, double *pllength, double *pistrong,
              double *piweak, double *pipolar, double *piwarp, double *pyield, double *pzstrong,
              double *pzweak, double *pc1, double *pc2, double *pc3, int *pmendrel, long *pminc);

/* This function computes the contribution to the generalized stiffness matrix from the
   frame elements and stores it as an array */
void stiff_fr (double *pss, double *pemod, double *pgmod, double *pcarea,
    double *poffset, int *posflag, double *pllength, double *pdefllen_ip,
    double *pistrong, double *piweak, double *pipolar, double *piwarp, int *pyldflag,
    double *pyield, double *pzstrong, double *pzweak, double *pc1_ip, double *pc2_ip,
    double *pc3_ip, double *pef_ip, double *pefFE_ip, int *pmendrel, long *pmaxa,
    long *pmcode);

// This function assigns non-zero elements of linear-elastic element stiffness matrix
void stiffe_fr (double *pk_fr, double *pemod, double *pgmod, double *pcarea,
    double *plength, double *pistrong, double *piweak, double *pipolar, double *piwarp,
    long n);

// This function adds contribution of geometric nonlinearity to element stiffness matrix
void stiffg_fr (double *pk_fr, double *peftot_ip, double *pdefllen_ip, double *pcarea,
    double *pipolar, long n);

// This function adds contribution of material nonlinearity to element stiffness matrix
void stiffm_fr (double *pk_fr, double *peftot_ip, int *pyldflag, double *pyield,
    double *pcarea, double *pzstrong, double *pzweak, long n);

/* This function modifies the frame element stiffness matrix to account for member end
   releases associated with bending degrees of freedom; this algorithm is similar to that
   implemented for the computation of the plastic reduction matrix */
void release (double *pk_fr, int *pmendrel, int n);

/* This function computes the contribution to the generalized internal force vector from
   the frame elements and frame element internal force vectors */
int forces_fr (double *pf_temp, double *pef_ip, double *pef_i, double *pefFE_ref,
    double *pefFE_ip, double *pefFE_i, int *pyldflag, double *pdd, double *pemod,
    double *pgmod, double *pcarea, double *poffset, int *posflag, double *pllength,
    double *pdefllen_ip, double *pistrong, double *piweak, double *pipolar,
    double *piwarp, double *pyield, double *pzstrong, double *pzweak, double *pc1_ip,
    double *pc2_ip, double *pc3_ip, double *pc1_i, double *pc2_i, double *pc3_i,
    int *pmendrel, long *pmcode, double *pdlpf, int *pitecnt);

/* This function computes the contribution to the generalized mass matrix from the
    frame elements and stores it as an array */
void mass_fr (double *psm, double *pcarea, double *pllength, double *pistrong, double *piweak, 
              double *pipolar, double *piwarp, double *pdens, int *posflag, double *poffset, 
              double *px, double *pxfr, long *pminc, long *pmcode, double *pjac);

/* This function, referred to as the "method of false position" or "regula falsi"
   (McGuire, Gallagher, and Ziemian, p. 356), computes the scalar multiplier to be
   applied to the current increment in force, such that the yield function to equal one,
   within a specified tolerance */
double regula_falsi (double *pp, double *pdp, double *pmy, double *pdmy, double *pmz,
    double *pdmz);

/* This function checks for the presence of unloading at hinges and returns appropriate
   flag for analysis to continue */
int unload (double *pphi, double *pp, double *pmy, double *pmz, double *pPy,
    double *pMpy, double *pMpz, double *pk_fr, double *pdd, long n);

/*
shell.c
*/

/* This function reads and echoes shell element properties, as well as computing element
   side-lengths, face-areas, and direction cosines */
void prop_sh (double *px, double *pemod, double *pnu, double *pxlocal, double *pthick,  double *pdens,
    double *pfarea, double *pslength, double *pyield, double *pc1, double *pc2,
    double *pc3, long *pminc);

/* This function computes the contribution to the generalized stiffness matrix from the
   shell elements and stores it as an array */
void stiff_sh (double *pss, double *pemod, double *pnu, double *px_temp, double *pxlocal,
    double *pthick, double *pfarea, double *pdeffarea_ip, double *pslength,
    double *pdefslen_ip, double *pyield, double *pc1_ip, double *pc2_ip, double *pc3_ip,
    double *pef_ip, double *pd_temp, double *pchi_temp, double *pefN_temp,
    double *pefM_temp, long *pmaxa, long *pminc, long *pmcode);

// This function assigns non-zero elements of linear-elastic element stiffness matrix
void stiffe_sh (double *pk_sh, double *pemod, double *pnu, double *pxlocal,
    double *pthick, double *pfarea, double *pslength, long ptr, long n);

// This function evaluates the linear-elastic membrane stiffness matrix
void stiffe_m_sh (double *pke_m_sh, double *pemod, double *pnu, double *pxlocal,
    double *pthick, double *pfarea, long ptr, long n);

/* This function employs an algorithm which explicitly evaluates the linear-elastic plate
   bending stiffness matrix (i.e. no numerical integration is performed) by taking
   advantage of the properties of the matrices involved in its evaluation. Details of the
   algorithm are found in Batoz (1982). */
void stiffe_b_sh (double *pke_b_sh, double *pemod, double *pnu, double *pxlocal,
    double *pthick, double *pfarea, double *pslength, long ptr, long n);

// This function adds contribution of geometric nonlinearity to element stiffness matrix
void stiffg_sh (double *pk_sh, double *pemod, double *pnu, double *pxlocal,
    double *pthick, double *pdeffarea_ip, double *pdm, long ptr, long n);

// This function assigns non-zero elements of elasto-plastic element stiffness matrix
void stiffm_sh (double *pk_sh, double *pemod, double *pnu, double *pxlocal,
    double *pthick, double *pdefffarea_ip, double *pdefslen_ip, double *pyield,
    double *pefN_temp, double *pefM_temp, double *pchi_temp, int yv, double *pNo,
    double *palpha, double *pMe, double *pNbar, double *pMbar, double *pMNbar,
    double *pq_fact, double *pr_fact, double *ps_fact, int *ph_fact, long ptr,
    long n);

// This function evaluates the elasto-plastic membrane stiffness matrix
void stiffm_m_sh (double *pkm_m_sh, double *pxlocal, double *pthick,
    double *pdeffarea_ip, double *pC, double *pfn, double *pj_fact, double *pk_fact,
    double *pB_fact, double *pdf_da, double *pda_dchi, long n);

/* This function evaluates the elasto-plastic plate bending stiffness matrix, employing
   a modified, but similar approach as in the stiff_b function */
void stiffm_b_sh (double *pkm_b_sh, double *pxlocal, double *pthick,
    double *pdeffarea_ip, double *pdefslen_ip, double *pC, double *pfm, double *pj_fact,
    double *pk_fact, double *pB_fact, double *pdf_da, double *pda_dchi, long n);

/* This function evaluates the coupled elasto-plastic stiffness matrix, employing a
   modified, but similar approach as in the stiff_b and stiffm_b functions */
void stiffm_mb_sh (double *pkm_mb_sh, double *pxlocal, double *pthick,
    double *pdeffarea_ip, double *pdefslen_ip, double *pC, double *pfn, double *pfm,
    double *pj_fact, double *pk_fact, double *pB_fact, double *pdf_da, double *pda_dchi,
    long n);

/* This function computes the contribution to the generalized mass matrix from the
   shell elements and stores it as an array */
void mass_sh (double *psm, double *pcarea, double *pdens, double *pthick,
              double *pfarea, double *pslength, 
              double *px, long *pminc, long *pmcode, double *pjac);

/* This function computes the contribution to the generalized internal force vector from
   the shell elements and shell element internal force vectors */
int forces_sh (double *pf_temp, double *pef_ip, double *pef_i, double *pefN_temp,
    double *pefM_temp, double *pdd, double *pd_temp, double *pchi_temp, double *px_temp,
    double *px_ip, double *pemod, double *pnu, double *pxlocal, double *pthick,
    double *pfarea, double *pdeffarea_ip, double *pslength, double *pdefslen_ip,
    double *pyield, double *pc1_ip, double *pc2_ip, double *pc3_ip, double *pc1_i,
    double *pc2_i, double *pc3_i, long *pminc, long *pmcode, long *pjcode);

// This function computes the local membrane coordinates for each element
void mem_coord (double *px_local, long n, long i, long j, long k, long l, double *px,
    double *pc1, double *pc2, double *pc3, long ptr);

/* This function computes the total (elastic plus plastic) strain and curvature
   increments */
void strn_curv (double *pstrn, double *pcurv, double *pddm, double *pddb,
    double *pxlocal, double *pdeffarea_ip, double *pdefslen_ip, long n);

/*
 brick.c
 */

void prop_br (double *pemod, double *pnu, double *pyield, double *pdens, double *pfdens, double *pbmod,
			  double *pfarea);

void stiff_br (double *pss, double *px, double *pemod, double *pnu, long *pminc, 
			   long *pmcode, long *pjcode, double *pJinv, double *pjac);

void mass_br (double *psm, double *pdens, double *px, long *pminc, long *pmcode, double *pjac);

void jacob (double *px, long *pminc, long *el, int *rval, int *sval, int *tval, double *pjac);

/*
 fsi.c
 */

void prop_fsi (double *px, double *pemod, double *pnu, double *pdens, double *pfdens, double *pbmod,
			   double *pfarea, 
			   double *pslength, double *pyield, long *pminc, long *pelface, long *pfsiinc, 
			   double *pnnorm, double *ptarea, double *pss, double *pss_fsi, double *psd_fsi, double *pabspt,
			   double *pnorpt, long *pmcode, long *pjcode, double *pL);

void stiff_fsi (long *pminc, long *pmcode, long *pjcode, double *pnnorm, double *ptarea, double *pfarea, double *pthick, 
				double *pdeffarea, double *pslength, double *pdefslen, double *pL, double *pA, double *pss, double *pss_fsi, 
				double *px, double *pxlocal, double *pemod, double *pnu, double *pJinv, double *pjac, double *pyield, 
				double *pc1, double *pc2, double *pc3, double *pef, double *pd, double *pchi, double *pefN, double *pefM, long *pmaxa);

void mass_fsi (long *pminc, long *pmcode, long *pjcode, double *pnnorm, double *ptarea, double *pcarea, double *pfarea, 
               double *pthick, double *pslength, double *pL, double *pLT, double *psm, double *psm_fsi, double *px, 
               double *pdens, double *pfdens, double *pJinv, double *pjac);

void L_br (long *pminc, long *pmcode, long *pjcode, long *pjcode_fsi, double *pnnorm, double *ptarea,
		   double *pL, double *pA, double *pG);

void load_fsi (long *pjcode_fsi, double *ptinpt, double *ppinpt, double *ppresinpt, double *paccinpt, double *pfdens,
			   double *pum, double *pvm, double *pam);

void q_fsi (long *pjcode_fsi, double *pqdyn, double *ptstps, double *papload, double *ppres, double *pacc, 
			double *pL, double *pA, double *pLp, double *pAu, double *pfdens, double *ptinpt, double *ppinpt, double *ppresinpt, double *paccinpt);

/*
solve.c
*/

/* This function determines the solution within the increment using the compact Gaussian
   elimination strategy (Holzer pp. 290, 296, 307) and based on subroutine COLSOL -
   "active column solution or the skyline (or column) reduction method"
   (Bathe 1976, p. 257) */
int solve (long *pjcode_fsi, double *pss, double *pss_fsi, double *psm, double *psm_fsi, double *psd_fsi, double *pr, double *pdd, long *pmaxa, double *pssd, int *pdet,
		   double *pum, double *pvm, double *pam, double *puc, double *pvc, double *pac, double *pqdyn, double *ptstps, 
		   double *pKeff, double *pReff, double *pMeff, int *pAp, int *pAi, double *pAx, double alpham, double alphaf, int *pipiv, int fact, double ddt, double *ppdisp, long *pkht, int *piter, int *pii, int *pij, int tstp);

/* This function performs LDL^t factorization of the stiffness matrix when using SLVFLAG == 0*/
int skyfact (long *pmaxa, double *pss_temp, double *pssd, double *pdd, int fact, int *pdet);

/*This function solves for x in Ax=b*/
int skysolve (long *pmaxa, double *pss, double *pssd, double *pdd, int fact, int *pdet);

/*This function solves for b in Ax=b*/
int skymult (long *pmaxa, double *pss_temp, double *pdd);

/*This function partitions matrix A into interior dofs and boundary dofs for solving systems with prescribed displacement boundary conditions*/
int matpart (long *pmaxa, long *pkht, double *pss, double *pqtot, double *puc, int *pii, int *pij);

/*
arc.c
*/

/* This function reads in the MSAL parameters and checks that the prescribed initial
   displacement does not correspond to a fixed degree of freedom */
int msal (double *pdk, long *pdkdof, long *pjnt, long *pjcode);

/* This function computes the roots of the quadratic equation, determines the correct
   root, and updates nodal displacements and load proportionality factor */
int quad (double *pa, double *pb, double *pc, double *pdt, double *pdp, double *pddr,
    double *pddq, double *pdd, double *pdlpf, double *plpft);

/*
misc.c
*/

/* This function transforms the element stiffness matrix from local into global
   coordinate system */
void transform (double *pk, double *pT_ip, double *pK, int n);

/* This function performs convergence check; assume that convergence is reached,
   i.e. convchk = 0, until proven otherwise */
int test (double *pd_temp, double *pdd, double *pf_temp, double *pfp, double *pqtot,
    double *pf_ip, double *pintener1, int *pconvchk, double *ptoldisp, double *ptolfor,
    double *ptolener);

// This function updates element geometry
void updatc (double *px_temp, double *px_ip, double *pxfr_temp, double *pdd,
    double *pdefllen_i, double *pdeffarea_i, double *pdefslen_i, double *poffset,
    int *posflag, double *pauxpt, double *pc1_i, double *pc2_i, double *pc3_i,
    long *pminc, long *pjcode);

// This function computes the dot product of two vectors
double dot (double *pa, double *pb, int n);

// This function computes the cross product of two vectors
void cross (double *pa, double *pb, double *pc, int flag);

// This function computes the inverse of a matrix through Gauss-Jordan elimination
void inverse (double *pGT_k_G, int n);

// This function outputs variables of interest to appropriate files
void output (double *plpf, int *pitecnt, double *pd, double *pef, int flag);

// This function closes the input and output files
int closeio (int flag);

// This function periodically backs up information for nonlinear dynamic analysis in the event of power outrages or hardware failures
void checkPoint(long tstep, long lss, double *puc, double *pvc, double *pac, double *pss, double *psm, double *pd,
                double *pf, double *pef, double *px, double *pc1, double *pc2, double *pc3, double *pdefllen,
                double *pllength, double *pefFE, double *pxfr, int *pyldflag, double *pdeffarea, double *pdefslen,
                double *pchi, double *pefN, double *pefM);

// The function restores the last saved data before power outrages or hardware failures
void restartStep(long lss, double *puc, double *pvc, double *pac, double *pss, double *psm, double *pd,
                 double *pf, double *pef, double *px, double *pc1, double *pc2, double *pc3, double *pdefllen,
                 double *pllength, double *pefFE, double *pxfr, int *pyldflag, double *pdeffarea, double *pdefslen,
                 double *pchi, double *pefN, double *pefM);

/*
memory.c
*/

// This function allocates memory for an array of type int
int * alloc_int (long arraylen);

// This function allocates memory for an array of type long
long * alloc_long (long arraylen);

// This function allocates memory for an array of type double
double * alloc_dbl (long arraylen);

// This function frees all allocated memory
int free_all (int **pp2p2i, int ni, long **pp2p2l, int nl, double **pp2p2d, int nd,
    int flag);
