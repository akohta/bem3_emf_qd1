/*
 * bem3_emf_qd1.h
 *
 *  Created on: Aug 21, 2019
 *      Author: ohta
 */

#ifndef BEM3_EMF_QD1_H_
#define BEM3_EMF_QD1_H_

#include "d3qd1_const.h"
#include "d3qd1_elem.h"

#define PREC_DEF_BE 2 // type settings for d3qd1_bv_solver.
#define PREC_DEF_FN 0 // type settings of field analysis functions in force_FN().

// -- d3qd1_setup.c --
void read_dqd1(int argc,char **argv,DQD1 *dq1); // read datafile for solver
void print_dqd1(DQD1 *qd);                      // print data
void print_dqd1_mksa(DQD1 *qd);                 // print data in MKSA system of units
void initialize_dqd1(DQD1 *qd);                 // memory allocation and initialize data for solver
void finalize_dqd1(DQD1 *qd);                   // memory free 
int domain_id(double *rt,DQD1 *qd);             // return domain id of point rt, return the main domain id if on boundary. for non periodic model.
int q0_domain_id(double *rt,DQD1 *qd);          // return domain id of point rt, return the main domain id if on boundary. for periodic model.
int q0_domain_id_l(double *rt,DQD1 *qd);        // ret periodic number l
void dat_write(char *filename,DQD1 *qd);        // output analysis result to binary file with specified filename
void dat_read(char *filename,DQD1 *qd);         // read datafile outputed by dat_write()


// -- d3qd1_solve_bieq.c --
void solve_bieq(DQD1 *md);                      // solve boundary integral equation  


// -- d3qd1_field.c --
// analysis modified electromagnetic potential by using boundary integral equations 
int mEMP_s(double complex *U,double *rt,int type,DQD1 *qd); // scattered or internal field
int mEMP_t(double complex *U,double *rt,int type,DQD1 *qd); // total field ( add incident field to scattered field )
int mEMP_i(double complex *U,double *rt,int type,DQD1 *qd); // incident field
// outputs
// U[0]=Ux, U[1]=Uy, U[2]=Uz ( vector potential ), U[3]=phi ( scalar potential ), return object id.
// inputs
// rt[0]=x, rt[1]=y, rt[2]=z, type, pointer of DQD1.
// type=0:4-point GL, type=1:9-pont or 7-point GL, type=2:GLN-point GL, type=3:GHN-point GL, type=4:DE integration.
// the higher the type numbers are, the lower error is, the slower calculation speed is. it is usually setted to 0 or 1. 
// the type must be set higher number when analyse near the boundary ( in the order of a element ).

// analysis electromagnetic field by using derivative boundary integral equations
int EH_mEMP_s(double complex *E,double complex *H,double *rt,int type,DQD1 *qd); // scattered or internal field 
int EH_mEMP_t(double complex *E,double complex *H,double *rt,int type,DQD1 *qd); // total field ( add incident field to scattered field )
int EH_mEMP_i(double complex *E,double complex *H,double *rt,int type,DQD1 *qd); // incident field
// outputs
// E[0]=Ex, E[1]=Ey, E[2]=Ez, H[0]=Hx, H[1]=Hy, H[2]=Hz, return object id.
// inputs
// rt[0]=x, rt[1]=y, rt[2]=z, type, pointer of DQD1.
// type is the same as mEMP_*()
// these fields are calculated from modified electromagnetic potential using derivative boundary integral equations.
// these are slower than EH_*() functions. the error is smaller than EH_*() functions in far-field.

// analysis electromagnetic field by using boundary integral equations
int EH_s(double complex *E,double complex *H,double *rt,int type,DQD1 *qd); 
int EH_t(double complex *E,double complex *H,double *rt,int type,DQD1 *qd);
int EH_i(double complex *E,double complex *H,double *rt,int type,DQD1 *qd);
// outputs and inputs are the same as EH_mEMP_*() functions.

// analysis electromagnetic field on the boundary 
void EH_s_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD1 *qd);
void EH_t_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD1 *qd);
void EH_i_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DQD1 *qd);
// outputs
// E[0]=Ex, E[1]=Ey, E[2]=Ez, H[0]=Hx, H[1]=Hy, H[2]=Hz.
// inputs
// did:domain id, t:element id defined in each domain, zeta_t,eta_t:point parameter on the element surface ( main domain side coordinate ), type, pointer of DQD1.
// type is the same as EH_*().

// -- d3qd1_force.c --
int force_FN(double *F,double *N,double *rc,int type,DQD1 *qd);
// outputs
// F[0]=Fx, F[1]=Fy, F[2]=Fz, N[0]=Nx, N[1]=Ny, N[2]=Nz.
// inputs
// rc[0]=x, rc[1]=y, rc[2]=z ( center of rotation ), type, pointer of DQD1.
// type=0:4-point GL, type!=0:9-point or 7-point GL.

#endif
