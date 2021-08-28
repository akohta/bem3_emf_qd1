/*
 * d3_qpgf_d1.h
 *
 *  Created on: Mar 17, 2019
 *      Author: ohta
 */

#ifndef SRC_D3_QPGF_D1_H_
#define SRC_D3_QPGF_D1_H_

#include "expint.h"
#include "chankel1_01.h"
#include "mod_besK01.h"
#include "Faddeeva.h"

// limit number
#define S_LIMIT 40

// Fourier domain summation
#define FD_LIMIT 100

// Ewald method
#define EW_LIMIT 100
#define EW_EXPM 3.0

// struct
typedef struct qp_data{
  double k;  // wave number
  double d;  // periodic number ( lattice vector is (d,0,0) )
  double kx; // x component of wave-vector

}QPDT;

int d3hm_qpgf_d1_fd(double complex *qG,double complex *dqG,double *r,double eps,QPDT *qd);
/* Fourier domain summation
    qG : quasi-periodic green function
  *dqG : derivative of quasi-periodic green function
         dqG[0]=d/dx(qG), dqG[1]=d/dy(qG), dqG[2]=d/dz(qG)
     r : (x,y,z)
         r[0]=x, r[1]=y, r[2]=z
   eps : requested relative error
   *qd : pointer of quasi-periodic data (struct QPDT)

   return code >0 : total summation number
               -1 : abnormal termination. reached to summation limit
               -2 : abnormal termination. \beta_n == 0.0 (divergence)
*/

int d3hm_qpgf_d1_ew(double complex *qG,double complex *dqG,double *r,double eps,QPDT *qd);
/* Ewald method
    qG : quasi-periodic green function
  *dqG : derivative of quasi-periodic green function
         dqG[0]=d/dx(qG), dqG[1]=d/dy(qG), dqG[2]=d/dz(qG)
     r : (x,y,z)
         r[0]=x, r[1]=y, r[2]=z
   eps : requested relative error
   *qd : pointer of quasi-periodic data (struct QPDT)

   return code >0 : total summation number
               -1 : abnormal termination. qG1 reached to summation limit
               -2 : abnormal termination. qG2 reached to summation limit
               -3 : abnormal termination. \beta_n == 0.0 (divergence)
*/

int d3hm_qpgf_d1(double complex *qG,double complex *dqG,double *r,double eps,QPDT *qd);
/* auto select
   return code >0 : total summation number
              -11 : Fourier domain summation, abnormal termination. reached to summation limit
              -12 : Fourier domain summation, abnormal termination. \beta_n == 0.0 (divergence)
              -21 : Ewald method, abnormal termination. qG1 reached to summation limit
              -22 : Ewald method, abnormal termination. qG2 reached to summation limit
              -23 : Ewald method, abnormal termination. \beta_n == 0.0 (divergence)
*/
int d3hm_qpgf_d1_qG (double complex  *qG,double *r,double eps,QPDT *qd);
int d3hm_qpgf_d1_dqG(double complex *dqG,double *r,double eps,QPDT *qd);

#endif /* SRC_D3_QPGF_D1_H_ */
