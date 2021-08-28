#ifndef SRC_D3_QPGF_D1_CS_H_
#define SRC_D3_QPGF_D1_CS_H_

#include "d3_qpgf_d1.h"

// quasi-periodic Green's functions of Helmholtz equation. 1-dimensional(x-axis) periodicity. cylindrical coordinate.
// definition : x=x, y=rho*cos(phi), z=rho*sin(phi), rho=sqrt(y^2+z^2)

int d3hm_qpgf_d1_cv1(double complex *qG,double complex *dqG,double *cr,double eps,QPDT *qd);
/* quasi-periodic Green's function and its first derivative
               qG : quasi-periodic Green's function.
             *dqG : first derivative of quasi-periodic Green's function.
                    dqG[0]=dqG/dx, dqG[1]=1/rho*dqG/drho, dqG[2]=dqG/dphi=0. 
               cr : argument in sylindrical coordinate system
                    r[0]=x, r[1]=rho, r[2]=phi ( phi is not needed )
              eps : requested relative error.
              *qd : pointer of quasi-periodic data (struct QPDT)

   return code >0 : total summation number
              -11 : Fourier domain summation, abnormal termination. reached to summation limit
              -12 : Fourier domain summation, abnormal termination. \beta_n == 0.0
              -21 : Ewald method, abnormal termination. qG1 reached to summation limit
              -22 : Ewald method, abnormal termination. qG2 reached to summation limit
              -23 : Ewald method, abnormal termination. \beta_n == 0.0 (divergence)
*/
int d3hm_qpgf_d1_cv1_qG(double complex *qG,double *cr,double eps,QPDT *qd);
int d3hm_qpgf_d1_cv1_dqG(double complex *dqG,double *cr,double eps,QPDT *qd);

int d3hm_qpgf_d1_cv2(double complex *qG,double complex *dqG,double complex *d2qG,double *cr,double eps,QPDT *qd);
/* quasi-periodic Green's function and its first and second derivative
               qG : quasi-periodic Green's function.
             *dqG : first derivative of quasi-periodic green's function.
                    dqG[0]=dqG/dx, dqG[1]=1/rho*dqG/drho, dqG[2]=dqG/dphi=0.
            *d2qG : second derivative of quasi-periodic green's function
                    d2qG[0]=d^2qG/dx^2, d2qG[1]=d^2qG/drhodx, d2qG[2]=d^2qG/drho^2.
               cr : argument in sylindrical coordinate system
                    r[0]=x, r[1]=rho, r[2]=phi ( phi is not needed )
              eps : requested relative error.
              *qd : pointer of quasi-periodic data (struct QPDT)

   return code >0 : total summation number
              -11 : Fourier domain summation, abnormal termination. reached to summation limit
              -12 : Fourier domain summation, abnormal termination. \beta_n == 0.0
              -21 : Ewald method, abnormal termination. qG1 reached to summation limit
              -22 : Ewald method, abnormal termination. qG2 reached to summation limit
              -23 : Ewald method, abnormal termination. \beta_n == 0.0
*/

#endif
