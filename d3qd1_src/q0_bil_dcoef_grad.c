/*
 * q0_bil_dcoef_grad.c
 *
 *  Created on: Sep 27, 2019
 *      Author: ohta
 */

#include "d3qd1_elem.h"

void q0_bil_dcoef_grad_4p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[3];
  double vR[3],vr[3],ro2,ro,aR,i_aR,aW,*vW,rdW,i_ro,vyro,vzro,rodW,rodWm,vp[3][3],rodv,rodvm;
  int i,l,err,p;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_grad_4p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_grad_4p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  CC[8]=0.0;
  for(p=0;p<3;p++){
    dC[p][8]=0.0;
    for(l=0;l<3;l++) vp[p][l]=0.0;
  }
  vp[0][0]=1.0;
  vp[1][1]=1.0;
  vp[2][2]=1.0;

  for(i=0;i<4;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];
    vW=bd->wen[s][i];
    ro2=vR[1]*vR[1]+vR[2]*vR[2];
    ro=sqrt(ro2);
    aR=sqrt(vR[0]*vR[0]+ro2);
    i_aR=1.0/aR;
    aW=vabs_d(vW);
    rdW=vdot_d(vR,vW)*i_aR;
    i_ro=1.0/ro;
    vyro=vR[1]*i_ro;
    vzro=vR[2]*i_ro;
    rodW=vyro*vW[1]+vzro*vW[2];
    rodWm=vyro*vW[1]-vzro*vW[2];
    vr[0]=vR[0];
    vr[1]=ro;
    vr[2]=0.0;
    err=d3hm_qpgf_d1_cv2(&qG,dqG,d2qG,vr,MEPS,&(bd->qd));
    if(err<0){
      printf("q0_bil_dcoef_grad.c, q0_bil_dcoef_grad_4p(), d3hm_qpgf_d1_cv2() error. err=%d. Exit...\n",err);
      exit(0);
    }

    CC[i  ]=qG*aW; // qG
    CC[i+4]=dqG[0]*vW[0]+dqG[1]*(vR[1]*vW[1]+vR[2]*vW[2]); // qH
    CC[8]+=bd->wt_44[i]*rdW*i_aR*i_aR; // F

    for(p=0;p<3;p++){
      rodv =vyro*vp[p][1]+vzro*vp[p][2];
      rodvm=vyro*vp[p][1]-vzro*vp[p][2];
      dC[p][i  ]=-(dqG[0]*vp[p][0]+dqG[1]*(vR[1]*vp[p][1]+vR[2]*vp[p][2]))*aW; // dqG/dv
      dC[p][i+4]=-(d2qG[0]*vW[0]*vp[p][0]+d2qG[1]*(rodW*vp[p][0]+rodv*vW[0])+d2qG[2]*rodW*rodv+dqG[1]*rodWm*rodvm);
      dC[p][8]+=bd->wt_44[i]*i_aR*i_aR*i_aR*(3.0*vR[p]*i_aR*rdW-vW[p]); // dF/dv
    }
  }
  CC[8]*=I_FP;
  for(p=0;p<3;p++) dC[p][8]*=I_FP;
}

void q0_bil_dcoef_grad_9p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[3];
  double vR[3],vr[3],ro2,ro,aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],Np,rdW,i_ro,vyro,vzro,rodW,rodv,rodWm,rodvm,vp[3][3];
  int i,l,ep,err,p;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_grad_9p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_grad_9p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n",creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(p=0;p<3;p++){
    for(l=0;l<9;l++) dC[p][l]=0.0;
    for(l=0;l<3;l++) vp[p][l]=0.0;
  }
  vp[0][0]=1.0;
  vp[1][1]=1.0;
  vp[2][2]=1.0;

  for(ep=0;ep<4;ep++){
    for(i=0;i<9;i++){
      zeta=bd->zt_49[i];
      eta =bd->et_49[i];
      bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      ro2=vR[1]*vR[1]+vR[2]*vR[2];
      ro=sqrt(ro2);
      aR=sqrt(vR[0]*vR[0]+ro2);
      i_aR=1.0/aR;
      aW=vabs_d(vW);
      rdW=vdot_d(vR,vW)*i_aR;
      Np=bil_Nn(ep,zeta,eta);
      i_ro=1.0/ro;
      vyro=vR[1]*i_ro;
      vzro=vR[2]*i_ro;
      rodW=vyro*vW[1]+vzro*vW[2];
      rodWm=vyro*vW[1]-vzro*vW[2];
      vr[0]=vR[0];
      vr[1]=ro;
      vr[2]=0.0;
      err=d3hm_qpgf_d1_cv2(&qG,dqG,d2qG,vr,MEPS,&(bd->qd));
      if(err<0){
        printf("q0_bil_dcoef_grad.c, q0_bil_dcoef_grad_9p(), d3hm_qpgf_d1_cv2() error. err=%d. Exit...\n",err);
        exit(0);
      }

      CC[ep  ]+=bd->wt_49[i]*qG*Np*aW; // qG
      CC[ep+4]+=bd->wt_49[i]*(dqG[0]*vW[0]+dqG[1]*(vR[1]*vW[1]+vR[2]*vW[2]))*Np; // qH

      for(p=0;p<3;p++){
        rodv =vyro*vp[p][1]+vzro*vp[p][2];
        rodvm=vyro*vp[p][1]-vzro*vp[p][2];
        dC[p][ep  ]+=bd->wt_49[i]*(dqG[0]*vp[p][0]+dqG[1]*(vR[1]*vp[p][1]+vR[2]*vp[p][2]))*aW*Np; // dqG/dv
        dC[p][ep+4]+=bd->wt_49[i]*(d2qG[0]*vW[0]*vp[p][0]+d2qG[1]*(rodW*vp[p][0]+rodv*vW[0])+d2qG[2]*rodW*rodv+dqG[1]*rodWm*rodvm)*Np; // dqH/dv
      }
      if(ep==0){
        CC[8]+=bd->wt_49[i]*rdW*i_aR*i_aR;
        for(p=0;p<3;p++) dC[p][8]+=bd->wt_49[i]*i_aR*i_aR*i_aR*(3.0*vR[p]*i_aR*rdW-vW[p]);
      }
    }
    for(p=0;p<3;p++){
      dC[p][ep  ]*=-1.0;
      dC[p][ep+4]*=-1.0;
    }
    if(ep==0){
      CC[8]*=I_FP;
      for(p=0;p<3;p++) dC[p][8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_grad_GL(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[3],tmpg,tmph,tmpdg[3],tmpdh[3];
  double vR[3],vr[3],ro2,ro,aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],rdW,Np,tmpf,tmpdf[3],i_ro,vyro,vzro,rodW,rodv,rodWm,rodvm,vp[3][3];
  int i,j,l,ep,err,p;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_grad_GL(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_grad_GL(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(p=0;p<3;p++){
    for(l=0;l<9;l++) dC[p][l]=0.0;
    for(l=0;l<3;l++) vp[p][l]=0.0;
  }
  vp[0][0]=1.0;
  vp[1][1]=1.0;
  vp[2][2]=1.0;

  for(ep=0;ep<4;ep++){
    for(i=0;i<GLN;i++){
      tmpg =0.0;
      tmph =0.0;
      tmpf =0.0;
      for(p=0;p<3;p++){
        tmpdg[p]=0.0;
        tmpdh[p]=0.0;
        tmpdf[p]=0.0;
      }
      for(j=0;j<GLN;j++){
        zeta=bd->xli[i];
        eta =bd->xli[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        ro2=vR[1]*vR[1]+vR[2]*vR[2];
        ro=sqrt(ro2);
        aR=sqrt(vR[0]*vR[0]+ro2);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        rdW=vdot_d(vR,vW)*i_aR;
        Np=bil_Nn(ep,zeta,eta);
        i_ro=1.0/ro;
        vyro=vR[1]*i_ro;
        vzro=vR[2]*i_ro;
        rodW=vyro*vW[1]+vzro*vW[2];
        rodWm=vyro*vW[1]-vzro*vW[2];
        vr[0]=vR[0];
        vr[1]=ro;
        vr[2]=0.0;
        err=d3hm_qpgf_d1_cv2(&qG,dqG,d2qG,vr,MEPS,&(bd->qd));
        if(err<0){
          printf("q0_bil_dcoef_grad.c, q0_bil_dcoef_grad_GL(), d3hm_qpgf_d1_cv2() error. err=%d. Exit...\n",err);
          exit(0);
        }

        tmpg +=bd->wli[j]*qG*Np*aW;
        tmph +=bd->wli[j]*(dqG[0]*vW[0]+dqG[1]*(vR[1]*vW[1]+vR[2]*vW[2]))*Np;
        for(p=0;p<3;p++){
          rodv =vyro*vp[p][1]+vzro*vp[p][2];
          rodvm=vyro*vp[p][1]-vzro*vp[p][2];
          tmpdg[p]+=bd->wli[j]*(dqG[0]*vp[p][0]+dqG[1]*(vR[1]*vp[p][1]+vR[2]*vp[p][2]))*aW*Np;
          tmpdh[p]+=bd->wli[j]*(d2qG[0]*vW[0]*vp[p][0]+d2qG[1]*(rodW*vp[p][0]+rodv*vW[0])+d2qG[2]*rodW*rodv+dqG[1]*rodWm*rodvm)*Np;
        }
        if(ep==0){
          tmpf +=bd->wli[j]*rdW*i_aR*i_aR;
          for(p=0;p<3;p++) tmpdf[p]+=bd->wli[j]*i_aR*i_aR*i_aR*(3.0*vR[p]*i_aR*rdW-vW[p]);
        }
      }
      CC[ep  ]+=bd->wli[i]*tmpg;
      CC[ep+4]+=bd->wli[i]*tmph;
      for(p=0;p<3;p++){
        dC[p][ep  ]+=bd->wli[i]*tmpdg[p];
        dC[p][ep+4]+=bd->wli[i]*tmpdh[p];
      }
      if(ep==0){
        CC[8] +=bd->wli[i]*tmpf;
        for(p=0;p<3;p++) dC[p][8]+=bd->wli[i]*tmpdf[p];
      }
    }
    for(p=0;p<3;p++){
      dC[p][ep  ]*=-1.0;
      dC[p][ep+4]*=-1.0;
    }
    if(ep==0){
      CC[8]*=I_FP;
      for(p=0;p<3;p++) dC[p][8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_grad_GH(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[3],tmpg,tmph,tmpdg[3],tmpdh[3];
  double vR[3],vr[3],ro2,ro,aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],rdW,Np,tmpf,tmpdf[3],i_ro,vyro,vzro,rodW,rodv,rodWm,rodvm,vp[3][3];
  int i,j,l,ep,err,p;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_grad_GH(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_grad_GH(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
          ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(p=0;p<3;p++){
    for(l=0;l<9;l++) dC[p][l]=0.0;
    for(l=0;l<3;l++) vp[p][l]=0.0;
  }
  vp[0][0]=1.0;
  vp[1][1]=1.0;
  vp[2][2]=1.0;

  for(ep=0;ep<4;ep++){
    for(i=0;i<GHN;i++){
      tmpg =0.0;
      tmph =0.0;
      tmpf =0.0;
      for(p=0;p<3;p++){
        tmpdg[p]=0.0;
        tmpdh[p]=0.0;
        tmpdf[p]=0.0;
      }
      for(j=0;j<GHN;j++){
        zeta=bd->xhi[i];
        eta =bd->xhi[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        ro2=vR[1]*vR[1]+vR[2]*vR[2];
        ro=sqrt(ro2);
        aR=sqrt(vR[0]*vR[0]+ro2);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        rdW=vdot_d(vR,vW)*i_aR;
        Np=bil_Nn(ep,zeta,eta);
        i_ro=1.0/ro;
        vyro=vR[1]*i_ro;
        vzro=vR[2]*i_ro;
        rodW=vyro*vW[1]+vzro*vW[2];
        rodWm=vyro*vW[1]-vzro*vW[2];
        vr[0]=vR[0];
        vr[1]=ro;
        vr[2]=0.0;
        err=d3hm_qpgf_d1_cv2(&qG,dqG,d2qG,vr,MEPS,&(bd->qd));
        if(err<0){
          printf("q0_bil_dcoef_grad.c, q0_bil_dcoef_grad_GH(), d3hm_qpgf_d1_cv2() error. err=%d. Exit...\n",err);
          exit(0);
        }

        tmpg +=bd->whi[j]*qG*Np*aW;
        tmph +=bd->whi[j]*(dqG[0]*vW[0]+dqG[1]*(vR[1]*vW[1]+vR[2]*vW[2]))*Np;

        for(p=0;p<3;p++){
          rodv =vyro*vp[p][1]+vzro*vp[p][2];
          rodvm=vyro*vp[p][1]-vzro*vp[p][2];
          tmpdg[p]+=bd->whi[j]*(dqG[0]*vp[p][0]+dqG[1]*(vR[1]*vp[p][1]+vR[2]*vp[p][2]))*aW*Np;
          tmpdh[p]+=bd->whi[j]*(d2qG[0]*vW[0]*vp[p][0]+d2qG[1]*(rodW*vp[p][0]+rodv*vW[0])+d2qG[2]*rodW*rodv+dqG[1]*rodWm*rodvm)*Np;
        }
        if(ep==0){
          tmpf +=bd->whi[j]*rdW*i_aR*i_aR;
          for(p=0;p<3;p++) tmpdf[p]+=bd->whi[j]*i_aR*i_aR*i_aR*(3.0*vR[p]*i_aR*rdW-vW[p]);
        }
      }
      CC[ep  ]+=bd->whi[i]*tmpg;
      CC[ep+4]+=bd->whi[i]*tmph;
      for(p=0;p<3;p++){
        dC[p][ep  ]+=bd->whi[i]*tmpdg[p];
        dC[p][ep+4]+=bd->whi[i]*tmpdh[p];
      }
      if(ep==0){
        CC[8] +=bd->whi[i]*tmpf;
        for(p=0;p<3;p++) dC[p][8]+=bd->whi[i]*tmpdf[p];
      }
    }
    for(p=0;p<3;p++){
      dC[p][ep  ]*=-1.0;
      dC[p][ep+4]*=-1.0;
    }
    if(ep==0){
      CC[8]*=I_FP;
      for(p=0;p<3;p++) dC[p][8]*=I_FP;
    }
  }
}
