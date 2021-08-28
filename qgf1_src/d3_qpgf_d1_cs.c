#include "d3_qpgf_d1_cs.h"

int d3hm_qpgf_d1_cv1(double complex *qG,double complex *dqG,double *cr,double eps,QPDT *qd)
{
  int d3hm_qpgf_d1_cv1_fd(double complex *qG,double complex *dqG,double *cr,double eps,QPDT *qd);
  int d3hm_qpgf_d1_cv1_ew(double complex *qG,double complex *dqG,double *cr,double eps,QPDT *qd);

  double kd,k2,bN2,aN,ap,bp2,am,bm2,FN,Fp,Fm,ro,abn;
  int N,err;

  N=(int)round(qd->kx*qd->d/(2.0*M_PI));

  kd=2.0*M_PI/qd->d;
  k2=qd->k*qd->k;

  ap=kd*(double)(N+S_LIMIT)-qd->kx;
  bp2=k2-ap*ap;
  am=kd*(double)(N-S_LIMIT)-qd->kx;
  bm2=k2-am*am;

  if(bp2>0.0 || bm2>0.0){ // Ewald method
    err=d3hm_qpgf_d1_cv1_ew(qG,dqG,cr,eps,qd);
    if(err<0) err-=20;
  }

  ro=cr[1];
  aN=kd*(double)N-qd->kx;
  bN2=k2-aN*aN;

  abn=sqrt(fabs(bN2));
  FN=0.5*M_PI*cabs(besj0(abn*ro)+I*besy0(abn*ro));
  abn=sqrt(fabs(bp2));
  Fp=fabs(bessk0(abn*ro));
  abn=sqrt(fabs(bm2));
  Fm=fabs(bessk0(abn*ro));

  if(Fp < eps*FN && Fm < eps*FN){ // Fourier
    err=d3hm_qpgf_d1_cv1_fd(qG,dqG,cr,eps,qd);
    if(err<0) err-=10;
  }
  else { // Ewald
    err=d3hm_qpgf_d1_cv1_ew(qG,dqG,cr,eps,qd);
    if(err<0) err-=20;
  }

  return err;
}

int d3hm_qpgf_d1_cv2(double complex *qG,double complex *dqG,double complex *d2qG,double *cr,double eps,QPDT *qd)
{
  int d3hm_qpgf_d1_cv2_fd(double complex *qG,double complex *dqG,double complex *d2qG,double *cr,double eps,QPDT *qd);
  int d3hm_qpgf_d1_cv2_ew(double complex *qG,double complex *dqG,double complex *d2qG,double *cr,double eps,QPDT *qd);

  double kd,k2,bN2,aN,ap,bp2,am,bm2,FN,Fp,Fm,ro,abn;
  int N,err;

  N=(int)round(qd->kx*qd->d/(2.0*M_PI));

  kd=2.0*M_PI/qd->d;
  k2=qd->k*qd->k;

  ap=kd*(double)(N+S_LIMIT)-qd->kx;
  bp2=k2-ap*ap;
  am=kd*(double)(N-S_LIMIT)-qd->kx;
  bm2=k2-am*am;

  if(bp2>0.0 || bm2>0.0){ // Ewald method
    err=d3hm_qpgf_d1_cv2_ew(qG,dqG,d2qG,cr,eps,qd);
    if(err<0) err-=20;
  }

  ro=cr[1];
  aN=kd*(double)N-qd->kx;
  bN2=k2-aN*aN;

  abn=sqrt(fabs(bN2));
  FN=0.5*M_PI*cabs(besj0(abn*ro)+I*besy0(abn*ro));
  abn=sqrt(fabs(bp2));
  Fp=fabs(bessk0(abn*ro));
  abn=sqrt(fabs(bm2));
  Fm=fabs(bessk0(abn*ro));

  if(Fp < eps*FN && Fm < eps*FN){ // Fourier
    err=d3hm_qpgf_d1_cv2_fd(qG,dqG,d2qG,cr,eps,qd);
    if(err<0) err-=10;
  }
  else { // Ewald
    err=d3hm_qpgf_d1_cv2_ew(qG,dqG,d2qG,cr,eps,qd);
    if(err<0) err-=20;
  }

  return err;
}

int d3hm_qpgf_d1_cv1_qG(double complex *qG,double *cr,double eps,QPDT *qd)
{
  int d3hm_qpgf_d1_cv1_fd_qG(double complex *qG,double *cr,double eps,QPDT *qd);
  int d3hm_qpgf_d1_cv1_ew_qG(double complex *qG,double *cr,double eps,QPDT *qd);

  double kd,k2,bN2,aN,ap,bp2,am,bm2,FN,Fp,Fm,ro,abn;
  int N,err;

  N=(int)round(qd->kx*qd->d/(2.0*M_PI));

  kd=2.0*M_PI/qd->d;
  k2=qd->k*qd->k;

  ap=kd*(double)(N+S_LIMIT)-qd->kx;
  bp2=k2-ap*ap;
  am=kd*(double)(N-S_LIMIT)-qd->kx;
  bm2=k2-am*am;

  if(bp2>0.0 || bm2>0.0){ // Ewald method
    err=d3hm_qpgf_d1_cv1_ew_qG(qG,cr,eps,qd);
    if(err<0) err-=20;
  }

  ro=cr[1];
  aN=kd*(double)N-qd->kx;
  bN2=k2-aN*aN;

  abn=sqrt(fabs(bN2));
  FN=0.5*M_PI*cabs(besj0(abn*ro)+I*besy0(abn*ro));
  abn=sqrt(fabs(bp2));
  Fp=fabs(bessk0(abn*ro));
  abn=sqrt(fabs(bm2));
  Fm=fabs(bessk0(abn*ro));

  if(Fp < eps*FN && Fm < eps*FN){ // Fourier
    err=d3hm_qpgf_d1_cv1_fd_qG(qG,cr,eps,qd);
    if(err<0) err-=10;
  }
  else { // Ewald
    err=d3hm_qpgf_d1_cv1_ew_qG(qG,cr,eps,qd);
    if(err<0) err-=20;
  }

  return err;
}

int d3hm_qpgf_d1_cv1_dqG(double complex *dqG,double *cr,double eps,QPDT *qd)
{
  int d3hm_qpgf_d1_cv1_fd_dqG(double complex *dqG,double *cr,double eps,QPDT *qd);
  int d3hm_qpgf_d1_cv1_ew_dqG(double complex *dqG,double *cr,double eps,QPDT *qd);

  double kd,k2,bN2,aN,ap,bp2,am,bm2,FN,Fp,Fm,ro,abn;
  int N,err;

  N=(int)round(qd->kx*qd->d/(2.0*M_PI));

  kd=2.0*M_PI/qd->d;
  k2=qd->k*qd->k;

  ap=kd*(double)(N+S_LIMIT)-qd->kx;
  bp2=k2-ap*ap;
  am=kd*(double)(N-S_LIMIT)-qd->kx;
  bm2=k2-am*am;

  if(bp2>0.0 || bm2>0.0){ // Ewald method
    err=d3hm_qpgf_d1_cv1_ew_dqG(dqG,cr,eps,qd);
    if(err<0) err-=20;
  }

  ro=cr[1];
  aN=kd*(double)N-qd->kx;
  bN2=k2-aN*aN;

  abn=sqrt(fabs(bN2));
  FN=0.5*M_PI*cabs(besj0(abn*ro)+I*besy0(abn*ro));
  abn=sqrt(fabs(bp2));
  Fp=fabs(bessk0(abn*ro));
  abn=sqrt(fabs(bm2));
  Fm=fabs(bessk0(abn*ro));

  if(Fp < eps*FN && Fm < eps*FN){ // Fourier
    err=d3hm_qpgf_d1_cv1_fd_dqG(dqG,cr,eps,qd);
    if(err<0) err-=10;
  }
  else { // Ewald
    err=d3hm_qpgf_d1_cv1_ew_dqG(dqG,cr,eps,qd);
    if(err<0) err-=20;
  }

  return err;  
}

/////////////////////////////////////////////////////////////
int d3hm_qpgf_d1_cv1_fd(double complex *qG,double complex *dqG,double *cr,double eps,QPDT *qd)
{ // Fourier domain summation
  int d3hm_qpgf_d1_cv1_fd_sd(double complex *qG,double complex *dqG,double *cr,double eps,QPDT *qd);

  double complex c;
  double rt[3],tmp;
  int err,l;

  rt[0]=cr[0];
  rt[1]=cr[1];
  rt[2]=cr[2];

  l=(int)floor(cr[0]/qd->d+0.5);
  rt[0]-=(double)l*qd->d;

  err=d3hm_qpgf_d1_cv1_fd_sd(qG,dqG,rt,eps,qd);
  if(l!=0){
    tmp=(double)l*qd->kx*qd->d;
    c=cos(tmp)-I*sin(tmp);
    *qG*=c;
    dqG[0]*=c;
    dqG[1]*=c;
  }

  return err;
}

int d3hm_qpgf_d1_cv1_fd_sd(double complex *qG,double complex *dqG,double *cr,double eps,QPDT *qd)
{

  double complex ekx,ten,i_ten,eiaxN,eiaxp,eiaxm,ch,tc,tf,oG,odG0,odG1;
  double arg,ro,k2,kd,kdx,bn2,aN,anp,anm,bn,bnr,ct,st;
  int i,N,flg[3],n,sig;

  ro=cr[1];
  k2=qd->k*qd->k;

  N=(int)round(qd->kx*qd->d/(2.0*M_PI));
  ch=M_PI*I/2.0;

  arg=-qd->kx*cr[0];
  ekx=cos(arg)+I*sin(arg);

  kd=2.0*M_PI/qd->d;
  kdx=kd*cr[0];
  ct=cos(kdx);
  st=sin(kdx);
  ten=ct+I*st;
  i_ten=ct-I*st;

  for(i=0;i<3;i++) flg[i]=0.0;

  n=N;
  aN=kd*(double)n-qd->kx;
  arg=kdx*(double)n;
  eiaxN=(cos(arg)+I*sin(arg))*ekx;
  bn2=k2-aN*aN;
  if(bn2>0.0){
    bn=sqrt(bn2);
    bnr=bn*ro;
    tc=eiaxN*ch;
    tf=tc*(besj0(bnr)+I*besy0(bnr));
    *qG=tf;
    dqG[0]=aN*tf;
    dqG[1]=-bn*tc*(besj1(bnr)+I*besy1(bnr));
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2);
    bnr=bn*ro;
    tf=eiaxN*bessk0(bnr);
    *qG=tf;
    dqG[0]=aN*tf;
    dqG[1]=-bn*eiaxN*bessk1(bnr);
  }
  else return -2; // beta_n=0.0;

  eiaxp=eiaxN;
  eiaxm=eiaxN;
  anp=aN;
  anm=aN;
  for(n=0;n<FD_LIMIT;n++){
    oG=*qG;
    odG0=dqG[0];
    odG1=dqG[1];
    sig=0;

    eiaxp*=ten;
    eiaxm*=i_ten;
    anp+=kd;
    anm-=kd;

    // +
    bn2=k2-anp*anp;
    if(bn2>0.0){
      bn=sqrt(bn2);
      bnr=bn*ro;
      tc=eiaxp*ch;
      tf=tc*(besj0(bnr)+I*besy0(bnr));
      *qG+=tf;
      dqG[0]+=anp*tf;
      dqG[1]+=-bn*tc*(besj1(bnr)+I*besy1(bnr));
      sig+=1;
    }
    else if(bn2<0.0){
      bn=sqrt(-bn2);
      bnr=bn*ro;
      tf=eiaxp*bessk0(bnr);
      *qG+=tf;
      dqG[0]+=anp*tf;
      dqG[1]+=-bn*eiaxp*bessk1(bnr);
    }
    else return -2; // beta_n=0.0;

    // -
    bn2=k2-anm*anm;
    if(bn2>0.0){
      bn=sqrt(bn2);
      bnr=bn*ro;
      tc=eiaxm*ch;
      tf=tc*(besj0(bnr)+I*besy0(bnr));
      *qG+=tf;
      dqG[0]+=anm*tf;
      dqG[1]+=-bn*tc*(besj1(bnr)+I*besy1(bnr));
      sig+=1;
    }
    else if(bn2<0.0){
      bn=sqrt(-bn2);
      bnr=bn*ro;
      tf=eiaxm*bessk0(bnr);
      *qG+=tf;
      dqG[0]+=anm*tf;
      dqG[1]+=-bn*eiaxm*bessk1(bnr);
    }
    else return -2; // beta_n=0.0;

    if(sig==0){
      if(flg[0]==0){
        if(  2.0*fabs(creal(oG)-creal(*qG)) < eps*(fabs(creal(oG))+fabs(creal(*qG)))
          && 2.0*fabs(cimag(oG)-cimag(*qG)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG))) ) flg[0]=1;
      }
      if(flg[1]==0){
        if(  2.0*fabs(creal(odG0)-creal(dqG[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG[0])))
          && 2.0*fabs(cimag(odG0)-cimag(dqG[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG[0]))) ) flg[1]=1;
      }
      if(flg[2]==0){
        if(  2.0*fabs(creal(odG1)-creal(dqG[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG[1])))
          && 2.0*fabs(cimag(odG1)-cimag(dqG[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG[1]))) ) flg[2]=1;
      }
    }

    if(flg[0]!=0 && flg[1]!=0 && flg[2]!=0) break;
  }

  if(n==FD_LIMIT) return -1;

  tc=1.0/(2.0*M_PI*qd->d);
  *qG*=tc;
  dqG[0]*=I*tc;

  tc/=ro;
  dqG[1]*=tc;
  dqG[2]=0.0;

  return 2*n+1;
}

int d3hm_qpgf_d1_cv1_ew(double complex *qG,double complex *dqG,double *cr,double eps,QPDT *qd)
{ // Ewald method
  int d3hm_ewald_qG1_cv1(double complex *qG1,double complex *dqG1,double *cr,double eps,double veps,QPDT *qd);
  int d3hm_ewald_qG2_cv1(double complex *qG1,double complex *dqG1,double *cr,double eps,double veps,QPDT *qd);

  double complex tqG,tdqG[3];
  double veps,r02;
  int err1,err2;

  veps=qd->d/(2.0*sqrt(M_PI));

  // veps
  r02=cr[0]*cr[0]+cr[1]*cr[1];
  if( (-r02/(4.0*veps*veps)+qd->k*qd->k*veps*veps) > EW_EXPM)
    veps=sqrt(EW_EXPM+sqrt(qd->k*qd->k*r02+EW_EXPM*EW_EXPM))/(sqrt(2.0)*qd->k);

  // qG1
  err1=d3hm_ewald_qG1_cv1(&tqG,tdqG,cr,eps,veps,qd);
  *qG=tqG;
  dqG[0]=tdqG[0];
  dqG[1]=tdqG[1];
  if(err1==-1) return -1;

  // qG2
  err2=d3hm_ewald_qG2_cv1(&tqG,tdqG,cr,eps,veps,qd);
  *qG+=tqG;
  dqG[0]+=tdqG[0];
  dqG[1]+=tdqG[1];
  if(err2==-1) return -2;
  if(err2==-2) return -3; // beta_n==0.0

  dqG[2]=0.0;

  return err1+err2;
}

int d3hm_ewald_qG1_cv1(double complex *qG1,double complex *dqG1,double *cr,double eps,double veps,QPDT *qd)
{
  double complex ekd,i_ekd,eLkd,tw,elkdp,elkdm,tce,tcf,oG,odG0,odG1;
  double kd,ro2,rl2,rl,i_rl,xLd,arg,ke,ke2,i2e,mi2e2,xldp,xldm,st,ct,td,cpe;
  int L,l,flg[3];

  kd=qd->kx*qd->d;
  st=sin(kd);
  ct=cos(kd);
  ekd=ct+I*st;
  i_ekd=ct-I*st;

  ro2=cr[1]*cr[1];
  ke=qd->k*veps;
  ke2=ke*ke;
  i2e=1.0/(2.0*veps);
  mi2e2=-i2e*i2e;
  cpe=1.0/(sqrt(M_PI)*veps);

  if(fabs(cr[0])<qd->d*sqrt(eps)) L=1; // avoid symmetric cancellation
  else L=0;
  flg[0]=0;  flg[1]=0;  flg[2]=0;

  arg=kd*(double)L;
  eLkd=cos(arg)+I*sin(arg);
  xLd=cr[0]+(double)L*qd->d;
  rl2=xLd*xLd+ro2;
  rl=sqrt(rl2);
  i_rl=1.0/rl;
  tw=Faddeeva_w(ke+I*rl*i2e,eps);
  tce=eLkd*i_rl*exp(rl2*mi2e2+ke2);
  tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
  *qG1=tce*creal(tw);
  dqG1[0]=xLd*tcf;
  dqG1[1]=tcf;

  elkdp=eLkd;
  elkdm=eLkd;
  xldp=xLd;
  xldm=xLd;
  for(l=0;l<EW_LIMIT;l++){
    oG=*qG1;
    odG0=dqG1[0];
    odG1=dqG1[1];
    elkdp*=ekd;
    elkdm*=i_ekd;
    xldp+=qd->d;
    xldm-=qd->d;

    // +
    rl2=xldp*xldp+ro2;
    rl=sqrt(rl2);
    i_rl=1.0/rl;
    tw=Faddeeva_w(ke+I*rl*i2e,eps);
    tce=elkdp*i_rl*exp(rl2*mi2e2+ke2);
    tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
    *qG1+=tce*creal(tw);
    dqG1[0]+=xldp*tcf;
    dqG1[1]+=tcf;

    //-
    rl2=xldm*xldm+ro2;
    rl=sqrt(rl2);
    i_rl=1.0/rl;
    tw=Faddeeva_w(ke+I*rl*i2e,eps);
    tce=elkdm*i_rl*exp(rl2*mi2e2+ke2);
    tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
    *qG1+=tce*creal(tw);
    dqG1[0]+=xldm*tcf;
    dqG1[1]+=tcf;

    if(flg[0]==0){
      if(  2.0*fabs(creal(oG)-creal(*qG1)) < eps*(fabs(creal(oG))+fabs(creal(*qG1)))
        && 2.0*fabs(cimag(oG)-cimag(*qG1)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG1))) ) flg[0]=1;
    }
    if(flg[1]==0){
      if(  2.0*fabs(creal(odG0)-creal(dqG1[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG1[0])))
        && 2.0*fabs(cimag(odG0)-cimag(dqG1[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG1[0]))) ) flg[1]=1;
    }
    if(flg[2]==0){
      if(  2.0*fabs(creal(odG1)-creal(dqG1[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG1[1])))
        && 2.0*fabs(cimag(odG1)-cimag(dqG1[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG1[1]))) ) flg[2]=1;
    }

    if(flg[0]!=0 && flg[1]!=0 && flg[2]!=0) break;
  }

  td=1.0/(4.0*M_PI);
  *qG1*=td;
  dqG1[0]*=-td;
  dqG1[1]*=-td;

  if(l!=EW_LIMIT) return 2*l+1;
  else return -1;
}

int d3hm_ewald_qG2_cv1(double complex *qG2,double complex *dqG2,double *cr,double eps,double veps,QPDT *qd)
{
  int sum_expintv(double complex *ret0,double complex *ret1,double ro2,double veps,double bn2,double eps);

  double complex ten,i_ten,eiaxN,tc0,tc1,eiaxp,eiaxm,oG,odG0,odG1;
  double ro2,kd,k2,ve2,arg,ct,st,aN,bn2,anp,anm;
  int N,n,err,flg[3],sig;

  flg[0]=0;  flg[1]=0;  flg[2]=0;
  N=(int)round(qd->kx*qd->d/(2.0*M_PI));

  ro2=cr[1]*cr[1];
  kd=2.0*M_PI/qd->d;
  k2=qd->k*qd->k;
  ve2=veps*veps;

  arg=kd*cr[0];
  st=sin(arg);  ct=cos(arg);
  ten=ct+I*st;
  i_ten=ct-I*st;

  n=N;
  aN=(double)n*kd-qd->kx;
  arg=aN*cr[0];
  eiaxN=cos(arg)+I*sin(arg);
  bn2=k2-aN*aN; if(bn2==0.0) return -2;
  err=sum_expintv(&tc0,&tc1,ro2,veps,bn2,eps);  if(err<0) return -1;
  *qG2=eiaxN*tc0;
  dqG2[0]=aN*eiaxN*tc0;
  dqG2[1]=eiaxN*tc1;

  eiaxp=eiaxN;
  eiaxm=eiaxN;
  anp=aN;
  anm=aN;
  for(n=0;n<EW_LIMIT;n++){
    oG=*qG2;
    odG0=dqG2[0];
    odG1=dqG2[1];
    eiaxp*=ten;
    eiaxm*=i_ten;
    anp+=kd;
    anm-=kd;
    sig=0;

    // +
    bn2=k2-anp*anp;  if(bn2==0.0) return -2;
    if(bn2<0.0) sig++;
    err=sum_expintv(&tc0,&tc1,ro2,veps,bn2,eps);  if(err<0) return -1;
    *qG2+=eiaxp*tc0;
    dqG2[0]+=anp*eiaxp*tc0;
    dqG2[1]+=eiaxp*tc1;

    // -
    bn2=k2-anm*anm;  if(bn2==0.0) return -2;
    if(bn2<0.0) sig++;
    err=sum_expintv(&tc0,&tc1,ro2,veps,bn2,eps);  if(err<0) return -1;
    *qG2+=eiaxm*tc0;
    dqG2[0]+=anm*eiaxm*tc0;
    dqG2[1]+=eiaxm*tc1;

    if(sig>0){
      if(flg[0]==0){
        if(  2.0*fabs(creal(oG)-creal(*qG2)) < eps*(fabs(creal(oG))+fabs(creal(*qG2)))
          && 2.0*fabs(cimag(oG)-cimag(*qG2)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG2))) ) flg[0]=1;
      }
      if(flg[1]==0){
        if(  2.0*fabs(creal(odG0)-creal(dqG2[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG2[0])))
          && 2.0*fabs(cimag(odG0)-cimag(dqG2[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG2[0]))) ) flg[1]=1;
      }
      if(flg[2]==0){
        if(  2.0*fabs(creal(odG1)-creal(dqG2[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG2[1])))
          && 2.0*fabs(cimag(odG1)-cimag(dqG2[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG2[1]))) ) flg[2]=1;
      }
    }

    if(flg[0]!=0 && flg[1]!=0 && flg[2]!=0) break;
  }

  tc0=1.0/(4.0*M_PI*qd->d);
  tc1=tc0/(2.0*ve2);
  *qG2*=tc0;
  dqG2[0]*=I*tc0;
  dqG2[1]*=tc1;

  if(n!=EW_LIMIT) return 2*n+1;
  else return -1;
}

int sum_expintv(double complex *ret0,double complex *ret1,double ro2,double veps,double bn2,double eps)
{
  double complex E1,E0,tG0,tG1,tG0o,tG1o;
  double sig,ipm,rpm,z,im,dp,rpm1,ipm1;
  int m,flg[2];

  flg[0]=0;  flg[1]=0;
  dp=ro2/(4.0*veps*veps);
  z=-bn2*veps*veps;

  sig=1.0;
  ipm=1.0;
  rpm=1.0;
  if(bn2>0.0) E1=I*M_PI-expint_i(-z);
  else if(bn2<0.0) E1=expint(1,z);
  else return -3;
  tG0=sig*ipm*rpm*E1;
  tG1=0.0;
  for(m=1;m<EW_LIMIT;m++){
    tG0o=tG0;
    tG1o=tG1;
    sig*=-1.0;
    im=1.0/(double)m;
    ipm1=ipm;
    ipm*=im;
    rpm1=rpm;
    rpm*=dp;
    E0=E1;
    E1=im*(exp(-z)-z*E0);
    tG0+=sig*ipm*rpm*E1;
    tG1+=sig*ipm1*rpm1*E1;

    if(bn2>0.0){
      if(flg[0]==0){
        if(  2.0*fabs(creal(tG0o)-creal(tG0)) < eps*(fabs(creal(tG0o))+fabs(creal(tG0)))
          && 2.0*fabs(cimag(tG0o)-cimag(tG0)) < eps*(fabs(cimag(tG0o))+fabs(cimag(tG0))) ) flg[0]=1;
      }
      if(flg[1]==0){
        if(  2.0*fabs(creal(tG1o)-creal(tG1)) < eps*(fabs(creal(tG1o))+fabs(creal(tG1)))
          && 2.0*fabs(cimag(tG1o)-cimag(tG1)) < eps*(fabs(cimag(tG1o))+fabs(cimag(tG1))) ) flg[1]=1;
      }
    }
    else {
      if(flg[0]==0){
        if(  2.0*fabs(creal(tG0o)-creal(tG0)) < eps*(fabs(creal(tG0o))+fabs(creal(tG0))) ) flg[0]=1;
      }
      if(flg[1]==0){
        if(  2.0*fabs(creal(tG1o)-creal(tG1)) < eps*(fabs(creal(tG1o))+fabs(creal(tG1))) ) flg[1]=1;
      }
    }
    if(flg[0]!=0 && flg[1]!=0) break;
  }

  *ret0=tG0;
  *ret1=tG1;

  if(m==EW_LIMIT) return -1;
  else return 2*m-1;
}

int d3hm_qpgf_d1_cv2_fd(double complex *qG,double complex *dqG,double complex *d2qG,double *cr,double eps,QPDT *qd)
{
  int d3hm_qpgf_d1_cv2_fd_sd(double complex *qG,double complex *dqG,double complex *d2qG,double *cr,double eps,QPDT *qd);

  double complex c;
  double rt[3],tmp;
  int err,l;

  rt[0]=cr[0];
  rt[1]=cr[1];
  rt[2]=cr[2];

  l=(int)floor(cr[0]/qd->d+0.5);
  rt[0]-=(double)l*qd->d;

  err=d3hm_qpgf_d1_cv2_fd_sd(qG,dqG,d2qG,rt,eps,qd);
  if(l!=0){
    tmp=(double)l*qd->kx*qd->d;
    c=cos(tmp)-I*sin(tmp);
    *qG*=c;
    dqG[0]*=c;
    dqG[1]*=c;

    d2qG[0]*=c;
    d2qG[1]*=c;
    d2qG[2]*=c;

  }

  return err;
}

int d3hm_qpgf_d1_cv2_fd_sd(double complex *qG,double complex *dqG,double complex *d2qG,double *cr,double eps,QPDT *qd)
{

  double complex ekx,ten,i_ten,eiaxN,eiaxp,eiaxm,ch,tc,tf0,tf1,tf2,oG,odG0,odG1,fh0,fh1,fh2,od2G0,od2G1,od2G2;
  double arg,ro,k2,kd,kdx,bn2,aN,anp,anm,bn,bnr,ct,st;
  int i,N,flg[6],n,sig;

  ro=cr[1];
  k2=qd->k*qd->k;

  N=(int)round(qd->kx*qd->d/(2.0*M_PI));
  ch=M_PI*I/2.0;

  arg=-qd->kx*cr[0];
  ekx=cos(arg)+I*sin(arg);

  kd=2.0*M_PI/qd->d;
  kdx=kd*cr[0];
  ct=cos(kdx);
  st=sin(kdx);
  ten=ct+I*st;
  i_ten=ct-I*st;

  for(i=0;i<6;i++) flg[i]=0.0;

  n=N;
  aN=kd*(double)n-qd->kx;
  arg=kdx*(double)n;
  eiaxN=(cos(arg)+I*sin(arg))*ekx;
  bn2=k2-aN*aN;
  if(bn2>0.0){
    bn=sqrt(bn2);
    bnr=bn*ro;
    tc=eiaxN*ch;
    fh0=besj0(bnr)+I*besy0(bnr);
    fh1=besj1(bnr)+I*besy1(bnr);
    fh2=fh0-fh1/bnr;
    tf0= tc*fh0;
    tf1=-tc*fh1;
    tf2=-tc*fh2;
    *qG=tf0;
    dqG[0]=aN*tf0;
    dqG[1]=bn*tf1;
    d2qG[0]=aN*aN*tf0;
    d2qG[1]=aN*bn*tf1;
    d2qG[2]=bn*bn*tf2;
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2);
    bnr=bn*ro;
    tc=eiaxN;
    fh0=bessk0(bnr);
    fh1=bessk1(bnr);
    fh2=-fh0-fh1/bnr;
    tf0= tc*fh0;
    tf1=-tc*fh1;
    tf2=-tc*fh2;
    *qG=tf0;
    dqG[0]=aN*tf0;
    dqG[1]=bn*tf1;
    d2qG[0]=aN*aN*tf0;
    d2qG[1]=aN*bn*tf1;
    d2qG[2]=bn*bn*tf2;
  }
  else return -2; // beta_n=0.0;

  eiaxp=eiaxN;
  eiaxm=eiaxN;
  anp=aN;
  anm=aN;
  for(n=0;n<FD_LIMIT;n++){
    oG=*qG;
    odG0=dqG[0];
    odG1=dqG[1];
    od2G0=d2qG[0];
    od2G1=d2qG[1];
    od2G2=d2qG[2];
    sig=0;

    eiaxp*=ten;
    eiaxm*=i_ten;
    anp+=kd;
    anm-=kd;

    // +
    bn2=k2-anp*anp;
    if(bn2>0.0){
      bn=sqrt(bn2);
      bnr=bn*ro;
      tc=eiaxp*ch;
      fh0=besj0(bnr)+I*besy0(bnr);
      fh1=besj1(bnr)+I*besy1(bnr);
      fh2=fh0-fh1/bnr;
      tf0= tc*fh0;
      tf1=-tc*fh1;
      tf2=-tc*fh2;
      *qG+=tf0;
      dqG[0] +=anp*tf0;
      dqG[1] +=bn*tf1;
      d2qG[0]+=anp*anp*tf0;
      d2qG[1]+=anp*bn*tf1;
      d2qG[2]+=bn*bn*tf2;
      sig+=1;
    }
    else if(bn2<0.0){
      bn=sqrt(-bn2);
      bnr=bn*ro;
      tc=eiaxp;
      fh0=bessk0(bnr);
      fh1=bessk1(bnr);
      fh2=-fh0-fh1/bnr;
      tf0= tc*fh0;
      tf1=-tc*fh1;
      tf2=-tc*fh2;
      *qG+=tf0;
      dqG[0]+=anp*tf0;
      dqG[1]+=bn*tf1;
      d2qG[0]+=anp*anp*tf0;
      d2qG[1]+=anp*bn*tf1;
      d2qG[2]+=bn*bn*tf2;
    }
    else return -2; // beta_n=0.0;


    // -
    bn2=k2-anm*anm;
    if(bn2>0.0){
      bn=sqrt(bn2);
      bnr=bn*ro;
      tc=eiaxm*ch;
      fh0=besj0(bnr)+I*besy0(bnr);
      fh1=besj1(bnr)+I*besy1(bnr);
      fh2=fh0-fh1/bnr;
      tf0= tc*fh0;
      tf1=-tc*fh1;
      tf2=-tc*fh2;
      *qG+=tf0;
      dqG[0] +=anm*tf0;
      dqG[1] +=bn*tf1;
      d2qG[0]+=anm*anm*tf0;
      d2qG[1]+=anm*bn*tf1;
      d2qG[2]+=bn*bn*tf2;
      sig+=1;
    }
    else if(bn2<0.0){
      bn=sqrt(-bn2);
      bnr=bn*ro;
      tc=eiaxm;
      fh0=bessk0(bnr);
      fh1=bessk1(bnr);
      fh2=-fh0-fh1/bnr;
      tf0= tc*fh0;
      tf1=-tc*fh1;
      tf2=-tc*fh2;
      *qG+=tf0;
      dqG[0]+=anm*tf0;
      dqG[1]+=bn*tf1;
      d2qG[0]+=anm*anm*tf0;
      d2qG[1]+=anm*bn*tf1;
      d2qG[2]+=bn*bn*tf2;
    }
    else return -2; // beta_n=0.0;


    if(sig==0){
      if(flg[0]==0){
        if(  2.0*fabs(creal(oG)-creal(*qG)) < eps*(fabs(creal(oG))+fabs(creal(*qG)))
          && 2.0*fabs(cimag(oG)-cimag(*qG)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG))) ) flg[0]=1;
      }
      if(flg[1]==0){
        if(  2.0*fabs(creal(odG0)-creal(dqG[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG[0])))
          && 2.0*fabs(cimag(odG0)-cimag(dqG[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG[0]))) ) flg[1]=1;
      }
      if(flg[2]==0){
        if(  2.0*fabs(creal(odG1)-creal(dqG[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG[1])))
          && 2.0*fabs(cimag(odG1)-cimag(dqG[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG[1]))) ) flg[2]=1;
      }
      if(flg[3]==0){
        if(  2.0*fabs(creal(od2G0)-creal(d2qG[0])) < eps*(fabs(creal(od2G0))+fabs(creal(d2qG[0])))
          && 2.0*fabs(cimag(od2G0)-cimag(d2qG[0])) < eps*(fabs(cimag(od2G0))+fabs(cimag(d2qG[0]))) ) flg[3]=1;
      }
      if(flg[4]==0){
        if(  2.0*fabs(creal(od2G1)-creal(d2qG[1])) < eps*(fabs(creal(od2G1))+fabs(creal(d2qG[1])))
          && 2.0*fabs(cimag(od2G1)-cimag(d2qG[1])) < eps*(fabs(cimag(od2G1))+fabs(cimag(d2qG[1]))) ) flg[4]=1;
      }
      if(flg[5]==0){
        if(  2.0*fabs(creal(od2G2)-creal(d2qG[2])) < eps*(fabs(creal(od2G2))+fabs(creal(d2qG[2])))
          && 2.0*fabs(cimag(od2G2)-cimag(d2qG[2])) < eps*(fabs(cimag(od2G2))+fabs(cimag(d2qG[2]))) ) flg[5]=1;
      }
    }

    if(flg[0]!=0 && flg[1]!=0 && flg[2]!=0 && flg[3]!=0 && flg[4]!=0 && flg[5]!=0) break;
  }

  if(n==FD_LIMIT) return -1;

  tc=1.0/(2.0*M_PI*qd->d);
  *qG*=tc;
  dqG[0]*=I*tc;

  dqG[1]*=tc/ro;
  dqG[2]=0.0;

  d2qG[0]*=-tc;
  d2qG[1]*=I*tc;
  d2qG[2]*=tc;

  return 2*n+1;
}

int d3hm_qpgf_d1_cv2_ew(double complex *qG,double complex *dqG,double complex *d2qG,double *cr,double eps,QPDT *qd)
{
  int d3hm_ewald_qG1_cv2(double complex *qG1,double complex *dqG1,double complex *d2qG1,double *cr,double eps,double veps,QPDT *qd);
  int d3hm_ewald_qG2_cv2(double complex *qG2,double complex *dqG2,double complex *d2qG2,double *cr,double eps,double veps,QPDT *qd);

  double complex tqG,tdqG[3],td2qG[3];
  double veps,r02;
  int err1,err2;

  veps=qd->d/(2.0*sqrt(M_PI));

  // veps
  r02=cr[0]*cr[0]+cr[1]*cr[1];
  if( (-r02/(4.0*veps*veps)+qd->k*qd->k*veps*veps) > EW_EXPM)
    veps=sqrt(EW_EXPM+sqrt(qd->k*qd->k*r02+EW_EXPM*EW_EXPM))/(sqrt(2.0)*qd->k);

  // qG1
  err1=d3hm_ewald_qG1_cv2(&tqG,tdqG,td2qG,cr,eps,veps,qd);
  *qG=tqG;
  dqG[0]=tdqG[0];
  dqG[1]=tdqG[1];
  d2qG[0]=td2qG[0];
  d2qG[1]=td2qG[1];
  d2qG[2]=td2qG[2];
  if(err1==-1) return -1;

  // qG2
  err2=d3hm_ewald_qG2_cv2(&tqG,tdqG,td2qG,cr,eps,veps,qd);
  *qG+=tqG;
  dqG[0]+=tdqG[0];
  dqG[1]+=tdqG[1];
  d2qG[0]+=td2qG[0];
  d2qG[1]+=td2qG[1];
  d2qG[2]+=td2qG[2];
  if(err2==-1) return -2;
  if(err2==-2) return -3; // beta_n==0.0

  dqG[2]=0.0;
  return err1+err2;
}

int d3hm_ewald_qG1_cv2(double complex *qG1,double complex *dqG1,double complex *d2qG1,double *cr,double eps,double veps,QPDT *qd)
{
  double complex ekd,i_ekd,eLkd,tw,elkdp,elkdm,tce,tce2,tcf,tcf2,oG,odG0,odG1,od2G0,od2G1,od2G2;
  double kd,ro,ro2,rl2,rl,i_rl,xLd,arg,ke,ke2,i2e,mi2e2,xldp,xldm,st,ct,td,cpe;
  int L,l,flg[6];

  kd=qd->kx*qd->d;
  st=sin(kd);
  ct=cos(kd);
  ekd=ct+I*st;
  i_ekd=ct-I*st;

  ro2=cr[1]*cr[1];
  ro=sqrt(ro2);
  ke=qd->k*veps;
  ke2=ke*ke;
  i2e=1.0/(2.0*veps);
  mi2e2=-i2e*i2e;
  cpe=1.0/(sqrt(M_PI)*veps);

  if(fabs(cr[0])<qd->d*sqrt(eps)) L=1; // avoid symmetric cancellation
  else L=0;
  flg[0]=0;  flg[1]=0;  flg[2]=0;
  flg[3]=0; flg[4]=0; flg[5]=0;

  arg=kd*(double)L;
  eLkd=cos(arg)+I*sin(arg);
  xLd=cr[0]+(double)L*qd->d;
  rl2=xLd*xLd+ro2;
  rl=sqrt(rl2);
  i_rl=1.0/rl;
  tw=Faddeeva_w(ke+I*rl*i2e,eps);
  tce=eLkd*i_rl*exp(rl2*mi2e2+ke2);
  tce2=i_rl*i_rl*i_rl*tce;
  tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
  tcf2=tce2*((3.0-qd->k*qd->k*rl2)*i_rl*creal(tw)-3.0*qd->k*cimag(tw)+cpe*(3.0-2.0*rl2*mi2e2));
  *qG1=tce*creal(tw);
  dqG1[0]=xLd*tcf;
  dqG1[1]=tcf;
  d2qG1[0]=xLd*xLd*tcf2;
  d2qG1[1]=xLd*tcf2;
  d2qG1[2]=tcf2;

  elkdp=eLkd;
  elkdm=eLkd;
  xldp=xLd;
  xldm=xLd;
  for(l=0;l<EW_LIMIT;l++){
    oG=*qG1;
    odG0=dqG1[0];
    odG1=dqG1[1];
    od2G0=d2qG1[0];
    od2G1=d2qG1[1];
    od2G2=d2qG1[2];

    elkdp*=ekd;
    elkdm*=i_ekd;
    xldp+=qd->d;
    xldm-=qd->d;

    // +
    rl2=xldp*xldp+ro2;
    rl=sqrt(rl2);
    i_rl=1.0/rl;
    tw=Faddeeva_w(ke+I*rl*i2e,eps);
    tce=elkdp*i_rl*exp(rl2*mi2e2+ke2);
    tce2=i_rl*i_rl*i_rl*tce;
    tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
    tcf2=tce2*((3.0-qd->k*qd->k*rl2)*i_rl*creal(tw)-3.0*qd->k*cimag(tw)+cpe*(3.0-2.0*rl2*mi2e2));
    *qG1+=tce*creal(tw);
    dqG1[0]+=xldp*tcf;
    dqG1[1]+=tcf;
    d2qG1[0]+=xldp*xldp*tcf2;
    d2qG1[1]+=xldp*tcf2;
    d2qG1[2]+=tcf2;

    //-
    rl2=xldm*xldm+ro2;
    rl=sqrt(rl2);
    i_rl=1.0/rl;
    tw=Faddeeva_w(ke+I*rl*i2e,eps);
    tce=elkdm*i_rl*exp(rl2*mi2e2+ke2);
    tce2=i_rl*i_rl*i_rl*tce;
    tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
    tcf2=tce2*((3.0-qd->k*qd->k*rl2)*i_rl*creal(tw)-3.0*qd->k*cimag(tw)+cpe*(3.0-2.0*rl2*mi2e2));
    *qG1+=tce*creal(tw);
    dqG1[0]+=xldm*tcf;
    dqG1[1]+=tcf;
    d2qG1[0]+=xldm*xldm*tcf2;
    d2qG1[1]+=xldm*tcf2;
    d2qG1[2]+=tcf2;

    if(flg[0]==0){
      if(  2.0*fabs(creal(oG)-creal(*qG1)) < eps*(fabs(creal(oG))+fabs(creal(*qG1)))
        && 2.0*fabs(cimag(oG)-cimag(*qG1)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG1))) ) flg[0]=1;
    }
    if(flg[1]==0){
      if(  2.0*fabs(creal(odG0)-creal(dqG1[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG1[0])))
        && 2.0*fabs(cimag(odG0)-cimag(dqG1[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG1[0]))) ) flg[1]=1;
    }
    if(flg[2]==0){
      if(  2.0*fabs(creal(odG1)-creal(dqG1[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG1[1])))
        && 2.0*fabs(cimag(odG1)-cimag(dqG1[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG1[1]))) ) flg[2]=1;
    }
    if(flg[3]==0){
      if(  2.0*fabs(creal(od2G0)-creal(d2qG1[0])) < eps*(fabs(creal(od2G0))+fabs(creal(d2qG1[0])))
        && 2.0*fabs(cimag(od2G0)-cimag(d2qG1[0])) < eps*(fabs(cimag(od2G0))+fabs(cimag(d2qG1[0]))) ) flg[3]=1;
    }
    if(flg[4]==0){
      if(  2.0*fabs(creal(od2G1)-creal(d2qG1[1])) < eps*(fabs(creal(od2G1))+fabs(creal(d2qG1[1])))
        && 2.0*fabs(cimag(od2G1)-cimag(d2qG1[1])) < eps*(fabs(cimag(od2G1))+fabs(cimag(d2qG1[1]))) ) flg[4]=1;
    }
    if(flg[5]==0){
      if(  2.0*fabs(creal(od2G2)-creal(d2qG1[2])) < eps*(fabs(creal(od2G2))+fabs(creal(d2qG1[2])))
        && 2.0*fabs(cimag(od2G2)-cimag(d2qG1[2])) < eps*(fabs(cimag(od2G2))+fabs(cimag(d2qG1[2]))) ) flg[5]=1;
    }

    if(flg[0]!=0 && flg[1]!=0 && flg[2]!=0 && flg[3]!=0 && flg[4]!=0 && flg[5]!=0) break;
  }

  td=1.0/(4.0*M_PI);
  *qG1*=td;
  dqG1[0]*=-td;
  dqG1[1]*=-td;
  d2qG1[0]=td*d2qG1[0]+dqG1[1];
  d2qG1[1]*=ro*td;
  d2qG1[2]=ro2*td*d2qG1[2]+dqG1[1];

  if(l!=EW_LIMIT) return 2*l+1;
  else return -1;
}

int d3hm_ewald_qG2_cv2(double complex *qG2,double complex *dqG2,double complex *d2qG2,double *cr,double eps,double veps,QPDT *qd)
{
  int sum_expintv2(double complex *ret0,double complex *ret1,double complex *ret2,double ro2,double veps,double bn2,double eps);

  double complex ten,i_ten,eiaxN,tc0,tc1,tc2,eiaxp,eiaxm,oG,odG0,odG1,od2G0,od2G1,od2G2;
  double ro2,kd,k2,ve2,arg,ct,st,aN,bn2,anp,anm,r_2e,i_ve;
  int N,n,err,flg[6],sig;

  flg[0]=0;  flg[1]=0;  flg[2]=0;
  flg[3]=0; flg[4]=0; flg[5]=0;
  N=(int)round(qd->kx*qd->d/(2.0*M_PI));

  ro2=cr[1]*cr[1];
  r_2e=sqrt(ro2)/(2.0*veps);
  kd=2.0*M_PI/qd->d;
  k2=qd->k*qd->k;
  ve2=veps*veps;
  i_ve=1.0/veps;

  arg=kd*cr[0];
  st=sin(arg);  ct=cos(arg);
  ten=ct+I*st;
  i_ten=ct-I*st;

  n=N;
  aN=(double)n*kd-qd->kx;
  arg=aN*cr[0];
  eiaxN=cos(arg)+I*sin(arg);
  bn2=k2-aN*aN; if(bn2==0.0) return -2;
  err=sum_expintv2(&tc0,&tc1,&tc2,ro2,veps,bn2,eps);  if(err<0) return -1;
  *qG2=eiaxN*tc0;
  dqG2[0]=aN*eiaxN*tc0;
  dqG2[1]=eiaxN*tc1;
  d2qG2[0]=aN*aN*eiaxN*tc0;
  d2qG2[1]=aN*eiaxN*tc1*r_2e;
  d2qG2[2]=eiaxN*tc2;

  eiaxp=eiaxN;
  eiaxm=eiaxN;
  anp=aN;
  anm=aN;
  for(n=0;n<EW_LIMIT;n++){
    oG=*qG2;
    odG0=dqG2[0];
    odG1=dqG2[1];
    od2G0=d2qG2[0];
    od2G1=d2qG2[1];
    od2G2=d2qG2[2];

    eiaxp*=ten;
    eiaxm*=i_ten;
    anp+=kd;
    anm-=kd;
    sig=0;

    // +
    bn2=k2-anp*anp;  if(bn2==0.0) return -2;
    if(bn2<0.0) sig++;
    err=sum_expintv2(&tc0,&tc1,&tc2,ro2,veps,bn2,eps);  if(err<0) return -1;
    *qG2+=eiaxp*tc0;
    dqG2[0]+=anp*eiaxp*tc0;
    dqG2[1]+=eiaxp*tc1;
    d2qG2[0]+=anp*anp*eiaxp*tc0;
    d2qG2[1]+=anp*eiaxp*tc1*r_2e;
    d2qG2[2]+=eiaxp*tc2;

    // -
    bn2=k2-anm*anm;  if(bn2==0.0) return -2;
    if(bn2<0.0) sig++;
    err=sum_expintv2(&tc0,&tc1,&tc2,ro2,veps,bn2,eps);  if(err<0) return -1;
    *qG2+=eiaxm*tc0;
    dqG2[0]+=anm*eiaxm*tc0;
    dqG2[1]+=eiaxm*tc1;
    d2qG2[0]+=anm*anm*eiaxm*tc0;
    d2qG2[1]+=anm*eiaxm*tc1*r_2e;
    d2qG2[2]+=eiaxm*tc2;

    if(sig>0){
      if(flg[0]==0){
        if(  2.0*fabs(creal(oG)-creal(*qG2)) < eps*(fabs(creal(oG))+fabs(creal(*qG2)))
          && 2.0*fabs(cimag(oG)-cimag(*qG2)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG2))) ) flg[0]=1;
      }
      if(flg[1]==0){
        if(  2.0*fabs(creal(odG0)-creal(dqG2[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG2[0])))
          && 2.0*fabs(cimag(odG0)-cimag(dqG2[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG2[0]))) ) flg[1]=1;
      }
      if(flg[2]==0){
        if(  2.0*fabs(creal(odG1)-creal(dqG2[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG2[1])))
          && 2.0*fabs(cimag(odG1)-cimag(dqG2[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG2[1]))) ) flg[2]=1;
      }
      if(flg[3]==0){
        if(  2.0*fabs(creal(od2G0)-creal(d2qG2[0])) < eps*(fabs(creal(od2G0))+fabs(creal(d2qG2[0])))
          && 2.0*fabs(cimag(od2G0)-cimag(d2qG2[0])) < eps*(fabs(cimag(od2G0))+fabs(cimag(d2qG2[0]))) ) flg[3]=1;
      }
      if(flg[4]==0){
        if(  2.0*fabs(creal(od2G1)-creal(d2qG2[1])) < eps*(fabs(creal(od2G1))+fabs(creal(d2qG2[1])))
          && 2.0*fabs(cimag(od2G1)-cimag(d2qG2[1])) < eps*(fabs(cimag(od2G1))+fabs(cimag(d2qG2[1]))) ) flg[4]=1;
      }
      if(flg[5]==0){
        if(  2.0*fabs(creal(od2G2)-creal(d2qG2[2])) < eps*(fabs(creal(od2G2))+fabs(creal(d2qG2[2])))
          && 2.0*fabs(cimag(od2G2)-cimag(d2qG2[2])) < eps*(fabs(cimag(od2G2))+fabs(cimag(d2qG2[2]))) ) flg[5]=1;
      }
    }

    if(flg[0]!=0 && flg[1]!=0 && flg[2]!=0 && flg[3]!=0 && flg[4]!=0 && flg[5]!=0) break;
  }

  tc0=1.0/(4.0*M_PI*qd->d);
  tc1=tc0/(2.0*ve2);
  *qG2*=tc0;
  dqG2[0]*=I*tc0;
  dqG2[1]*=tc1;
  d2qG2[0]*=-tc0;
  d2qG2[1]*=I*tc0*i_ve;
  d2qG2[2]*=tc1;

  if(n!=EW_LIMIT) return 2*n+1;
  else return -1;
}

int sum_expintv2(double complex *ret0,double complex *ret1,double complex *ret2,double ro2,double veps,double bn2,double eps)
{
  double complex E1,E0,tG0,tG1,tG2,tG0o,tG1o,tG2o;
  double sig,ipm,rpm,z,im,dp,rpm1,ipm1;
  int m,flg[3];

  flg[0]=0; flg[1]=0; flg[2]=0;
  dp=ro2/(4.0*veps*veps);
  z=-bn2*veps*veps;

  sig=1.0;
  ipm=1.0;
  rpm=1.0;
  if(bn2>0.0) E1=I*M_PI-expint_i(-z);
  else if(bn2<0.0) E1=expint(1,z);
  else return -3;
  tG0=sig*ipm*rpm*E1;
  tG1=0.0;
  tG2=0.0;
  for(m=1;m<EW_LIMIT;m++){
    tG0o=tG0;
    tG1o=tG1;
    tG2o=tG2;

    sig*=-1.0;
    im=1.0/(double)m;
    ipm1=ipm;
    ipm*=im;
    rpm1=rpm;
    rpm*=dp;
    E0=E1;
    E1=im*(exp(-z)-z*E0);
    tG0+=sig*ipm*rpm*E1;
    tG1+=sig*ipm1*rpm1*E1;
    tG2+=sig*ipm1*rpm1*E1*(double)(2*m-1);

    if(bn2>0.0){
      if(flg[0]==0){
        if(  2.0*fabs(creal(tG0o)-creal(tG0)) < eps*(fabs(creal(tG0o))+fabs(creal(tG0)))
          && 2.0*fabs(cimag(tG0o)-cimag(tG0)) < eps*(fabs(cimag(tG0o))+fabs(cimag(tG0))) ) flg[0]=1;
      }
      if(flg[1]==0){
        if(  2.0*fabs(creal(tG1o)-creal(tG1)) < eps*(fabs(creal(tG1o))+fabs(creal(tG1)))
          && 2.0*fabs(cimag(tG1o)-cimag(tG1)) < eps*(fabs(cimag(tG1o))+fabs(cimag(tG1))) ) flg[1]=1;
      }
      if(flg[2]==0){
        if(  2.0*fabs(creal(tG2o)-creal(tG2)) < eps*(fabs(creal(tG2o))+fabs(creal(tG2)))
          && 2.0*fabs(cimag(tG2o)-cimag(tG2)) < eps*(fabs(cimag(tG2o))+fabs(cimag(tG2))) ) flg[2]=1;
      }
    }
    else{
      if(flg[0]==0){
        if(  2.0*fabs(creal(tG0o)-creal(tG0)) < eps*(fabs(creal(tG0o))+fabs(creal(tG0))) ) flg[0]=1;
      }
      if(flg[1]==0){
        if(  2.0*fabs(creal(tG1o)-creal(tG1)) < eps*(fabs(creal(tG1o))+fabs(creal(tG1))) ) flg[1]=1;
      }
      if(flg[2]==0){
        if(  2.0*fabs(creal(tG2o)-creal(tG2)) < eps*(fabs(creal(tG2o))+fabs(creal(tG2))) ) flg[2]=1;
      }
    }
    if(flg[0]!=0 && flg[1]!=0 && flg[2]!=0) break;
  }

  *ret0=tG0;
  *ret1=tG1;
  *ret2=tG2;

  if(m==EW_LIMIT) return -1;
  else return 2*m-1;
}

int d3hm_qpgf_d1_cv1_fd_qG(double complex *qG,double *cr,double eps,QPDT *qd)
{ // Fourier domain summation
  int d3hm_qpgf_d1_cv1_fd_qG_sd(double complex *qG,double *cr,double eps,QPDT *qd);

  double complex c;
  double rt[3],tmp;
  int err,l;

  rt[0]=cr[0];
  rt[1]=cr[1];
  rt[2]=cr[2];

  l=(int)floor(cr[0]/qd->d+0.5);
  rt[0]-=(double)l*qd->d;

  err=d3hm_qpgf_d1_cv1_fd_qG_sd(qG,rt,eps,qd);
  if(l!=0){
    tmp=(double)l*qd->kx*qd->d;
    c=cos(tmp)-I*sin(tmp);
    *qG*=c;
  }

  return err;
}

int d3hm_qpgf_d1_cv1_fd_qG_sd(double complex *qG,double *cr,double eps,QPDT *qd)
{
  double complex ekx,ten,i_ten,eiaxN,eiaxp,eiaxm,ch,tc,tf,oG;
  double arg,ro,k2,kd,kdx,bn2,aN,anp,anm,bn,bnr,ct,st;
  int N,n,sig;

  ro=cr[1];
  k2=qd->k*qd->k;

  N=(int)round(qd->kx*qd->d/(2.0*M_PI));
  ch=M_PI*I/2.0;

  arg=-qd->kx*cr[0];
  ekx=cos(arg)+I*sin(arg);

  kd=2.0*M_PI/qd->d;
  kdx=kd*cr[0];
  ct=cos(kdx);
  st=sin(kdx);
  ten=ct+I*st;
  i_ten=ct-I*st;

  n=N;
  aN=kd*(double)n-qd->kx;
  arg=kdx*(double)n;
  eiaxN=(cos(arg)+I*sin(arg))*ekx;
  bn2=k2-aN*aN;
  if(bn2>0.0){
    bn=sqrt(bn2);
    bnr=bn*ro;
    tc=eiaxN*ch;
    tf=tc*(besj0(bnr)+I*besy0(bnr));
    *qG=tf;
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2);
    bnr=bn*ro;
    tf=eiaxN*bessk0(bnr);
    *qG=tf;
  }
  else return -2; // beta_n=0.0;

  eiaxp=eiaxN;
  eiaxm=eiaxN;
  anp=aN;
  anm=aN;
  for(n=0;n<FD_LIMIT;n++){
    oG=*qG;
    sig=0;

    eiaxp*=ten;
    eiaxm*=i_ten;
    anp+=kd;
    anm-=kd;

    // +
    bn2=k2-anp*anp;
    if(bn2>0.0){
      bn=sqrt(bn2);
      bnr=bn*ro;
      tc=eiaxp*ch;
      tf=tc*(besj0(bnr)+I*besy0(bnr));
      *qG+=tf;
      sig+=1;
    }
    else if(bn2<0.0){
      bn=sqrt(-bn2);
      bnr=bn*ro;
      tf=eiaxp*bessk0(bnr);
      *qG+=tf;
    }
    else return -2; // beta_n=0.0;

    // -
    bn2=k2-anm*anm;
    if(bn2>0.0){
      bn=sqrt(bn2);
      bnr=bn*ro;
      tc=eiaxm*ch;
      tf=tc*(besj0(bnr)+I*besy0(bnr));
      *qG+=tf;
      sig+=1;
    }
    else if(bn2<0.0){
      bn=sqrt(-bn2);
      bnr=bn*ro;
      tf=eiaxm*bessk0(bnr);
      *qG+=tf;
    }
    else return -2; // beta_n=0.0;

    if(sig==0){
      if(  2.0*fabs(creal(oG)-creal(*qG)) < eps*(fabs(creal(oG))+fabs(creal(*qG)))
        && 2.0*fabs(cimag(oG)-cimag(*qG)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG))) ) break;
    }
  }

  if(n==FD_LIMIT) return -1;

  tc=1.0/(2.0*M_PI*qd->d);
  *qG*=tc;

  return 2*n+1;
}

int d3hm_qpgf_d1_cv1_ew_qG(double complex *qG,double *cr,double eps,QPDT *qd)
{ // Ewald method
  int d3hm_ewald_qG1_cv1_qG(double complex *qG1,double *cr,double eps,double veps,QPDT *qd);
  int d3hm_ewald_qG2_cv1_qG(double complex *qG1,double *cr,double eps,double veps,QPDT *qd);

  double complex tqG;
  double veps,r02;
  int err1,err2;

  veps=qd->d/(2.0*sqrt(M_PI));

  // veps
  r02=cr[0]*cr[0]+cr[1]*cr[1];
  if( (-r02/(4.0*veps*veps)+qd->k*qd->k*veps*veps) > EW_EXPM)
    veps=sqrt(EW_EXPM+sqrt(qd->k*qd->k*r02+EW_EXPM*EW_EXPM))/(sqrt(2.0)*qd->k);

  // qG1
  err1=d3hm_ewald_qG1_cv1_qG(&tqG,cr,eps,veps,qd);
  *qG=tqG;
  if(err1==-1) return -1;

  // qG2
  err2=d3hm_ewald_qG2_cv1_qG(&tqG,cr,eps,veps,qd);
  *qG+=tqG;
  if(err2==-1) return -2;
  if(err2==-2) return -3; // beta_n==0.0

  return err1+err2;
}

int d3hm_ewald_qG1_cv1_qG(double complex *qG1,double *cr,double eps,double veps,QPDT *qd)
{
  double complex ekd,i_ekd,eLkd,tw,elkdp,elkdm,tce,oG;
  double kd,ro2,rl2,rl,i_rl,xLd,arg,ke,ke2,i2e,mi2e2,xldp,xldm,st,ct,td;
  int L,l;

  kd=qd->kx*qd->d;
  st=sin(kd);
  ct=cos(kd);
  ekd=ct+I*st;
  i_ekd=ct-I*st;

  ro2=cr[1]*cr[1];
  ke=qd->k*veps;
  ke2=ke*ke;
  i2e=1.0/(2.0*veps);
  mi2e2=-i2e*i2e;

  if(fabs(cr[0])<qd->d*sqrt(eps)) L=1; // avoid symmetric cancellation
  else L=0;

  arg=kd*(double)L;
  eLkd=cos(arg)+I*sin(arg);
  xLd=cr[0]+(double)L*qd->d;
  rl2=xLd*xLd+ro2;
  rl=sqrt(rl2);
  i_rl=1.0/rl;
  tw=Faddeeva_w(ke+I*rl*i2e,eps);
  tce=eLkd*i_rl*exp(rl2*mi2e2+ke2);
  *qG1=tce*creal(tw);

  elkdp=eLkd;
  elkdm=eLkd;
  xldp=xLd;
  xldm=xLd;
  for(l=0;l<EW_LIMIT;l++){
    oG=*qG1;
    elkdp*=ekd;
    elkdm*=i_ekd;
    xldp+=qd->d;
    xldm-=qd->d;

    // +
    rl2=xldp*xldp+ro2;
    rl=sqrt(rl2);
    i_rl=1.0/rl;
    tw=Faddeeva_w(ke+I*rl*i2e,eps);
    tce=elkdp*i_rl*exp(rl2*mi2e2+ke2);
    *qG1+=tce*creal(tw);

    //-
    rl2=xldm*xldm+ro2;
    rl=sqrt(rl2);
    i_rl=1.0/rl;
    tw=Faddeeva_w(ke+I*rl*i2e,eps);
    tce=elkdm*i_rl*exp(rl2*mi2e2+ke2);
    *qG1+=tce*creal(tw);

    if(  2.0*fabs(creal(oG)-creal(*qG1)) < eps*(fabs(creal(oG))+fabs(creal(*qG1)))
      && 2.0*fabs(cimag(oG)-cimag(*qG1)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG1))) ) break;
  }

  td=1.0/(4.0*M_PI);
  *qG1*=td;

  if(l!=EW_LIMIT) return 2*l+1;
  else return -1;
}

int d3hm_ewald_qG2_cv1_qG(double complex *qG2,double *cr,double eps,double veps,QPDT *qd)
{
  int sum_expintv_qG(double complex *ret0,double ro2,double veps,double bn2,double eps);

  double complex ten,i_ten,eiaxN,tc0,eiaxp,eiaxm,oG;
  double ro2,kd,k2,arg,ct,st,aN,bn2,anp,anm;
  int N,n,err,sig;

  N=(int)round(qd->kx*qd->d/(2.0*M_PI));

  ro2=cr[1]*cr[1];
  kd=2.0*M_PI/qd->d;
  k2=qd->k*qd->k;

  arg=kd*cr[0];
  st=sin(arg);  ct=cos(arg);
  ten=ct+I*st;
  i_ten=ct-I*st;

  n=N;
  aN=(double)n*kd-qd->kx;
  arg=aN*cr[0];
  eiaxN=cos(arg)+I*sin(arg);
  bn2=k2-aN*aN; if(bn2==0.0) return -2;
  err=sum_expintv_qG(&tc0,ro2,veps,bn2,eps);  if(err<0) return -1;
  *qG2=eiaxN*tc0;

  eiaxp=eiaxN;
  eiaxm=eiaxN;
  anp=aN;
  anm=aN;
  for(n=0;n<EW_LIMIT;n++){
    oG=*qG2;
    eiaxp*=ten;
    eiaxm*=i_ten;
    anp+=kd;
    anm-=kd;
    sig=0;

    // +
    bn2=k2-anp*anp;  if(bn2==0.0) return -2;
    if(bn2<0.0) sig++;
    err=sum_expintv_qG(&tc0,ro2,veps,bn2,eps);  if(err<0) return -1;
    *qG2+=eiaxp*tc0;

    // -
    bn2=k2-anm*anm;  if(bn2==0.0) return -2;
    if(bn2<0.0) sig++;
    err=sum_expintv_qG(&tc0,ro2,veps,bn2,eps);  if(err<0) return -1;
    *qG2+=eiaxm*tc0;

    if(sig>0){
      if(  2.0*fabs(creal(oG)-creal(*qG2)) < eps*(fabs(creal(oG))+fabs(creal(*qG2)))
        && 2.0*fabs(cimag(oG)-cimag(*qG2)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG2))) ) break;
    }
  }

  tc0=1.0/(4.0*M_PI*qd->d);
  *qG2*=tc0;

  if(n!=EW_LIMIT) return 2*n+1;
  else return -1;
}

int sum_expintv_qG(double complex *ret0,double ro2,double veps,double bn2,double eps)
{
  double complex E1,E0,tG0,tG0o;
  double sig,ipm,rpm,z,im,dp;
  int m;

  dp=ro2/(4.0*veps*veps);
  z=-bn2*veps*veps;

  sig=1.0;
  ipm=1.0;
  rpm=1.0;
  if(bn2>0.0) E1=I*M_PI-expint_i(-z);
  else if(bn2<0.0) E1=expint(1,z);
  else return -3;
  tG0=sig*ipm*rpm*E1;
  for(m=1;m<EW_LIMIT;m++){
    tG0o=tG0;
    sig*=-1.0;
    im=1.0/(double)m;
    ipm*=im;
    rpm*=dp;
    E0=E1;
    E1=im*(exp(-z)-z*E0);
    tG0+=sig*ipm*rpm*E1;

    if(bn2>0.0){
      if(  2.0*fabs(creal(tG0o)-creal(tG0)) < eps*(fabs(creal(tG0o))+fabs(creal(tG0)))
        && 2.0*fabs(cimag(tG0o)-cimag(tG0)) < eps*(fabs(cimag(tG0o))+fabs(cimag(tG0))) ) break;
    }
    else {
      if(  2.0*fabs(creal(tG0o)-creal(tG0)) < eps*(fabs(creal(tG0o))+fabs(creal(tG0))) ) break;
    }
  }

  *ret0=tG0;

  if(m==EW_LIMIT) return -1;
  else return 2*m-1;
}

int d3hm_qpgf_d1_cv1_fd_dqG(double complex *dqG,double *cr,double eps,QPDT *qd)
{ // Fourier domain summation
  int d3hm_qpgf_d1_cv1_fd_sd_dqG(double complex *dqG,double *cr,double eps,QPDT *qd);

  double complex c;
  double rt[3],tmp;
  int err,l;

  rt[0]=cr[0];
  rt[1]=cr[1];
  rt[2]=cr[2];

  l=(int)floor(cr[0]/qd->d+0.5);
  rt[0]-=(double)l*qd->d;

  err=d3hm_qpgf_d1_cv1_fd_sd_dqG(dqG,rt,eps,qd);
  if(l!=0){
    tmp=(double)l*qd->kx*qd->d;
    c=cos(tmp)-I*sin(tmp);
    dqG[0]*=c;
    dqG[1]*=c;
  }

  return err;
}

int d3hm_qpgf_d1_cv1_fd_sd_dqG(double complex *dqG,double *cr,double eps,QPDT *qd)
{
  double complex ekx,ten,i_ten,eiaxN,eiaxp,eiaxm,ch,tc,tf,odG0,odG1;
  double arg,ro,k2,kd,kdx,bn2,aN,anp,anm,bn,bnr,ct,st;
  int i,N,flg[2],n,sig;

  ro=cr[1];
  k2=qd->k*qd->k;

  N=(int)round(qd->kx*qd->d/(2.0*M_PI));
  ch=M_PI*I/2.0;

  arg=-qd->kx*cr[0];
  ekx=cos(arg)+I*sin(arg);

  kd=2.0*M_PI/qd->d;
  kdx=kd*cr[0];
  ct=cos(kdx);
  st=sin(kdx);
  ten=ct+I*st;
  i_ten=ct-I*st;

  for(i=0;i<2;i++) flg[i]=0.0;

  n=N;
  aN=kd*(double)n-qd->kx;
  arg=kdx*(double)n;
  eiaxN=(cos(arg)+I*sin(arg))*ekx;
  bn2=k2-aN*aN;
  if(bn2>0.0){
    bn=sqrt(bn2);
    bnr=bn*ro;
    tc=eiaxN*ch;
    tf=tc*(besj0(bnr)+I*besy0(bnr));
    dqG[0]=aN*tf;
    dqG[1]=-bn*tc*(besj1(bnr)+I*besy1(bnr));
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2);
    bnr=bn*ro;
    tf=eiaxN*bessk0(bnr);
    dqG[0]=aN*tf;
    dqG[1]=-bn*eiaxN*bessk1(bnr);
  }
  else return -2; // beta_n=0.0;

  eiaxp=eiaxN;
  eiaxm=eiaxN;
  anp=aN;
  anm=aN;
  for(n=0;n<FD_LIMIT;n++){
    odG0=dqG[0];
    odG1=dqG[1];
    sig=0;

    eiaxp*=ten;
    eiaxm*=i_ten;
    anp+=kd;
    anm-=kd;

    // +
    bn2=k2-anp*anp;
    if(bn2>0.0){
      bn=sqrt(bn2);
      bnr=bn*ro;
      tc=eiaxp*ch;
      tf=tc*(besj0(bnr)+I*besy0(bnr));
      dqG[0]+=anp*tf;
      dqG[1]+=-bn*tc*(besj1(bnr)+I*besy1(bnr));
      sig+=1;
    }
    else if(bn2<0.0){
      bn=sqrt(-bn2);
      bnr=bn*ro;
      tf=eiaxp*bessk0(bnr);
      dqG[0]+=anp*tf;
      dqG[1]+=-bn*eiaxp*bessk1(bnr);
    }
    else return -2; // beta_n=0.0;

    // -
    bn2=k2-anm*anm;
    if(bn2>0.0){
      bn=sqrt(bn2);
      bnr=bn*ro;
      tc=eiaxm*ch;
      tf=tc*(besj0(bnr)+I*besy0(bnr));
      dqG[0]+=anm*tf;
      dqG[1]+=-bn*tc*(besj1(bnr)+I*besy1(bnr));
      sig+=1;
    }
    else if(bn2<0.0){
      bn=sqrt(-bn2);
      bnr=bn*ro;
      tf=eiaxm*bessk0(bnr);
      dqG[0]+=anm*tf;
      dqG[1]+=-bn*eiaxm*bessk1(bnr);
    }
    else return -2; // beta_n=0.0;

    if(sig==0){
      if(flg[0]==0){
        if(  2.0*fabs(creal(odG0)-creal(dqG[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG[0])))
          && 2.0*fabs(cimag(odG0)-cimag(dqG[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG[0]))) ) flg[0]=1;
      }
      if(flg[1]==0){
        if(  2.0*fabs(creal(odG1)-creal(dqG[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG[1])))
          && 2.0*fabs(cimag(odG1)-cimag(dqG[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG[1]))) ) flg[1]=1;
      }
    }

    if(flg[0]!=0 && flg[1]!=0) break;
  }

  if(n==FD_LIMIT) return -1;

  tc=1.0/(2.0*M_PI*qd->d);
  dqG[0]*=I*tc;

  tc/=ro;
  dqG[1]*=tc;
  dqG[2]=0.0;

  return 2*n+1;
}

int d3hm_qpgf_d1_cv1_ew_dqG(double complex *dqG,double *cr,double eps,QPDT *qd)
{ // Ewald method
  int d3hm_ewald_qG1_cv1_dqG(double complex *dqG1,double *cr,double eps,double veps,QPDT *qd);
  int d3hm_ewald_qG2_cv1_dqG(double complex *dqG1,double *cr,double eps,double veps,QPDT *qd);

  double complex tdqG[3];
  double veps,r02;
  int err1,err2;

  veps=qd->d/(2.0*sqrt(M_PI));

  // veps
  r02=cr[0]*cr[0]+cr[1]*cr[1];
  if( (-r02/(4.0*veps*veps)+qd->k*qd->k*veps*veps) > EW_EXPM)
    veps=sqrt(EW_EXPM+sqrt(qd->k*qd->k*r02+EW_EXPM*EW_EXPM))/(sqrt(2.0)*qd->k);

  // qG1
  err1=d3hm_ewald_qG1_cv1_dqG(tdqG,cr,eps,veps,qd);
  dqG[0]=tdqG[0];
  dqG[1]=tdqG[1];
  if(err1==-1) return -1;

  // qG2
  err2=d3hm_ewald_qG2_cv1_dqG(tdqG,cr,eps,veps,qd);
  dqG[0]+=tdqG[0];
  dqG[1]+=tdqG[1];
  if(err2==-1) return -2;
  if(err2==-2) return -3; // beta_n==0.0

  dqG[2]=0.0;
  return err1+err2;
}

int d3hm_ewald_qG1_cv1_dqG(double complex *dqG1,double *cr,double eps,double veps,QPDT *qd)
{
  double complex ekd,i_ekd,eLkd,tw,elkdp,elkdm,tce,tcf,odG0,odG1;
  double kd,ro2,rl2,rl,i_rl,xLd,arg,ke,ke2,i2e,mi2e2,xldp,xldm,st,ct,td,cpe;
  int L,l,flg[2];

  kd=qd->kx*qd->d;
  st=sin(kd);
  ct=cos(kd);
  ekd=ct+I*st;
  i_ekd=ct-I*st;

  ro2=cr[1]*cr[1];
  ke=qd->k*veps;
  ke2=ke*ke;
  i2e=1.0/(2.0*veps);
  mi2e2=-i2e*i2e;
  cpe=1.0/(sqrt(M_PI)*veps);

  if(fabs(cr[0])<qd->d*sqrt(eps)) L=1; // avoid symmetric cancellation
  else L=0;
  flg[0]=0;  flg[1]=0;

  arg=kd*(double)L;
  eLkd=cos(arg)+I*sin(arg);
  xLd=cr[0]+(double)L*qd->d;
  rl2=xLd*xLd+ro2;
  rl=sqrt(rl2);
  i_rl=1.0/rl;
  tw=Faddeeva_w(ke+I*rl*i2e,eps);
  tce=eLkd*i_rl*exp(rl2*mi2e2+ke2);
  tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
  dqG1[0]=xLd*tcf;
  dqG1[1]=tcf;

  elkdp=eLkd;
  elkdm=eLkd;
  xldp=xLd;
  xldm=xLd;
  for(l=0;l<EW_LIMIT;l++){
    odG0=dqG1[0];
    odG1=dqG1[1];
    elkdp*=ekd;
    elkdm*=i_ekd;
    xldp+=qd->d;
    xldm-=qd->d;

    // +
    rl2=xldp*xldp+ro2;
    rl=sqrt(rl2);
    i_rl=1.0/rl;
    tw=Faddeeva_w(ke+I*rl*i2e,eps);
    tce=elkdp*i_rl*exp(rl2*mi2e2+ke2);
    tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
    dqG1[0]+=xldp*tcf;
    dqG1[1]+=tcf;

    //-
    rl2=xldm*xldm+ro2;
    rl=sqrt(rl2);
    i_rl=1.0/rl;
    tw=Faddeeva_w(ke+I*rl*i2e,eps);
    tce=elkdm*i_rl*exp(rl2*mi2e2+ke2);
    tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
    dqG1[0]+=xldm*tcf;
    dqG1[1]+=tcf;

    if(flg[0]==0){
      if(  2.0*fabs(creal(odG0)-creal(dqG1[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG1[0])))
        && 2.0*fabs(cimag(odG0)-cimag(dqG1[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG1[0]))) ) flg[0]=1;
    }
    if(flg[1]==0){
      if(  2.0*fabs(creal(odG1)-creal(dqG1[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG1[1])))
        && 2.0*fabs(cimag(odG1)-cimag(dqG1[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG1[1]))) ) flg[1]=1;
    }

    if(flg[0]!=0 && flg[1]!=0) break;
  }

  td=1.0/(4.0*M_PI);
  dqG1[0]*=-td;
  dqG1[1]*=-td;

  if(l!=EW_LIMIT) return 2*l+1;
  else return -1;
}

int d3hm_ewald_qG2_cv1_dqG(double complex *dqG2,double *cr,double eps,double veps,QPDT *qd)
{
  int sum_expintv(double complex *ret0,double complex *ret1,double ro2,double veps,double bn2,double eps);

  double complex ten,i_ten,eiaxN,tc0,tc1,eiaxp,eiaxm,odG0,odG1;
  double ro2,kd,k2,ve2,arg,ct,st,aN,bn2,anp,anm;
  int N,n,err,flg[2],sig;

  flg[0]=0;  flg[1]=0;
  N=(int)round(qd->kx*qd->d/(2.0*M_PI));

  ro2=cr[1]*cr[1];
  kd=2.0*M_PI/qd->d;
  k2=qd->k*qd->k;
  ve2=veps*veps;

  arg=kd*cr[0];
  st=sin(arg);  ct=cos(arg);
  ten=ct+I*st;
  i_ten=ct-I*st;

  n=N;
  aN=(double)n*kd-qd->kx;
  arg=aN*cr[0];
  eiaxN=cos(arg)+I*sin(arg);
  bn2=k2-aN*aN; if(bn2==0.0) return -2;
  err=sum_expintv(&tc0,&tc1,ro2,veps,bn2,eps);  if(err<0) return -1;
  dqG2[0]=aN*eiaxN*tc0;
  dqG2[1]=eiaxN*tc1;

  eiaxp=eiaxN;
  eiaxm=eiaxN;
  anp=aN;
  anm=aN;
  for(n=0;n<EW_LIMIT;n++){
    odG0=dqG2[0];
    odG1=dqG2[1];
    eiaxp*=ten;
    eiaxm*=i_ten;
    anp+=kd;
    anm-=kd;
    sig=0;

    // +
    bn2=k2-anp*anp;  if(bn2==0.0) return -2;
    if(bn2<0.0) sig++;
    err=sum_expintv(&tc0,&tc1,ro2,veps,bn2,eps);  if(err<0) return -1;
    dqG2[0]+=anp*eiaxp*tc0;
    dqG2[1]+=eiaxp*tc1;

    // -
    bn2=k2-anm*anm;  if(bn2==0.0) return -2;
    if(bn2<0.0) sig++;
    err=sum_expintv(&tc0,&tc1,ro2,veps,bn2,eps);  if(err<0) return -1;
    dqG2[0]+=anm*eiaxm*tc0;
    dqG2[1]+=eiaxm*tc1;

    if(sig>0){
      if(flg[0]==0){
        if(  2.0*fabs(creal(odG0)-creal(dqG2[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG2[0])))
          && 2.0*fabs(cimag(odG0)-cimag(dqG2[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG2[0]))) ) flg[0]=1;
      }
      if(flg[1]==0){
        if(  2.0*fabs(creal(odG1)-creal(dqG2[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG2[1])))
          && 2.0*fabs(cimag(odG1)-cimag(dqG2[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG2[1]))) ) flg[1]=1;
      }
    }

    if(flg[0]!=0 && flg[1]!=0) break;
  }

  tc0=1.0/(4.0*M_PI*qd->d);
  tc1=tc0/(2.0*ve2);
  dqG2[0]*=I*tc0;
  dqG2[1]*=tc1;

  if(n!=EW_LIMIT) return 2*n+1;
  else return -1;
}
