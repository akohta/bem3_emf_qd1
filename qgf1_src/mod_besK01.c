#include "mod_besK01.h"

double bessk0(double x)
{
  double polyn(double *f,int nsize,int nstart,double x);
  
  static double k0pi[5]={1.425433617130587e-6, 2.150707366040937e-4, 1.187082088663404e-2, 2.346487949187396e-1,1.0};
  static double k0qi[3]={8.362215678646257e-5, 1.518396076767770e-2, 9.847324170755358e-1};
  static double k0p [5]={3.454715527986737e-6, 4.574734709978264e-4, 2.066458134619875e-2, 2.770731240515333e-1, 1.159315156584126e-1};
  static double k0q [3]={9.809660603621949e-5, 1.627693622304549e-2, 9.836249671709183e-1};
  static double k0pp[8]={6.160228690102976e-2, 3.595376024148513e+0, 3.198289277679660e+1, 9.285288485892228e+1,
      1.121012633939949e+2, 6.123767403223466e+1, 1.475731032429900e+1, 1.253314137315499e+0};
  static double k0qq[8]={1.408295601966600e-1, 4.443672926432041e+0, 3.181399777449301e+1, 8.318077493230258e+1,
      9.496513373427093e+1, 5.027773590829784e+1, 1.189963006673403e+1, 1.0};

  double tx,tm;
  if (x <= 1.0) {
    tx=x*x;
    tm = polyn(k0pi,5,0,tx)*log(x)/polyn(k0qi,3,0,1.0-tx);
    return polyn(k0p,5,0,tx)/polyn(k0q,3,0,1.0-tx)-tm;
  } else {
    tx=1.0/x;
    return exp(-x)*polyn(k0pp,8,0,tx)/(polyn(k0qq,8,0,tx)*sqrt(x));
  }
}

double bessk1(double x)
{
  double polyn(double *f,int nsize,int nstart,double x);
  
  static double k1pi[5]={1.239567816344855e-7, 2.397509908859959e-5, 1.818666382168295e-3, 5.598072040178741e-2, 0.5};
  static double k1qi[3]={5.881933053917096e-5, 1.292092053534579e-2, 9.870202601341150e-1};
  static double k1p [5]={-3.110372465429008e-7, -5.385594871975406e-5, -3.477550948593604e-3, -8.109417631822442e-2, -3.079657578292062e-1};
  static double k1q [3]={6.774221332947002e-5, 1.375094061153160e-2, 9.861813171751389e-1};
  static double k1pp[8]={3.144418558991021e-1, 7.265230396353690, 4.356596656837691e1, 1.040442011439181e2,
      1.147386690867892e2, 6.063161173098803e1, 1.457171340220454e1, 1.253314137315502};
  static double k1qq[8]={2.538540887654872e-2, 1.857244676566022, 1.850303673841586e1, 5.863377227890893e1,
      7.616113213117645e1, 4.427488496597630e1, 1.125154514806458e1, 1.0};

  double tx,tm;
  if (x <= 1.0) {
    tx=x*x;
    tm = polyn(k1pi,5,0,tx)*log(x)/polyn(k1qi,3,0,1.0-tx);
    return x*(polyn(k1p,5,0,tx)/polyn(k1q,3,0,1.0-tx)+tm)+1.0/x;
  } else {
    tx=1.0/x;
    return exp(-x)*polyn(k1pp,8,0,tx)/(polyn(k1qq,8,0,tx)*sqrt(x));
  }
}

double polyn(double *f,int nsize,int nstart,double x)
{
  double ret;
  int n;

  ret=f[nstart];
  for(n=nstart+1;n<nsize;n++){
    ret=ret*x+f[n];
  }
  return ret;
}

