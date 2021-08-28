#include "bem3_emf_qd1.h"

int main(int argc,char **argv)
{
  DQD1 qd1;
  
  read_dqd1(argc,argv,&qd1);
  print_dqd1(&qd1);
  //print_dqd1_mksa(&qd1);
  initialize_dqd1(&qd1);

  solve_bieq(&qd1);
  printf("\noutput data file : %s\n",argv[5]);
  dat_write(argv[5],&qd1);
  finalize_dqd1(&qd1);
  
  return 0;
}
