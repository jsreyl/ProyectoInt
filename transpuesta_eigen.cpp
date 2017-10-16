#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include "papi.h"
const int N = 1024;

int code_to_be_measured(const Eigen::MatrixXd & M, Eigen::MatrixXd & MT);

int main(int argc, char **argv)
{
  int N = atoi(argv[1]);
  
  // Matrix declaration
  Eigen::MatrixXd A(N, N), AT(N, N);
  // initialize matrices
  for (int ii =0; ii < N; ++ii) {
    for (int jj =0; jj < N; ++jj) {
      A(ii, jj) = ii + jj + 1;
      AT(ii, jj) = 0.0;
    }
  }
  
  // PAPI vars
  float real_time, proc_time,mflops;
  long long flpops;
  float ireal_time, iproc_time, imflops;
  long long iflpops;
  int retval;
  
  // PERFOMANCE MEASURE
  // start PAPI counters
  if((retval=PAPI_flops(&ireal_time,&iproc_time,&iflpops,&imflops)) < PAPI_OK)
    {
      printf("Could not initialise PAPI_flops \n");
      printf("Your platform may not support floating point operation event.\n");
      printf("retval: %d\n", retval);
      exit(1);
    }
  code_to_be_measured(A, AT);
  if((retval=PAPI_flops( &real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
    {
      printf("retval: %d\n", retval);
      exit(1);
    }
  //  printf("Real_time: %f Proc_time: %f Total flpops: %lld MFLOPS: %f\n", real_time, proc_time,flpops,mflops);
  std::cout<< N<<" "<<proc_time<<" "<<mflops<<std::endl;
  return 0;
}

int code_to_be_measured(const Eigen::MatrixXd & M, Eigen::MatrixXd & MT)
{
  MT = M.transpose().eval();
  return 0;
}
