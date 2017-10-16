#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "papi.h"

int code_to_be_measured(const double * M, double * MT);

int main(int argc, char **argv)
{
  int N = atoi(argv[1]); //Run N in a loop in bash
  //As int for M in 2 4 8 16 32..16384; do ./a.out $M; done
  // Matrix declaration : Modeled as 1D array
  // Declare as pointers and ask for memory to use the heap
  double *A = new double [N*N], *AT = new double [N*N];
  // initialize matrices
  for (int ii =0; ii < N; ++ii) {
    for (int jj =0; jj < N; ++jj) {
      A[ii*N + jj] = ii + jj + 1;
      AT[ii*N + jj] = 0.0;
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
  std::cout<<N<<" "<<proc_time<<" "<<mflops<<std::endl;
  delete [] A;
  delete [] AT;
  return 0;
}

int code_to_be_measured(const double * M, double * MT)
{
  // simple matrix transpose
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      MT[ii*N +jj] = M[jj*N + ii];
    }
  }
  return 0;
}
