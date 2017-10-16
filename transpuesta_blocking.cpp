#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "papi.h"

const int csize = 32; //The optimal csize from the previous study

int code_to_be_measured(const double * M, double * MT);

int main(int argc, char **argv)
{
  int N = atoi(argv[1]); //We run N with a for loop in bash
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
  // matrix transpose with blocking
  for (int ii = 0; ii < N; ii+=csize) {
    for (int jj = 0; jj < N; jj+=csize) {
      for(int i=ii; i<min(N, ii+csize-1);++i){
	for(int j = jj; j < min(N;jj+csi<e-1);++j){
	  MT[i*N +j] = M[j*N + i];
	}
      }
    }
  }
  return 0;
}
