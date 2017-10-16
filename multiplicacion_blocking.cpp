#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "papi.h"

const int csize = 32;

int code_to_be_measured(const double * M1, const double * M2, double * M);

int main(int argc, char **argv)
{
  int N = atoi(argv[1]); //We run N with a for loop in bash
  
  // Matrix declaration : Modeled as 1D array
  // Declare as pointers and ask for memory to use the heap
  double *A = new double [N*N], *B = new double [N*N], *C = new double [N*N];
  // initialize matrices
  for (int ii =0; ii < N; ++ii) {
    for (int jj =0; jj < N; ++jj) {
      A[ii*N + jj] = ii + jj + 1;
      B[ii*N + jj] = ii + jj +2; //We could use any other deffinition tho
      C[ii*N + jj] = 0.0;
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
  code_to_be_measured(A, B, C);
  if((retval=PAPI_flops( &real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
    {
      printf("retval: %d\n", retval);
      exit(1);
    }
  //printf("Real_time: %f Proc_time: %f Total flpops: %lld MFLOPS: %f\n", real_time, proc_time,flpops,mflops);
  std::cout<<N<<" "<<proc_time <<" "<<mflops<<std::endl; //Printing matrix size, cpu time and mfolps.
  //In theory the Mflops should be normalized as we're not using divison nor any complicated operations.
  
  delete [] A;
  delete [] AT;
  return 0;
}

int code_to_be_measured(const double * M1, const double * M2, double * M)
{
  // matrix multiplication with blocking
  for (int ii = 0; ii < N; ii+=csize)
    for (int jj = 0; jj < N; jj+=csize)
      for (int kk = 0; kk < N; kk+=csize)
	for(int i=ii; i<min(N, ii+csize-1);++i)
	  for(int j = jj; j < min(N;jj+csi<e-1);++j)
	    for(int j = jj; j < min(N;jj+csi<e-1);++j)
	      M[i*N +j] += M1[i*N + k]*M2[k*N + j];
  return 0;
}
