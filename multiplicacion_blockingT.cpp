#include <iostream>
#include <cstdlib>
#include<vector>

const int csize = 20;

std::vector<double> code_to_be_measured(const std::vector<double> M1, const std::vector<double> M2, int N);

int main(int argc, char **argv)
{
  int N = atoi(argv[1]); //We run N with a for loop in bash
  
  // Matrix declaration : Modeled as 1D array
  std::vector<double> A(N*N), B(N*N), C(N*N);
  // initialize matrices
  for (int ii =0; ii < N; ++ii) {
    for (int jj =0; jj < N; ++jj) {
      A[ii*N + jj] = ii + jj + 1;
      B[ii*N + jj] = ii + jj +2; //We could use any other deffinition tho
      C[ii*N + jj] = 0.0;
    }
  }
  
  C = code_to_be_measured(A, B, N);
  
  std::cout<<"Begin test"<<std::endl;
  std::cout<<A[1*N+0]<<A[1*N+1]<<A[1*N+2]<<A[1*N+3]<<std::endl;
  std::cout<<B[0*N+2]<<B[1*N+2]<<B[2*N+2]<<B[3*N+2]<<std::endl;
  std::cout<<C[1*N+2]<<std::endl;
  std::cout<<"End test"<<std::endl;  

  
  return 0;
}

std::vector<double> code_to_be_measured(const std::vector<double> M1, const std::vector<double> M2, int N)
{
  std::vector<double> M(N*N);
  //simple matrix multiplication
  for (int ii = 0; ii < N; ii+=csize)
    for (int jj = 0; jj < N; jj+=csize)
      for (int kk = 0; kk < N; kk+=csize)
	for(int i=ii; i < std::min(N, ii+csize);++i)
	  for(int j = jj; j < std::min(N, jj+csize);++j)
	    for(int k = kk; k < std::min(N, kk+csize);++k)
	      M[i*N +j] += M1[i*N + k]*M2[k*N + j];
  return M;
}
