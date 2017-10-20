#include <iostream>
#include <cstdlib>
#include<vector>

const int csize = 32; //The optimal csize from the previous study

std::vector<double> code_to_be_measured(const std::vector<double> M, int N);

int main(int argc, char **argv)
{
  int N = atoi(argv[1]); //We run N with a for loop in bash
  // Matrix declaration : Modeled as 1D array

  std::vector<double> A(N*N), AT(N*N);
  // initialize matrices
  for (int ii =0; ii < N; ++ii) {
    for (int jj =0; jj < N; ++jj) {
      A[ii*N + jj] = ii + jj + 1;
      AT[ii*N + jj] = 0.0;
    }
  }
  AT = code_to_be_measured(A, N);

  std::cout<<"Begin test"<<std::endl;
  std::cout<<A[1*N+1]<<AT[1*N+1]<<std::endl;
  std::cout<<A[2*N+1]<<AT[1*N+2]<<std::endl;
  std::cout<<A[1*N+2]<<AT[2*N+1]<<std::endl;
  std::cout<<"End test"<<std::endl;  
  
  return 0;
}

std::vector<double> code_to_be_measured(const std::vector<double> M, int N)
{
  // matrix transpose with blocking
  std::vector<double> MT(N*N);
  for (int ii = 0; ii < N; ii+=csize)
    for (int jj = 0; jj < N; jj+=csize)
      for(int i=ii; i<std::min(N, ii+csize);++i)
	for(int j = jj; j < std::min(N,jj+csize);++j)
	  MT[i*N +j] = 2.0*M[j*N + i]; //Multiply by 2.0 so we can calculate floating point operations
  
  return MT;
}
