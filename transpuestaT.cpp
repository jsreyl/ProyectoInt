#include <iostream>
#include <cstdio>
#include <cstdlib>
#include<vector>

std::vector<double> code_to_be_measured(const std::vector<double> M, int N);

int main(int argc, char **argv)
{
  int N = atoi(argv[1]); //Run N in a loop in bash
  //As int for M in 2 4 8 16 32..16384; do ./a.out $M; done
  // Matrix declaration : Modeled as 1D array
  // Declare as pointers and ask for memory to use the heap
  std::vector<double>A(N*N), AT(N*N);
  
  // initialize matrices
  for (int ii =0; ii < N; ++ii) {
    for (int jj =0; jj < N; ++jj) {
      A[ii*N + jj] = ii + 1;
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
  // simple matrix transpose
  std::vector<double> MT(N*N);
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      MT[ii*N +jj] = 2.0*M[jj*N + ii]; //Multiply by 2.0 so we can calculate floating point operations
    }
  }
  return MT;
}
