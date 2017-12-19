/*
 * main.cpp
 *
 *  Created on: Dec 14, 2017
 *      Author: diego
 */

#include <iostream>
#include <omp.h>

#include "numeric_type.h"
#include "Hermitian3x3.h"
#include "Matrix3x3.h"

#define EIGEN_DONT_PARALLELIZE

//Number of replicates
const unsigned int  replicates = 10000;
//Number of matrices in each replicate
const unsigned int N = 1000000;

int main(int argc, char* argv[]){

	double cholesky_time_omp = 0.0, fast_time_omp = 0.0;

	//Random generator seed
	srand((unsigned int) time(0));

	//Allocate memory for original matrix and inverses
	numeric::Hermitian3x3* matrices = new numeric::Hermitian3x3[N];
	numeric::Hermitian3x3* matricesCholesky = new numeric::Hermitian3x3[N];
	numeric::Hermitian3x3* matricesFast = new numeric::Hermitian3x3[N];

	//Allocate memory for the determinants
	numeric::real* detCholesky = new numeric::real[N];
	numeric::real* detFast = new numeric::real[N];

	for (unsigned int rep = 0; rep < replicates; rep++) {
		//Generate matrices
		for(unsigned int n = 0; n < N; n++){
			Eigen::Matrix<numeric::complex, 3, 3, Eigen::RowMajor> A = Eigen::Matrix<numeric::complex, 3, 3, Eigen::RowMajor>::Random(3,3); A = A*A.conjugate().transpose();

			numeric::real a = A(0,0).real(), d = A(1,1).real(), f = A(2,2).real();
			numeric::complex b = std::complex<numeric::real>(A(0,1)), c =  std::complex<numeric::real>(A(0,2)), e = std::complex<numeric::real>(A(1,2));

			matrices[n].setDiagonal(0, a);
			matrices[n].setDiagonal(1, d);
			matrices[n].setDiagonal(2, f);

			matrices[n].setOffDiagonal(0,b);
			matrices[n].setOffDiagonal(1,c);
			matrices[n].setOffDiagonal(2,e);

			matricesCholesky[n] = matrices[n];
			matricesFast[n] = matrices[n];
		}


		//Perform Cholesky based inversion
		cholesky_time_omp = omp_get_wtime();
		#pragma omp parallel
		{
			#pragma omp for
			for(unsigned int n = 0; n < N; n++){
				matricesCholesky[n].inverseCholesky(&detCholesky[n]);
			}
		}
		cholesky_time_omp = omp_get_wtime()-cholesky_time_omp;

		//Perform Cholesky based inversion
		fast_time_omp = omp_get_wtime();
		#pragma omp parallel
		{
			#pragma omp for
			for(unsigned int n = 0; n < N; n++){
				matricesFast[n].inverseFast(&detFast[n]);
			}
		}
		fast_time_omp = omp_get_wtime()-fast_time_omp;
	}

	std::cout << "The times with OMP are (ms): " << std::endl;
	std::cout << "Cholesky: " << 1000*cholesky_time_omp/replicates << "." << std::endl;
	std::cout << "Fast: " << 1000*fast_time_omp/replicates << "." << std::endl;

	delete[] matrices;
	delete[] matricesCholesky;
	delete[] matricesFast;

	return 0;
}


