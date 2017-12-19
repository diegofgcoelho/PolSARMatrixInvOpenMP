/*
 * Sinclair_test.cpp
 *
 *  Created on: Dec 13, 2017
 *      Author: diego
 */

#include <iostream>
#include "../Sources/Hermitian3x3.h"
#include "/home/diego/softwares/eigen3.3/Eigen/Dense"

namespace numeric {

void testprint(){
	std::cout << "*****************************Testing the << operator*****************************" << std::endl;
	real a = 1.0, d = 2.3, f = 5;
	complex b = std::complex<real>(3.1,2), c = std::complex<real>(-1.1,3.5), e = std::complex<real>(-.1,0.7);
	Hermitian3x3 H(a, d, f, b, c, e);

	std::cout << "Printing Hermitian matrix" << std::endl;
	std::cout << H << std::endl;
}

void testCholesky(){
	std::cout << "*****************************Testing the Cholesky Factorization*****************************" << std::endl;
	srand((unsigned int) time(0));
	Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(3,3); A = A*A.conjugate().transpose();
	Eigen::LLT<Eigen::MatrixXcd> lltOfA(A);
	Eigen::MatrixXcd lA = lltOfA.matrixL();

	real a = A(0,0).real(), d = A(1,1).real(), f = A(2,2).real();
	complex b = std::complex<real>(A(0,1)), c =  std::complex<real>(A(0,2)), e = std::complex<real>(A(1,2));
	//Creating Hermitian3x3 object
	Hermitian3x3 H(a, d, f, b, c, e);

	real determinant = 0;

	std::cout << "Eigen matrix: \n" << A << std::endl;
	std::cout << "Hermitian3x3 matrix: " << H << std::endl;

	H.inverseCholesky(&determinant);

	std::cout << "Eigen Cholesky: \n" << lA << std::endl;
}

void testInverseLCholesky(){
	std::cout << "*****************************Testing the Inverse of Cholesky L Factor*****************************" << std::endl;
	srand((unsigned int) time(0));
	Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(3,3); A = A*A.conjugate().transpose();
	Eigen::LLT<Eigen::MatrixXcd> lltOfA(A);
	Eigen::MatrixXcd lA = lltOfA.matrixL();

	real a = A(0,0).real(), d = A(1,1).real(), f = A(2,2).real();
	complex b = std::complex<real>(A(0,1)), c =  std::complex<real>(A(0,2)), e = std::complex<real>(A(1,2));
	//Creating Hermitian3x3 object
	Hermitian3x3 H(a, d, f, b, c, e);

	real determinant = 0;

	std::cout << "Eigen matrix: \n" << A << std::endl;
	std::cout << "Hermitian3x3 matrix: " << H << std::endl;

	H.inverseCholesky(&determinant);

	std::cout << "Eigen L Inverse: \n" << lA.inverse() << std::endl;
}

void testInverseCholesky(){
	std::cout << "*****************************Testing the Inverse based on Cholesky*****************************" << std::endl;
	srand((unsigned int) time(0));
	Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(3,3); A = A*A.conjugate().transpose();

	real a = A(0,0).real(), d = A(1,1).real(), f = A(2,2).real();
	complex b = std::complex<real>(A(0,1)), c =  std::complex<real>(A(0,2)), e = std::complex<real>(A(1,2));
	//Creating Hermitian3x3 object
	Hermitian3x3 H(a, d, f, b, c, e);

	real determinant = 0;

	std::cout << "Eigen matrix: \n" << A << std::endl;
	std::cout << "Hermitian3x3 matrix: " << H << std::endl;

	H.inverseCholesky(&determinant);

	std::cout << "Eigen Inverse: \n" << A.inverse() << std::endl;
	std::cout << "Eigen Determinant: " << A.determinant() << std::endl;
	std::cout << "Hermitian3x3 Cholesky Inverse: \n" << H << std::endl;
	std::cout << "Hermitian3x3 Cholesky Determinant: " << determinant << std::endl;
}

void testInverseFast(){
	std::cout << "*****************************Testing the Fast Inversion*****************************" << std::endl;
	srand((unsigned int) time(0));
	Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(3,3); A = A*A.conjugate().transpose();

	real a = A(0,0).real(), d = A(1,1).real(), f = A(2,2).real();
	complex b = std::complex<real>(A(0,1)), c =  std::complex<real>(A(0,2)), e = std::complex<real>(A(1,2));
	//Creating Hermitian3x3 object
	Hermitian3x3 H(a, d, f, b, c, e);

	real determinant = 0;

	std::cout << "Eigen matrix: \n" << A << std::endl;
	std::cout << "Hermitian3x3 matrix: " << H << std::endl;

	H.inverseFast(&determinant);

	std::cout << "Eigen Inverse: \n" << A.inverse() << std::endl;
	std::cout << "Eigen Determinant: " << A.determinant() << std::endl;
	std::cout << "\nHermitian3x3 Fast Inverse: " << H << std::endl;
	std::cout << "Hermitian3x3 Fast Determinant: " << determinant << std::endl;
}

void testMultiplyOperatorHermitian3x3andHermitian3x3(){
	std::cout << "*****************************Testing the * operator*****************************" << std::endl;
	srand((unsigned int) time(0));
	Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(3,3); A = A*A.conjugate().transpose();

	real a = A(0,0).real(), d = A(1,1).real(), f = A(2,2).real();
	complex b = std::complex<real>(A(0,1)), c =  std::complex<real>(A(0,2)), e = std::complex<real>(A(1,2));

	Hermitian3x3 HA(a, d, f, b, c, e);

	Eigen::MatrixXcd B = Eigen::MatrixXcd::Random(3,3); B = B*B.conjugate().transpose();

	a = B(0,0).real(); d = B(1,1).real(); f = B(2,2).real();
	b = std::complex<real>(B(0,1)); c =  std::complex<real>(B(0,2)); e = std::complex<real>(B(1,2));

	//Creating Hermitian3x3 object
	Hermitian3x3 HB(a, d, f, b, c, e);

	std::cout << "Eigen matrix A: \n" << A << std::endl;
	std::cout << "Eigen matrix B: \n" << B << std::endl;
	std::cout << "Eigen matrix A*B: \n" << A*B << std::endl;

	std::cout << "Hermitian3x3 matrix A: \n" << HA << std::endl;
	std::cout << "Hermitian3x3 matrix B: \n" << HB << std::endl;
	std::cout << "Hermitian3x3 matrix A*B: \n" << HA*HB << std::endl;

}

void testCopyOperator(){
	std::cout << "*****************************Testing the = operator*****************************" << std::endl;
	srand((unsigned int) time(0));
	Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(3,3); A = A*A.conjugate().transpose();

	real a = A(0,0).real(), d = A(1,1).real(), f = A(2,2).real();
	complex b = std::complex<real>(A(0,1)), c =  std::complex<real>(A(0,2)), e = std::complex<real>(A(1,2));

	Hermitian3x3 HA(a, d, f, b, c, e);
	Hermitian3x3 HB;

	std::cout << "Eigen matrix A: \n" << A << std::endl;

	std::cout << "Hermitian3x3 matrix A: \n" << HA << std::endl;
	std::cout << "Hermitian3x3 matrix B before copy (=): \n" << HB << std::endl;
	HB = HA;
	std::cout << "Hermitian3x3 matrix B after copy (=): \n" << HB << std::endl;


}

void testMultiplyOperatorHermitian3x3andMatrix3x3() {
	std::cout << "*****************************Testing the * operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(3,3); A = A*A.conjugate().transpose();
	real a = A(0,0).real(), d = A(1,1).real(), f = A(2,2).real();
	complex b = std::complex<real>(A(0,1)), c =  std::complex<real>(A(0,2)), e = std::complex<real>(A(1,2));
	Hermitian3x3 HA(a, d, f, b, c, e);

	Eigen::MatrixXcd B = Eigen::MatrixXcd::Random(3,3); B = B*B.conjugate().transpose();
	a = B(0,0).real(); d = B(1,1).real(); f = B(2,2).real();
	b = std::complex<real>(B(0,1)); c =  std::complex<real>(B(0,2)); e = std::complex<real>(B(1,2));
	Hermitian3x3 HB(a, d, f, b, c, e);

	Matrix3x3 C = HA*HB;
	C = C.inverse();

	Hermitian3x3 Arec;
	Arec = HB*C;

	real det;
	HA.inverseFast(&det);

	std::cout << "Hermitian3x3 HA: \n" << HA << std::endl;
	std::cout << "Recovered A: \n" << Arec << std::endl;

//	for (int var = 0; var < 3; ++var) {
//		Eigen::Matrix<numeric::complex, 3, 3, Eigen::RowMajor> A = Eigen::Matrix<numeric::complex, 3, 3, Eigen::RowMajor>::Random(3,3); A = A*A.conjugate().transpose();
//		Eigen::Matrix<numeric::complex, 3, 3, Eigen::RowMajor> B = Eigen::Matrix<numeric::complex, 3, 3, Eigen::RowMajor>::Random(3,3); B = B*B.conjugate().transpose();
//		C = A*B;
//
//
//		Matrix3x3 D; D = C.inverse();
//		complex dett; bool flag;
//		Matrix3x3 E; D.computeInverseAndDetWithCheck(E, dett, flag);
//		D.computeInverseAndDetWithCheck(D, dett, flag);
//
//		std::cout << "C: \n" << C << std::endl;
//		std::cout << "D: \n" << D << std::endl;
//		std::cout << "E: \n" << E << std::endl;
//		std::cout << "det = " << dett << "\n\n\n\n" << std::endl;
//	}

}

void testDetFast(){
	std::cout << "*****************************Testing the detFast*****************************" << std::endl;
	srand((unsigned int) time(0));

	Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(3,3); A = A*A.conjugate().transpose();
	real a = A(0,0).real(), d = A(1,1).real(), f = A(2,2).real();
	complex b = std::complex<real>(A(0,1)), c =  std::complex<real>(A(0,2)), e = std::complex<real>(A(1,2));
	Hermitian3x3 HA(a, d, f, b, c, e);

	std::cout << "Eigen determinant: " << A.determinant().real() << std::endl;
	std::cout << "Hermitian3x3 determinant: " << HA.detFast() << std::endl;

}


} /* namespace numeric */
