/*
 * Matrix3x3_test.cpp
 *
 *  Created on: Dec 17, 2017
 *      Author: diego
 */

#include "../Sources/Matrix3x3.h"
#include "/home/diego/softwares/eigen3.3/Eigen/Dense"

namespace numeric {

void testMultiplyOperatorMatrix3x3andHermitian3x3(){
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

	Hermitian3x3 Brec;
	Brec = C*HA;

	real det;
	HB.inverseFast(&det);

	std::cout << "Hermitian3x3 HB: \n" << HB << std::endl;
	std::cout << "Recovered B: \n" << Brec << std::endl;

}

} /* namespace numeric */
