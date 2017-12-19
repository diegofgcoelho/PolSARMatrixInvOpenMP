/*
 * Sinclair.cpp
 *
 *  Created on: Dec 13, 2017
 *      Author: diego
 */

#include "Hermitian3x3.h"

namespace numeric {

Hermitian3x3::Hermitian3x3() {
	for(unsigned int i = 0; i < 3; i++){
		this->diagonal[i] = 0.0;
		this->offdiagonal[i] = std::complex<real>(0.0, 0.0);
	}
}

Hermitian3x3::Hermitian3x3(Hermitian3x3& A) {
	for(unsigned int i = 0; i < 3; i++){
		this->diagonal[i] = A.getDiagonal(i);
		this->offdiagonal[i] = A.getOffDiagonal(i);
	}
}

Hermitian3x3::Hermitian3x3(real a, real b, real c, complex ca, complex cb, complex cc) {
	this->diagonal[0] = a;
	this->diagonal[1] = b;
	this->diagonal[2] = c;
	this->offdiagonal[0] = ca;
	this->offdiagonal[1] = cb;
	this->offdiagonal[2] = cc;
}

Hermitian3x3::~Hermitian3x3() {
}

real Hermitian3x3::detFast(){
	real determinant;
	//compute the cofactors
	//a
	real a = this->diagonal[1]*this->diagonal[2]-std::norm(this->offdiagonal[2]);
	//b
	complex b = this->offdiagonal[1]*std::conj(this->offdiagonal[2])-this->offdiagonal[0]*this->diagonal[2];
	//c
	complex c = this->offdiagonal[0]*this->offdiagonal[2]-this->offdiagonal[1]*this->diagonal[1];

	//Compute determinant
	determinant = this->diagonal[0]*a+this->offdiagonal[0].real()*b.real()+
			this->offdiagonal[0].imag()*b.imag()+this->offdiagonal[1].real()*c.real()+
			this->offdiagonal[1].imag()*c.imag();

	return determinant;
}

void Hermitian3x3::inverseCholesky(real* determinant) {
	//Compute Cholesky factorization of A = L^TL
	//l00
	this->diagonal[0] = sqrt(this->diagonal[0]);
	//v0
	real v0 = 1/this->diagonal[0];
	//l10
	this->offdiagonal[0] = std::conj(this->offdiagonal[0])*v0;
	//l20
	this->offdiagonal[1] = std::conj(this->offdiagonal[1])*v0;
	//l11
	this->diagonal[1] = sqrt(this->diagonal[1]-std::norm(this->offdiagonal[0]));
	//v1
	real v1 = 1/this->diagonal[1];
	//l21
	this->offdiagonal[2] = std::conj((this->offdiagonal[2]-this->offdiagonal[0]*std::conj(this->offdiagonal[1]))*v1);
	//l22
	this->diagonal[2] = sqrt(this->diagonal[2]-std::norm(this->offdiagonal[1])-std::norm(this->offdiagonal[2]));
	//v2
	real v2 = 1/this->diagonal[2];

	//std::cout << "Cholesky factorization is\n" << *this << std::endl;

	//Compute the determinant
	*determinant = pow(this->diagonal[0]*this->diagonal[1]*this->diagonal[2],2);

	//Compute the elements of L factor inverse, [L^-1]_ij = tij
	//t00
	this->diagonal[0] = v0;
	//t11
	this->diagonal[1] = v1;
	//t20
	this->offdiagonal[1] = v2*v0*(this->offdiagonal[2]*this->offdiagonal[0]*v1-this->offdiagonal[1]);
	//t10
	this->offdiagonal[0] = -v1*this->offdiagonal[0]*v0;
	//this->offdiagonal[0].real(-v1*this->offdiagonal[0].real()*v0);
	//this->offdiagonal[0].imag(-v1*this->offdiagonal[0].imag()*v0);
	//t21
	this->offdiagonal[2] = -v1*v2*this->offdiagonal[2];
	//t22
	this->diagonal[2] = v2;

	//std::cout << "Inverse of L factor is\n" << *this << std::endl;

	//Compute the elements of A^-1, [A^-1]_ij = aij
	//a00
	this->diagonal[0] = this->diagonal[0]*this->diagonal[0]+std::norm(this->offdiagonal[0])+std::norm(this->offdiagonal[1]);
	//a01
	this->offdiagonal[0] = std::conj(this->offdiagonal[0])*this->diagonal[1]+std::conj(this->offdiagonal[1])*this->offdiagonal[2];
	//a02
	this->offdiagonal[1] = std::conj(this->offdiagonal[1])*this->diagonal[2];
	//a11
	this->diagonal[1] = this->diagonal[1]*this->diagonal[1]+std::norm(this->offdiagonal[2]);
	//a12
	this->offdiagonal[2] = std::conj(this->offdiagonal[2])*this->diagonal[2];
	//a22
	this->diagonal[2] = this->diagonal[2]*this->diagonal[2];
}

void Hermitian3x3::inverseFast(real* determinant) {
	//compute the cofactors
	//a
	real a = this->diagonal[1]*this->diagonal[2]-std::norm(this->offdiagonal[2]);
	//b
	complex b = this->offdiagonal[1]*std::conj(this->offdiagonal[2])-this->offdiagonal[0]*this->diagonal[2];
	//c
	complex c = this->offdiagonal[0]*this->offdiagonal[2]-this->offdiagonal[1]*this->diagonal[1];
	//d
	real d = this->diagonal[0]*this->diagonal[2]-std::norm(this->offdiagonal[1]);
	//e
	complex e = this->offdiagonal[1]*std::conj(this->offdiagonal[0])-this->diagonal[0]*this->offdiagonal[2];
	//f
	real f = this->diagonal[0]*this->diagonal[1]-std::norm(this->offdiagonal[0]);

	//Compute determinant
	*determinant = this->diagonal[0]*a+this->offdiagonal[0].real()*b.real()+
			this->offdiagonal[0].imag()*b.imag()+this->offdiagonal[1].real()*c.real()+
			this->offdiagonal[1].imag()*c.imag();

	//Compute determinant inverse
	real t = 1/(*determinant);

	//Compute the elements of the inverse
	this->diagonal[0] = t*a;
	this->offdiagonal[0] = t*b;
	this->offdiagonal[1] = t*c;
	this->diagonal[1] = t*d;
	this->offdiagonal[2] = t*e;
	this->diagonal[2] = t*f;
}

Matrix3x3 Hermitian3x3::operator* (const Hermitian3x3& obj) {
	complex output[9];
	//Compute the elements at position ij
	//00
	output[0] = this->diagonal[0]*obj.diagonal[0]+
			this->offdiagonal[0]*std::conj(obj.offdiagonal[0])+
			this->offdiagonal[1]*std::conj(obj.offdiagonal[1]);
	//01
	output[1] = this->diagonal[0]*obj.offdiagonal[0]+
			this->offdiagonal[0]*obj.diagonal[1]+
			this->offdiagonal[1]*std::conj(obj.offdiagonal[2]);
	//02
	output[2] = this->diagonal[0]*obj.offdiagonal[1]+
			this->offdiagonal[0]*obj.offdiagonal[2]+
			this->offdiagonal[1]*obj.diagonal[2];
	//10
	output[3] = std::conj(this->offdiagonal[0])*obj.diagonal[0]+
			this->diagonal[1]*std::conj(obj.offdiagonal[0])+
			this->offdiagonal[2]*std::conj(obj.offdiagonal[1]);
	//11
	output[4] = std::conj(this->offdiagonal[0])*obj.offdiagonal[0]+
			this->diagonal[1]*obj.diagonal[1]+
			this->offdiagonal[2]*std::conj(obj.offdiagonal[2]);
	//12
	output[5] = std::conj(this->offdiagonal[0])*obj.offdiagonal[1]+
				this->diagonal[1]*obj.offdiagonal[2]+
				this->offdiagonal[2]*obj.diagonal[2];
	//20
	output[6] = std::conj(this->offdiagonal[1])*obj.diagonal[0]+
				std::conj(this->offdiagonal[2])*std::conj(obj.offdiagonal[0])+
				this->diagonal[2]*std::conj(obj.offdiagonal[1]);
	//21
	output[7] = std::conj(this->offdiagonal[1])*obj.offdiagonal[0]+
				std::conj(this->offdiagonal[2])*obj.diagonal[1]+
				this->diagonal[2]*std::conj(obj.offdiagonal[2]);
	//22
	output[8] = std::conj(this->offdiagonal[1])*obj.offdiagonal[1]+
					std::conj(this->offdiagonal[2])*obj.offdiagonal[2]+
					this->diagonal[2]*obj.diagonal[2];
	Matrix3x3 Mat;
	Mat << output[0], output[1], output[2], output[3], output[4], output[5], output[6], output[7], output[8];
	return Mat;
}

Hermitian3x3 Hermitian3x3::operator*(const Matrix3x3& obj){

	Hermitian3x3 H;
	//Compute the elements at position ij
	//00
	H.diagonal[0] = (this->diagonal[0]*obj(0,0)+this->offdiagonal[0]*obj(1,0)+this->offdiagonal[1]*obj(2,0)).real();
	//01
	H.offdiagonal[0] = this->diagonal[0]*obj(0,1)+this->offdiagonal[0]*obj(1,1)+this->offdiagonal[1]*obj(2,1);
	//02
	H.offdiagonal[1] = this->diagonal[0]*obj(0,2)+this->offdiagonal[0]*obj(1,2)+this->offdiagonal[1]*obj(2,2);
	//11
	H.diagonal[1] = (std::conj(this->offdiagonal[0])*obj(0,1)+this->diagonal[1]*obj(1,1)+this->offdiagonal[2]*obj(2,1)).real();
	//12
	H.offdiagonal[2] = std::conj(this->offdiagonal[0])*obj(0,2)+this->diagonal[1]*obj(1,2)+this->offdiagonal[2]*obj(2,2);
	//22
	H.diagonal[2] = (std::conj(this->offdiagonal[1])*obj(0,2)+std::conj(this->offdiagonal[2])*obj(1,2)+this->diagonal[2]*obj(2,2)).real();

	return H;
}

} /* namespace numeric */
