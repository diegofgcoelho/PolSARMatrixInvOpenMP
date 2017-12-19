/*
 * Sinclair.h
 *
 *  Created on: Dec 13, 2017
 *      Author: diego
 */

#ifndef HERMITIAN3X3_H_
#define HERMITIAN3X3_H_

#include <cmath>
#include <iostream>
//#include "support.h"
#include "numeric_type.h"
#include "Matrix3x3.h"


namespace numeric {

class Matrix3x3;

class Hermitian3x3 {
public:
	Hermitian3x3();
	Hermitian3x3(Hermitian3x3& A);
	Hermitian3x3(real a, real b, real c, complex ca, complex cb, complex cc);
	virtual ~Hermitian3x3();
	friend class Matrix3x3;
	inline real getDiagonal(unsigned int i){
		return this->diagonal[i];
	}
	inline void setDiagonal(unsigned int i, real a){
		this->diagonal[i] = a;
	}
	inline complex getOffDiagonal(unsigned int i){
		return this->offdiagonal[i];
	}
	inline void setOffDiagonal(unsigned int i, complex a){
		this->offdiagonal[i] = a;
	}
	/*
	 * This method computes only the matrix determinant using the cofactors.
	 */
	real detFast();
	/*
	 * This method computes the inverse matrix in place using Cholesky factorization.
	 * It also returns the determinant value in the address pointed by the input argument.
	 */
	void inverseCholesky(real* determinant);
	/*
	 * This method computes the inverse matrix in place using the proposed fast algorithm.
	 * It also returns the determinant value in the address pointed by the input argument.
	 */
	void inverseFast(real* determinant);

	friend std::ostream& operator<< (std::ostream& os, Hermitian3x3& obj){
		os << "\n[" << obj.diagonal[0] << "     " << obj.offdiagonal[0] << " " << obj.offdiagonal[1] << "\n";
		os << "x.xx    " << obj.diagonal[1] << "   " << obj.offdiagonal[2] << "\n";
		os << "x.xx    " << "x.xx       " << obj.diagonal[2] << "]\n";
		return os;
	}

	Matrix3x3 operator*(const Hermitian3x3& obj);
	Hermitian3x3 operator*(const Matrix3x3& obj);

	Hermitian3x3& operator= (const Hermitian3x3& obj){
			if(this != &obj){
				for(unsigned int i = 0; i < 3; i++){
					this->diagonal[i] = obj.diagonal[i];
					this->offdiagonal[i] = obj.offdiagonal[i];
				}
			}
			return *this;
		}
private:
	/*
	 * The elements of the diagonal array represents, respectively,
	 * the elements along the main diagonal: (0,0), (1,1), and (2,2).
	 */
	real diagonal[3];
	/*
	 * The elements of the offdiagonal array represents, respectively,
	 * the elements off the diagonal in the following order: (0,1), (0,2), and (1,2).
	 */
	complex offdiagonal[3];
};

//Test functions
void testprint();
void testCholesky();
void testInverseLCholesky();
void testInverseCholesky();
void testInverseFast();
void testMultiplyOperatorHermitian3x3andHermitian3x3();
void testCopyOperator();
void testMultiplyOperatorHermitian3x3andMatrix3x3();
void testDetFast();
} /* namespace numeric */

#endif /* HERMITIAN3X3_H_ */
