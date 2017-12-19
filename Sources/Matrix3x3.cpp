/*
 * Matrix3x3.cpp
 *
 *  Created on: Dec 17, 2017
 *      Author: diego
 */

#include "Matrix3x3.h"

namespace numeric {

Matrix3x3::Matrix3x3() {
	// TODO Auto-generated constructor stub

}

Matrix3x3::~Matrix3x3() {
	// TODO Auto-generated destructor stub
}

Hermitian3x3 Matrix3x3::operator *(const Hermitian3x3& obj) {
	Hermitian3x3 H;

	//00
	H.setDiagonal(0, ((*this)(0,0)*obj.diagonal[0]+(*this)(0,1)*std::conj(obj.offdiagonal[0])+(*this)(0,2)*std::conj(obj.offdiagonal[1])).real());
	//01
	H.setOffDiagonal(0, (*this)(0,0)*obj.offdiagonal[0]+(*this)(0,1)*obj.diagonal[1]+(*this)(0,2)*std::conj(obj.offdiagonal[2]));
	//02
	H.setOffDiagonal(1, (*this)(0,0)*obj.offdiagonal[1]+(*this)(0,1)*obj.offdiagonal[2]+(*this)(0,2)*obj.diagonal[2]);
	//11
	H.setDiagonal(1, ((*this)(1,0)*(obj.offdiagonal[0])+(*this)(1,1)*obj.diagonal[1]+(*this)(1,2)*std::conj(obj.offdiagonal[2])).real());
	//12
	H.setOffDiagonal(2, ((*this)(1,0)*(obj.offdiagonal[1])+(*this)(1,1)*obj.offdiagonal[2]+(*this)(1,2)*obj.diagonal[2]));
	//22
	H.setDiagonal(2, ((*this)(2,0)*(obj.offdiagonal[1])+(*this)(2,1)*obj.offdiagonal[2]+(*this)(2,2)*obj.diagonal[2]).real());

	return H;
}

} /* namespace numeric */
