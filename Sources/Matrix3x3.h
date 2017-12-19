/*
 * Matrix3x3.h
 *
 *  Created on: Dec 17, 2017
 *      Author: diego
 */

#ifndef MATRIX3X3_H_
#define MATRIX3X3_H_

//#include "support.h"
#include "numeric_type.h"
#include "Hermitian3x3.h"


namespace numeric {

class Hermitian3x3;

class Matrix3x3:public Eigen::Matrix<complex, 3, 3, Eigen::RowMajor> {
public:
	friend class Hermitian3x3;
	Matrix3x3();
	virtual ~Matrix3x3();
	/*
	 * Defining the * operator that will denote matrix multiplication of a
	 * Matrix3x3 by a Hermitian3x3 matrix.
	 */
	Hermitian3x3 operator*(const Hermitian3x3& obj);
	/*
	 * Defining the = operators because they are not inhereted. These operators will
	 * be used for copying the rhs, either it is Eigen::Matrix or Matrix3x3. In the first case,
	 * the input Eigen::Matrix is assumed to have size 3x3, exactly.
	 */
	inline Matrix3x3 operator=(const Eigen::Matrix<complex, 3, 3, Eigen::RowMajor>& obj){
		for(unsigned int i = 0; i < 3; i++){
			for(unsigned int j = 0; j < 3; j++){
				(*this)(i,j) = static_cast<complex>(obj(i,j));
			}
		}
		return *this;
	}
	inline Matrix3x3 operator=(const Matrix3x3& obj){
		for(unsigned int i = 0; i < 3; i++){
			for(unsigned int j = 0; j < 3; j++){
				(*this)(i,j) = obj(i,j);
			}
		}
		return *this;
	}

};

void testMultiplyOperatorMatrix3x3andHermitian3x3();

} /* namespace numeric */

#endif /* MATRIX3X3_H_ */
