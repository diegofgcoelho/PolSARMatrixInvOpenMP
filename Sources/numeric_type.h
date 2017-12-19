/*
 * numeric_type.h
 *
 *  Created on: Dec 13, 2017
 *      Author: diego
 */

#ifndef NUMERIC_TYPE_H_
#define NUMERIC_TYPE_H_

#include <complex>
#include "/home/diego/softwares/eigen3.3/Eigen/Dense"

namespace numeric {

/*
 * These are masks for the numeric type that we want to use in the simulations.
 */

typedef float real;

typedef std::complex<real> complex;

}  // namespace numeric

#endif /* NUMERIC_TYPE_H_ */
