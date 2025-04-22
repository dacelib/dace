/******************************************************************************
*                                                                             *
* DIFFERENTIAL ALGEBRA CORE ENGINE                                            *
*                                                                             *
*******************************************************************************
*                                                                             *
* Copyright 2016 Politecnico di Milano (2014 Dinamica Srl)                    *
* Licensed under the Apache License, Version 2.0 (the "License");             *
* you may not use this file except in compliance with the License.            *
* You may obtain a copy of the License at                                     *
*                                                                             *
*    http://www.apache.org/licenses/LICENSE-2.0                               *
*                                                                             *
* Unless required by applicable law or agreed to in writing, software         *
* distributed under the License is distributed on an "AS IS" BASIS,           *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    *
* See the License for the specific language governing permissions and         *
* limitations under the License.                                              *
*                                                                             *
*******************************************************************************/

/*
 * MathExtension.cpp
 *
 *  Created on: Sep. 22, 2014
 *      Author: Dinamica Srl
 */

// C++ stdlib classes
#include <cmath>

// DACE classes
#include "dace/config.h"
#include "dace/MathExtension.h"

namespace DACE{

double absolute(const double x){
/*! Absolute value.
   \param[in] x Function argument.
 */
    return std::abs(x);
}

double cons(const double x){
/*! Constant part. For double type this is just x.
   \param[in] x Function argument.
 */
    return x;
}

double logb(const double x, const double b){
/*! Logarithm relative to base b.
   \param[in] x Function argument.
   \param[in] b Base of the logarithm (must be positive).
 */
    return std::log(x)/std::log(b);
}

double isrt(const double x){
/*! Inverse square root 1/sqrt(x).
   \param[in] x Function argument.
 */
    return 1.0/std::sqrt(x);
}

double sqr(const double x){
/*! Square of x.
   \param[in] x Function argument.
 */
    return x*x;
}

double minv(const double x){
/*! Multiplicative inverse 1/x.
   \param[in] x Function argument.
 */
    return 1.0/x;
}

double root(const double x, const int p){
/*! p-th root of x.
   \param[in] x Function argument.
   \param[in] p Root to take.
 */
    return std::pow(x, 1.0/p);
}

}
