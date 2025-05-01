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
 * MathExtension.h
 *
 *  Created on: Sep. 22, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_MATHEXTENSION_H_
#define DINAMICA_MATHEXTENSION_H_

namespace DACE{

DACE_API double absolute(const double x);                        //!< Absolute value
DACE_API double cons(const double x);                            //!< Constant part (i.e. the value x)
DACE_API double logb(const double x, const double b = 10.0);     //!< Logarithm relative to base b
DACE_API double isrt(const double x);                            //!< Inverse square root
DACE_API double icbrt(const double x);                           //!< Inverse cube root
DACE_API double sqr(const double x);                             //!< Square
DACE_API double minv(const double x);                            //!< Multiplicative inverse
DACE_API double root(const double x, const int p = 2);           //!< p-th root
DACE_API double norm(const double x, const int type = 0);        //!< norm (same as abs)

}

#endif /* DINAMICA_MATHEXTENSION_H_ */
