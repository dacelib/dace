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
 * DA_t.h
 *
 *  Created on: Apr 07, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_DA_T_H_
#define DINAMICA_DA_T_H_

// DACE classes
#include "DA/compiledDA.h"
#include "DA/DA.h"

namespace DACE{

/********************************************************************************
*     DA polynomial evaluation routines
*********************************************************************************/
template<class T> T DA::eval(const std::vector<T> &args) const{
/*! Generic evaluation of the DA with a vector of arithmetic type T arguments.
   \param[in] args std::vector<T> of arithmetic type T with which the DA vector is evaluated
   \return The result of the evaluation
   \throw DACE::DACEException
   \note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
   \sa compiledDA::eval
 */
    compiledDA cda(*this);
    return cda.eval(args)[0];
}

template<class T> T DA::eval(const T args[], const unsigned int length) const{
/*! Generic evaluation of the DA with an array of arithmetic type T arguments.
   \param[in] args array of arithmetic type T with which the DA vector is evaluated.
   \param[in] length number of elements in the array args.
   \return The result of the evaluation.
   \throw DACE::DACEException
   \note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
   \sa compiledDA::eval
 */
    compiledDA cda(*this);
    return cda.eval(args,length)[0];
}

template<class T> T DA::evalScalar(const T &arg) const{
/*! Generic evaluation of the DA with a single arithmetic type T argument.
   \param[in] arg single variable of arithmetic type T of the first independent DA variable.
   \return The result of the evaluation.
   \note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
   \sa compiledDA::evalScalar
 */
    compiledDA cda(*this);
    return cda.evalScalar(arg)[0];
}

template<class T> T eval(const DA &da, const std::vector<T> &args) {
/*! Generic evaluation of the DA with a vector of arithmetic type T arguments.
   \param[in] da a given DA object.
   \param[in] args std::vector<T> of arithmetic type T with which the DA vector is evaluated.
   \return The result of the evaluation.
   \throw DACE::DACEException
   \note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
   \sa compiledDA
 */
    return da.eval(args);
}

template<class T> T eval(const DA &da, const T args[], const unsigned int length) {
/*! Generic evaluation of the DA with an array of arithmetic type T arguments.
   \param[in] da a given DA object.
   \param[in] args array of arithmetic type T with which the DA vector is evaluated.
   \param[in] length number of elements in the array args.
   \return The result of the evaluation.
   \throw DACE::DACEException
   \note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
   \sa compiledDA::eval
 */
    return da.eval(args,length);
}

template<class T> T evalScalar(const DA &da, const T &arg) {
/*! Generic evaluation of the DA with a single arithmetic type T arguments.
   \param[in] da a given DA object.
   \param[in] arg single variable of arithmetic type T of the first independent DA variable.
   \return The result of the evaluation.
   \throw DACE::DACEException
   \note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
   \sa compiledDA::evalScalar
 */
    return da.evalScalar(arg);
}

}
#endif /* DINAMICA_DA_T_H_ */
