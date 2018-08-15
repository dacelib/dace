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
 * compiledDA.h
 *
 *  Created on: Mar 01, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_COMPILEDDA_H_
#define DINAMICA_COMPILEDDA_H_

// C++ stdlib classes used in this public interface
#include <vector>

namespace DACE{

class DA;   // forward declaration

/*! compiledDA class representing a precomputed representation of a polynomial for efficient evaluation. */
class DACE_API compiledDA
{
private:
    double *ac;             //!< Compiled polynomial evaluation data
    unsigned int dim;       //!< Number of polynomials (dimension)
    unsigned int ord;       //!< Maximum order of the polynomial
    unsigned int vars;      //!< Number of variables in the polynomial
    unsigned int terms;     //!< Number of terms in the polynomial

public:
    /********************************************************************************
    *     Constructors & Destructors
    *********************************************************************************/
    compiledDA(const compiledDA &cda);                      //!< Copy constructor.
    compiledDA(const DA &da);                               //!< Constructor from a single DA.
    compiledDA(const std::vector<DA> &da);                  //!< Constructor from a vector of DA.
    ~compiledDA() throw();                                  //!< Default destructor.

    /********************************************************************************
    *     Assignments
    *********************************************************************************/
    compiledDA& operator=(const compiledDA &cda);           //!< Copy assignment.

    /********************************************************************************
    *     Evaluation
    *********************************************************************************/
    template<class V> V eval(const V &args) const;                                          //!< Evaluate the compiled polynomial with a vector of any arithmetic type and return vector of results.
    template<class T> std::vector<T> eval(const T args[], const unsigned int length) const; //!< Evaluate the compiled polynomial with an array of any arithmetic type and return vector of results.
    template<class T> std::vector<T> evalScalar(const T &arg) const;                        //!< Evaluate the compiled polynomial with a single variable of arithmetic type and return vector of results.
    template<class T> void eval(const std::vector<T> &args, std::vector<T> &res) const;     //!< Evaluate the compiled polynomial with a vector of any arithmetic type and store results in provided vector.

    /********************************************************************************
    *     Member access routines
    *********************************************************************************/
    const double* getAc() const;                            //!< Return a pointer to the internal coefficient array
    unsigned int getDim() const;                            //!< Return the number of polynomials
    unsigned int getOrd() const;                            //!< Return the maximum order in all polynomials
    unsigned int getVars() const;                           //!< Return the number of independent variables in all polynomials
    unsigned int getTerms() const;                          //!< Return the number of terms
};

// specializations for particularly efficient evaluation with double and DA arguments implemented in the library
template<> DACE_API void compiledDA::eval(const std::vector<DA> &args, std::vector<DA> &res) const;
template<> DACE_API void compiledDA::eval(const std::vector<double> &args, std::vector<double> &res) const;

}
#endif /* DINAMICA_COMPILEDDA_H_ */
