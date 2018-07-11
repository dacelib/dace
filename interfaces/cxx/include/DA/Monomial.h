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
 * Monomial.h
 *
 *  Created on: Mar 10, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_MONOMIAL_H_
#define DINAMICA_MONOMIAL_H_

// C++ stdlib classes used in this public interface
#include <vector>
#include <string>
#include <ostream>

#include "DA/Def.h"

namespace DACE{

/*! Monomial class */
class DACE_API Monomial
{
public:
    std::vector<unsigned int> m_jj;     /*!< Vector of exponents.               */
    double m_coeff;                     /*!< Coefficient.                       */

    Monomial();                         /*!< Default constructor.               */

    unsigned int order() const;         /*!< Return the order of the Monomial.  */
    std::string toString() const;       /*!< Convert current monomial to string.*/
};

DACE_API std::ostream& operator<< (std::ostream &out, const Monomial &m);    /*!< Overload output stream operator. */

}

#endif /* DINAMICA_MONOMIAL_H_ */
