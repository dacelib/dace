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
 * Monomial.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Dinamica Srl
 */

// C++ stdlib classes used only internally in the implementation
#include <sstream>
#include <iomanip>

// DACE classes
#include "dace/Monomial.h"
#include "dace/DA.h"

namespace DACE{

/********************************************************************************
*     Constructors & Destructors
*********************************************************************************/
Monomial::Monomial() : m_jj(DA::getMaxVariables()), m_coeff(0.0) {
/*! Create a Monomial object large enough to hold all current monomials.
 */
}

unsigned int Monomial::order() const{
/*! Compute the order of the monomial.
   \return Order of the monomial
 */
    unsigned int ord = 0;

    for(unsigned int i = 0; i < m_jj.size(); i++)
        ord += m_jj[i];

    return ord;
}

std::string Monomial::toString() const{
/*! Convert monomial to string.
   \return A string representing the monomial in human readable form.
 */
    std::ostringstream oss;

    oss << "     I  COEFFICIENT              ORDER EXPONENTS" << std::endl;

    // value and order
    oss << "     1  ";
    oss << std::setiosflags(std::ios::uppercase) << std::setprecision(16) << std::scientific << std::setw(24) << m_coeff;
    oss << std::setw(4) << order() << std::setw(1) << ' ';

    // exponents
    for(unsigned int i = 0; i < m_jj.size(); i++){
        oss << std::setw(1) << ' ' << std::setw(2) << m_jj[i];}

    oss << std::endl;
    oss << "------------------------------------------------" << std::endl;

    return oss.str();
}

std::ostream& operator<< (std::ostream &out, const Monomial &m){
/*! Overload of std::operator<< in iostream.
   \param[in] out standard output stream.
   \param[in] m Monomial vector to be printed in the stream
   \return The output stream out.
 */
    out << m.toString();
    return out;
}

}
