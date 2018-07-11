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
 * DAFormatter.h
 *
 *  Created on: Oct 18, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_DAFORMATTER_H_
#define DINAMICA_DAFORMATTER_H_

// C++ stdlib classes used in this public interface
#include <vector>
#include <string>

#include "DA/Def.h"

namespace DACE{

// forward declaration
class DA;

/*! Abstract class providing a DA formatter to output DA vectors in some advanced format. */
class DACE_API DAFormatter
{
public:
    virtual std::string format(const DA &da) = 0;
    virtual std::string format(const std::vector<DA> &da) = 0;
};

/*! Class containing the elements of a simple format as used by the DASimpleFormatter.
   \sa DASimpleFormatter
*/
struct DASimpleFormat {
    std::string pos, neg, mul, pre_pow, var, pre_var, post_var, pow, post_pow, linebreak;
    int first_var, first_pow;
    unsigned int monperline;
    bool shorten;
};

/*! DASimpleFormatter class which formats a DA vector using simple rules to output code suitable for various programming languages. */
class DACE_API DASimpleFormatter : public DAFormatter
{
public:
    static const DASimpleFormat C;
    static const DASimpleFormat C_POW;
    static const DASimpleFormat FORTRAN;
    static const DASimpleFormat FORTRAN_POW;
    static const DASimpleFormat MATLAB;
    static const DASimpleFormat MATLAB_POW;
    static const DASimpleFormat LATEX;

    DASimpleFormat sf;

    DASimpleFormatter() : sf(C) {};
    DASimpleFormatter(const DASimpleFormat& isf) : sf(isf) {};

    std::string format(const DA &da);
    std::string format(const std::vector<DA> &da);
};

}
#endif /* DINAMICA_DAFORMATTER_H_ */
