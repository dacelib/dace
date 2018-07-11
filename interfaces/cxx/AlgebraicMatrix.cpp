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
 * AlgebraicMatrix.cpp
 *
 *  Created on: January 21, 2015
 *      Author: Dinamica Srl
 */

// C++ stdlib classes used only internally in this implementation
#include <string>
#include <sstream>

// DACE classes
#include "DA/DA.h"
#include "DA/AlgebraicMatrix.h"
#include "DA/AlgebraicMatrix_t.h"

namespace DACE {

/***********************************************************************************
*     Input/Output routines
************************************************************************************/
template<> std::ostream& operator<<(std::ostream &out, const AlgebraicMatrix<DA> &obj)
{
/*! Specialized stream output operator for DA matrices.
   \param[in] out Standard output stream.
   \param[in] obj AlgebraicMatrix<DA> to be printed in the stream
   \return Reference to the standard output stream.
 */
    const unsigned int nrows = obj.nrows();
    const unsigned int ncols = obj.ncols();

    out << "[[[ " << nrows << "x" << ncols << " matrix" << std::endl;
    for(unsigned int j = 0; j<ncols; j++) {
        out << "    Column " << j+1 << std::endl;
        for(unsigned int i = 0; i<nrows; i++) {
            out << obj.at(i,j); }
    }
    out << "]]]" << std::endl;

    return out;
}

template<> std::istream& operator>>(std::istream &in, AlgebraicMatrix<DA> &obj){
/*! Specialized stream input operator for DA matrices.
   \param[in] in Standard input stream.
   \param[out] obj AlgebraicMatrix<DA> to be read from the stream
   \return Reference to the standard input stream.
*/
    // read the first line
    std::string line;
    std::getline(in, line);

    unsigned int n_rows = 0;
    unsigned int n_cols = 0;

    if(in.good()){
        // Find the number of rows
        std::size_t found = line.find_first_of('x');
        std::string size_str;

        for( unsigned int j = 4; j<found; j++ )
            size_str += line[j];

        if(!(std::istringstream(size_str) >> n_rows))
            n_rows = 0;

        // Find the number of columns (stop when 'm' of "matrix" is met)
        std::size_t found2 = line.find_first_of('m', found);
        size_str.clear();

        for( unsigned int j = found+1; j<found2; j++)
            size_str += line[j];

        if(!(std::istringstream(size_str) >> n_cols))
            n_cols = 0;

        // resize the object to meet the size of the vector to be read
        obj.resize(n_rows, n_cols);

        // fill the AlgebraicMatrix
        for( unsigned int j = 0; j<n_cols; j++){
            // skip the line at the beginning of the column (containing Column X)
            std::getline(in, line);

            // Read all columns in current row j
            for(unsigned int i=0; i<n_cols; i++)
                in >> obj.at(i,j);
        }

        // skip the line at the end of a AlgebraicMatrix (containing ]]])
        std::getline(in, line);
    } else {
        obj.resize(0, 0); }

    return in;
}

}
