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
 * AlgebraicVector.cpp
 *
 *  Created on: January 23, 2015
 *      Author: Dinamica Srl
 */

// C++ stdlib classes used only internally in this implementation
#include <algorithm>

// DACE classes
#include "dace/config.h"
#include "dace/DA.h"
#include "dace/compiledDA.h"
#include "dace/compiledDA_t.h"
#include "dace/AlgebraicVector.h"
#include "dace/AlgebraicVector_t.h"

namespace DACE {

/***********************************************************************************
*     Coefficient access routines
************************************************************************************/
#ifdef WITH_ALGEBRAICMATRIX
template<> AlgebraicMatrix<double> AlgebraicVector<DA>::linear() const {
/*! Return the linear part of a AlgebraicVector<T>. NOT DEFINED FOR TYPES OTHER THAN DA.
   \return A AlgebraicMatrix<double> of dimension size by nvar, where size is the
    size of the AlgebraicVector<T> considered and nvar is the number of variables defined
    during the DACE initialization. Each row contains the linear part of the corresponding
    DA included in the original AlgebraicVector<T>.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
 */
    const size_t size = this->size();
    const int nvar = DA::getMaxVariables();

    AlgebraicMatrix<double> out(size, nvar);
    for(size_t i=0; i<size; i++) {
          out.setrow( i, (*this)[i].linear() );
    }
    return out;
}
#else
template<> std::vector< std::vector<double> > AlgebraicVector<DA>::linear() const {
/*! Return the linear part of a AlgebraicVector<T>. NOT DEFINED FOR TYPES OTHER THAN DA.
   \return A std::vector< std::vector<double> >, where each std::vector<double> contains
    the linear part of the corresponding DA included in the original AlgebraicVector<T>.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
 */
    const size_t size = this->size();

    std::vector< std::vector<double> > out(size);
    for(size_t i=0; i<size; i++) {
          out[i] = (*this)[i].linear();
    }
    return out;
}
#endif /* WITH_ALGEBRAICMATRIX */

template<> AlgebraicVector<DA> AlgebraicVector<DA>::trim(const unsigned int min, const unsigned int max) const {
/*! Returns an AlgebraicVector<DA> with all monomials of order less than min and greater than max removed (trimmed). The result is copied in a new AlgebraicVector<DA>.
   \param[in] min minimum order to be preserved.
   \param[in] max maximum order to be preserved.
   \return A new AlgebraicVector<DA> containing the result of the trimming.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
*/
    AlgebraicVector<DA> tmp(this->size());

    for(size_t i=0; i<this->size(); i++)
        tmp[i] = (*this)[i].trim(min, max);

    return tmp;
}

/***********************************************************************************
*     Math routines
************************************************************************************/
template<> AlgebraicVector<DA> AlgebraicVector<DA>::deriv(const unsigned int p) const {
/*! Compute the derivative of a AlgebraicVector<T> with respect to variable p.
    The result is copied in a new AlgebraicVector<T>. NOT DEFINED FOR TYPES OTHER THAN DA.
   \param[in] p variable with respect to which the derivative is calculated.
   \return A new AlgebraicVector<T> containing the result of the derivation.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
 */
    const size_t size = this->size();
    AlgebraicVector<DA> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = (*this)[i].deriv(p);}

    return temp;
}

template<> AlgebraicVector<DA> AlgebraicVector<DA>::integ(const unsigned int p) const {
/*! Compute the integral of a AlgebraicVector<T> with respect to variable p.
    The result is copied in a new AlgebraicVector<T>. NOT DEFINED FOR TYPES OTHER THAN DA.
   \param[in] p variable with respect to which the integral is calculated.
   \return A new AlgebraicVector<T> containing the result of the integration.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
 */
    const size_t size = this->size();
    AlgebraicVector<DA> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = (*this)[i].integ(p);}

    return temp;
}

/***********************************************************************************
*     Polynomial evaluation routines
************************************************************************************/
template<> compiledDA AlgebraicVector<DA>::compile() const {
/*! Compile vector of polynomials and create a compiledDA object.
   \return The compiled DA object.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
 */
    return compiledDA(*this);
}

template<> AlgebraicVector<DA> AlgebraicVector<DA>::plug(const unsigned int var, const double val) const {
/*! Partial evaluation of vector of polynomials. In each element of the vector,
    variable var is replaced by the value val. The resulting vector of DAs
    is returned.
   \param[in] var variable number to be replaced
   \param[in] val value by which to replace the variable
   \return A new DA object containing the resulting DA object.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
 */
    const size_t size = this->size();
    AlgebraicVector<DA> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = (*this)[i].plug(var,val);}

    return temp;
}

#ifndef WITH_ALGEBRAICMATRIX
/*! \cond */
/* Internal routine to compute a matrix inverse of a double precision matrix.
   Algorithm based on the Gauss elimination with full pivot (from the Numerical
   Cookbook) adapted for C++. This is the same algorithm but a different
   implementation than in the AlgebraicMatrix class.
   This is NOT intended for public use. Limited error checking is performed
   in accordance with the exclusive use of this routine in map inversion.
 */
template<> void AlgebraicVector<DA>::matrix_inverse(std::vector< std::vector<double> > &A) {
    using std::abs;

    const size_t n = A.size();
    std::vector<size_t> indexc(n), indexr(n), ipiv(n, 0);

    for (size_t i=0; i<n; i++) {
        size_t icol = 0, irow = 0;
        double big = 0.0;
        for (size_t j=0; j<n; j++)
            if (ipiv[j] == 0)
                for (size_t k=0; k<n; k++)
                    if (ipiv[k] == 0)
                        if (abs(A[j][k]) >= big) {
                            big = abs(A[j][k]);
                            irow = j;
                            icol = k;}
        ipiv[icol] = 1;
        if (irow != icol)
            for (size_t l=0; l<n; l++) std::swap(A[irow][l], A[icol][l]);
        indexr[i] = irow;
        indexc[i] = icol;
        if (A[icol][icol] == 0.0) throw std::runtime_error("DACE::AlgebraicVector<DA>::inverse: linear matrix inverse does not exist.");
        const double pivinv = 1.0/A[icol][icol];
        A[icol][icol] = 1.0;
        for (size_t l=0; l<n; l++) A[icol][l] *= pivinv;
        for (size_t ll=0; ll<n; ll++)
            if (ll != icol) {
                const double temp = A[ll][icol];
                A[ll][icol] = 0.0;
                for (size_t l=0; l<n; l++) A[ll][l] -= A[icol][l]*temp;}
    }

    for (size_t i=n; i>0; i--)
        if (indexr[i-1] != indexc[i-1])
            for (size_t k=0; k<n; k++)
                std::swap(A[k][indexr[i-1]], A[k][indexc[i-1]]);
}
/*! \endcond */
#endif /* WITH_ALGEBRAICMATRIX */

template<> AlgebraicVector<DA> AlgebraicVector<DA>::invert() const {
/*! Invert the polynomials map given by the AlgebraicVector<DA>.
   \return the inverted polynomials
   \throw std::runtime_error
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
*/
    const unsigned int ord = DA::getTO();
    const size_t nvar = this->size();

    if(nvar>DA::getMaxVariables())
        throw std::runtime_error("DACE::AlgebraicVector<DA>::inverse: dimension of vector exceeds maximum number of DA variables.");

    // Create DA identity
    AlgebraicVector<DA> DDA = AlgebraicVector<DA>::identity(nvar);

    // Split map into constant part AC, non-constant part M, and non-linear part AN
    AlgebraicVector<double> AC = this->cons();
    AlgebraicVector<DA> M = this->trim(1);
    AlgebraicVector<DA> AN = M.trim(2);

#ifdef WITH_ALGEBRAICMATRIX
    // Extract the linear coefficients matrix
    AlgebraicMatrix<double> AL = M.linear();

    // Compute the inverse of linear coefficients matrix
    AlgebraicMatrix<double> AI = AL.inv();

    // Compute DA representation of the inverse of the linear part of the map and its composition with non-linear part AN
    compiledDA AIoAN(AI*AN);
    AlgebraicVector<DA> Linv = AI*DDA;
#else
    // Compute the inverse of linear coefficients matrix
    std::vector< std::vector<double> > AI = M.linear();
    matrix_inverse(AI);

    // Compute DA representation of the inverse of the linear part of the map and its composition with non-linear part AN
    AlgebraicVector<DA> Linv(nvar);
    // Linv = AI*AN
    for(size_t i=0; i<nvar; i++) {
        Linv[i] = 0.0;
        for(size_t j=0; j<nvar; j++)
            Linv[i] += AI[i][j]*AN[j];}
    compiledDA AIoAN(Linv);
    // Linv = AI*DDA
    for(size_t i=0; i<nvar; i++) {
        Linv[i] = 0.0;
        for(size_t j=0; j<nvar; j++)
            Linv[i] += AI[i][j]*DDA[j];}
#endif /* WITH_ALGEBRAICMATRIX */

    // Iterate to obtain the inverse map
    AlgebraicVector<DA> MI = Linv;
    for(unsigned int i=1; i<ord; i++) {
        DA::setTO(i+1);
        MI = Linv - AIoAN.eval(MI);}

    return MI.eval(DDA-AC);
}

/********************************************************************************
*     Static factory routines
*********************************************************************************/
template<> AlgebraicVector<DA> AlgebraicVector<DA>::identity(const size_t n) {
/*! Return the DA identity of dimension n.
   \param[in] n The dimendion of the identity.
   \return AlgebraicVector<DA> containing the DA identity in n dimensions
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
 */
    AlgebraicVector<DA> temp(n);
    for(size_t i=0; i < n; i++) {
        temp[i] = DA((int)(i+1));}

    return temp;
}

/***********************************************************************************
*     Non-member functions
************************************************************************************/
#ifdef WITH_ALGEBRAICMATRIX
template<> AlgebraicMatrix<double> linear(const AlgebraicVector<DA> &obj) {
/*! Return the linear part of a AlgebraicVector<T>.
   \param[in] obj AlgebraicVector<T> to extract linear part from
   \return An AlgebraicMatrix<double> of dimensions size by nvar, where
    size is the size of the AlgebraicVector<T> considered and nvar is the
    number of variables defined during the DACE initialization.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
   \sa AlgebraicVector<T>::linear
 */
#else
template<> std::vector< std::vector<double> > linear(const AlgebraicVector<DA> &obj) {
/*! Return the linear part of a AlgebraicVector<T>. Only defined for AlgebraicVector<DA>.
   \param[in] obj AlgebraicVector<T> to extract linear part from
   \return A std::vector< std::vector<double> > containing the linear parts of
    each component of the AlgebraicVector<DA> obj.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
   \sa AlgebraicVector<T>::linear
 */
#endif /* WITH_ALGEBRAICMATRIX */
    return obj.linear();
}

template<> AlgebraicVector<DA> trim(const AlgebraicVector<DA> &obj, unsigned int min, unsigned int max) {
/*! Returns an AlgebraicVector<DA> with all monomials of order less than min and greater than max removed (trimmed). The result is copied in a new AlgebraicVector<DA>.
   \param[in] obj the AlgebraicVector<DA> to be trimmed.
   \param[in] min minimum order to be preserved.
   \param[in] max maximum order to be preserved.
   \return A new AlgebraicVector<DA> containing the result of the trimming.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
   \sa AlgebraicVector<T>::trim
*/
    return obj.trim(min, max);
}

template<> AlgebraicVector<DA> deriv(const AlgebraicVector<DA> &obj, const unsigned int p) {
/*! Compute the derivative of a AlgebraicVector<T> with respect to variable p.
    The result is copied in a new AlgebraicVector<T>.
   \param[in] obj AlgebraicVector<T>.
   \param[in] p variable with respect to which the derivative is calculated.
   \return A new AlgebraicVector<T> containing the result of the derivation.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
   \sa AlgebraicVector<T>::deriv
 */
    return obj.deriv(p);
}

template<> AlgebraicVector<DA> integ(const AlgebraicVector<DA> &obj, const unsigned int p) {
/*! Compute the integral of a AlgebraicVector<T> with respect to variable p.
    The result is copied in a new AlgebraicVector<T>.
   \param[in] obj AlgebraicVector<T>.
   \param[in] p variable with respect to which the integral is calculated.
   \return A new AlgebraicVector<T> containing the result of the integration.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
   \sa AlgebraicVector<T>::integ
 */
    return obj.integ(p);
}

template<> compiledDA compile(const AlgebraicVector<DA> &obj) {
/*! Compile vector of polynomials and create a compiledDA object.
   \param[in] obj The AlgebraicVector to compile.
   \return The compiled DA object.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
 */
    return obj.compile();
}

template<> AlgebraicVector<DA> plug(const AlgebraicVector<DA> &obj, const unsigned int var, const double val) {
/*! Partial evaluation of vector of polynomials. In each element of the vector,
    variable var is replaced by the value val. The resulting vector of DAs
    is returned.
   \param[in] obj The vector to partially evaluate.
   \param[in] var Variable number to be replaced.
   \param[in] val Value by which to replace the variable.
   \return A new DA object containing the resulting DA object.
   \note This DA specific function is only available in AlgebraicVector<DA>.
    When called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
 */
    return obj.plug(var, val);
}

}
