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
 * AlgebraicVector_t.h
 *
 *  Created on: Sep. 10, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_ALGEBRAICVECTOR_T_H_
#define DINAMICA_ALGEBRAICVECTOR_T_H_

// C++ stdlib classes used only internally in this implementation
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

// DACE classes
#include "dace/PromotionTrait.h"
#include "dace/MathExtension.h"
#include "dace/compiledDA.h"
#include "dace/AlgebraicVector.h"
#ifdef WITH_ALGEBRAICMATRIX
#include "dace/AlgebraicMatrix.h"
#include "dace/AlgebraicMatrix_t.h"
#endif /* WITH_ALGEBRAICMATRIX */

namespace DACE {

/***********************************************************************************
*     Constructors
************************************************************************************/
template<typename T> AlgebraicVector<T>::AlgebraicVector() : std::vector<T>() {
/*! Default Constructor to create empty AlgebraicVector
 */
}

template<typename T> AlgebraicVector<T>::AlgebraicVector(size_t size) : std::vector<T>(size) {
/*! Constructor with size to allocate a vector of the given size with elements initialized using their default constructor.
   \param[in] size length of AlgebraicVector.
 */
}

template<typename T> AlgebraicVector<T>::AlgebraicVector(size_t size, const T &d) : std::vector<T>(size, d) {
/*! Constructor with size and elements value to allocate a vector of the given size with elements initialized as copies of d.
   \param[in] size length of AlgebraicVector
   \param[in] d    initial value for the elements
 */
}

template<typename T> AlgebraicVector<T>::AlgebraicVector(const std::vector<T> &v) : std::vector<T>(v) {
/*! Copy constructor to create a copy of any existing vector.
   \param[in] v vector to be copied into AlgebraicVector
 */
}

template<typename T> AlgebraicVector<T>::AlgebraicVector(std::initializer_list<T> l) : std::vector<T>(l) {
/*! Constructor to create a vector from an initializer list.
   \param[in] l braced initializer list to be copied into the AlgebraicVector
 */
}

template<typename T> AlgebraicVector<T>::AlgebraicVector(const std::vector<T> &v, size_t first, size_t last) : std::vector<T>(v.begin()+first, v.begin()+last+1) {
/*! Extraction constructor to copy only a given range of elements from vector v.
   \param[in] v vector to be copied into AlgebraicVector
   \param[in] first index of the first element to be copied
   \param[in] last index of the last element to be copied
   \note The constructor does not perform any range checking for the extraction.
   \sa AlgebraicVector<T>::extract
 */
    // Notice that the range in the std::vector constructor above includes all
    // elements between first and last, including the first excluding the last.
    // Hence the +1.
}

/***********************************************************************************
*     Coefficient access routines
************************************************************************************/
template<typename T> AlgebraicVector<double> AlgebraicVector<T>::cons() const {
/*! Return the constant parts of each element.
   \return A AlgebraicVector<double> containing the constant part of each element.
 */
    using DACE::cons;

    const size_t size = this->size();
    AlgebraicVector<double> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = cons((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::extract(const size_t first, const size_t last) const {
/*! Extracts elements from AlgebraicVector.
   \param[in] first index of first element to be extracted
   \param[in] last  index of last element to be extracted
   \return A new AlgebraicVector<T> with elements from position first to last.
   \throw std::runtime_error
   \note This routine performs range checking and throws an error if the indices are out of range.
*/
    if(first>=this->size() || last>=this->size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::take: Indices out of bounds.");

    return AlgebraicVector<T>(*this, first, last);
}

template<typename T> template<typename V> AlgebraicVector<typename PromotionTrait< T, V >::returnType> AlgebraicVector<T>::concat(const std::vector<V> &obj) const {
/*! Append an AlgebraicVector to the end of the current one and return the new vector.
   \param[in] obj The AlgebraicVector to be appended.
   \return A new AlgebraicVector containing the elements of both vectors,
    cast upwards if necessary.
*/
    const size_t size1 = this->size();
    const size_t size2 = obj.size();
    AlgebraicVector<typename PromotionTrait< T, V >::returnType> res(size1+size2);

    for(size_t i=0; i<size1; i++)
        res[i] = (*this)[i];
    for(size_t i=0; i<size2; i++)
        res[i+size1] = obj[i];

    return res;
}

/***********************************************************************************
*     Algebraic operations
************************************************************************************/
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::operator-() const {
/*! Returns the additive inverse of the vector.
   \return A new AlgebraicVector, with the opposite sign.
 */
    return -1.0*(*this);
}

template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator+=(const AlgebraicVector<U> &obj) {
/*! Add the given AlgebraicVector to ourselves.
   \param[in] obj An AlgebraicVector.
   \return A reference to ourselves.
   \throw std::runtime_error
 */
    const size_t size = this->size();
    if(size != obj.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator+=: Vectors must have the same length.");

    for(size_t i=0; i<size; i++) {
        (*this)[i] += obj[i];}
    return *this;
}

template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator+=(const U &obj) {
/*! Add the given scalar to ourselves componentwise.
   \param[in] obj A scalar value.
   \return A reference to ourselves.
 */
    const size_t size = this->size();
    for(size_t i=0; i<size; i++) {
        (*this)[i] += obj;}
    return *this;
}

template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator-=(const AlgebraicVector<U> &obj) {
/*! Subtract the given AlgebraicVector from ourselves.
   \param[in] obj An AlgebraicVector.
   \return A reference to ourselves.
   \throw std::runtime_error
 */
    const size_t size = this->size();
    if(size != obj.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator-=: Vectors must have the same length.");

    for(size_t i=0; i<size; i++) {
        (*this)[i] -= obj[i];}
    return *this;
}

template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator-=(const U &obj) {
/*! Subtract the given scalar from ourselves componentwise.
   \param[in] obj A scalar value.
   \return A reference to ourselves.
 */
    const size_t size = this->size();
    for(size_t i=0; i<size; i++) {
        (*this)[i] -= obj;}
    return *this;
}

template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator*=(const AlgebraicVector<U> &obj) {
/*! Multiply the given AlgebraicVector with ourselves componentwise.
   \param[in] obj An AlgebraicVector.
   \return A reference to ourselves.
 */
    const size_t size = this->size();
    if(size != obj.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator*=: Vectors must have the same length.");

    for(size_t i=0; i<size; i++) {
        (*this)[i] *= obj[i];}
    return *this;
}

template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator*=(const U &obj) {
/*! Multiply the given scalar with ourselves.
   \param[in] obj A scalar value.
   \return A reference to ourselves.
 */
    const size_t size = this->size();
    for(size_t i=0; i<size; i++) {
        (*this)[i] *= obj;}
    return *this;
}

template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator/=(const AlgebraicVector<U> &obj) {
/*! Divide ourselves by the given AlgebraicVector componentwise.
   \param[in] obj An AlgebraicVector.
   \return A reference to ourselves.
   \throw std::runtime_error
 */
    const size_t size = this->size();
    if(size != obj.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator/=: Vectors must have the same length.");

    for(size_t i=0; i<size; i++) {
        (*this)[i] /= obj[i];}
    return *this;
}

template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator/=(const U &obj) {
/*! Divide ourselves by the given scalar.
   \param[in] obj A scalar value.
   \return A reference to ourselves.
 */
    const size_t size = this->size();
    for(size_t i=0; i<size; i++) {
        (*this)[i] /= obj;}
    return *this;
}

template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator<<(const std::vector<U> &obj) {
/*! Append elements of vector obj to the end of ourself, converting the type to match ours if necessary.
   \param[in] obj Vector of elements to append.
   \return A reference to ourselves.
 */
    const size_t size = obj.size();
    for(size_t i=0; i<size; i++) {
        (*this).push_back((T)obj[i]);}
    return *this;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator+(const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2) {
/*! Compute the addition between two AlgebraicVectors.
   \param[in] obj1 first AlgebraicVector.
   \param[in] obj2 second AlgebraicVector.
   \return A new AlgebraicVector, (obj1+obj2).
   \throw std::runtime_error
 */
    if(obj1.size() != obj2.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator+: Vectors must have the same length.");

    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] + obj2[i];}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator+(const AlgebraicVector<U> &obj1, const V &obj2) {
/*! Compute the addition between a AlgebraicVector and a scalar value.
   \param[in] obj1 a AlgebraicVector.
   \param[in] obj2 a scalar value.
   \return A new AlgebraicVector, (obj1+obj2).
 */
    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] + obj2;}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator+(const U &obj1, const AlgebraicVector<V> &obj2) {
/*! Compute the addition between a scalar value and a AlgebraicVector.
   \param[in] obj1 a scalar value.
   \param[in] obj2 a AlgebraicVector.
   \return A new AlgebraicVector, (obj1+obj2).
 */
    const size_t size = obj2.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1 + obj2[i];}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator-(const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2) {
/*! Compute the subtraction between two AlgebraicVectors.
   \param[in] obj1 first AlgebraicVector.
   \param[in] obj2 second AlgebraicVector.
   \return A new AlgebraicVector, (obj1-obj2).
   \throw std::runtime_error
 */
    if(obj1.size() != obj2.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator-: Vectors must have the same length.");

    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] - obj2[i];}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator-(const AlgebraicVector<U> &obj1, const V &obj2) {
/*! Compute the subtraction between a AlgebraicVector and a scalar value.
   \param[in] obj1 a AlgebraicVector.
   \param[in] obj2 a scalar value.
   \return A new AlgebraicVector, (obj1-obj2).
 */
    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] - obj2;}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator-(const U &obj1, const AlgebraicVector<V> &obj2) {
/*! Compute the subtraction between a scalar value and a AlgebraicVector.
   \param[in] obj1 a scalar value.
   \param[in] obj2 a AlgebraicVector.
   \return A new AlgebraicVector, (obj1-obj2).
 */
    const size_t size = obj2.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1 - obj2[i];}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator*(const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2) {
/*! Compute the element-wise multiplication between two AlgebraicVectors.
   \param[in] obj1 first AlgebraicVector.
   \param[in] obj2 second AlgebraicVector.
   \return A new AlgebraicVector, (obj1*obj2).
   \throw std::runtime_error
 */
    if(obj1.size() != obj2.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator*: Vectors must have the same length.");

    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] * obj2[i];}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator*(const AlgebraicVector<U> &obj1, const V &obj2) {
/*! Compute the multiplication between a AlgebraicVector and a scalar value.
   \param[in] obj1 a AlgebraicVector.
   \param[in] obj2 a scalar value.
   \return A new AlgebraicVector, (obj1*obj2).
 */
    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] * obj2;}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator*(const U &obj1, const AlgebraicVector<V> &obj2) {
/*! Compute the multiplication between a scalar value and a AlgebraicVector.
   \param[in] obj1 a scalar value.
   \param[in] obj2 a AlgebraicVector.
   \return A new AlgebraicVector, (obj1*obj2).
 */
    const size_t size = obj2.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1 * obj2[i];}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator/(const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2) {
/*! Compute the element-wise division between two AlgebraicVectors.
   \param[in] obj1 first AlgebraicVector.
   \param[in] obj2 second AlgebraicVector.
   \return A new AlgebraicVector, (obj1/obj2).
   \throw std::runtime_error
 */
    if(obj1.size() != obj2.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator/: Vectors must have the same length.");

    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] / obj2[i];}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator/(const AlgebraicVector<U> &obj1, const V &obj2) {
/*! Compute the division between a AlgebraicVector and a scalar value.
   \param[in] obj1 a AlgebraicVector.
   \param[in] obj2 a scalar value.
   \return A new AlgebraicVector, (obj1/obj2).
 */
    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] / obj2;}
    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator/(const U &obj1, const AlgebraicVector<V> &obj2) {
/*! Compute the division between a scalar value and a AlgebraicVector.
   \param[in] obj1 a scalar value.
   \param[in] obj2 a AlgebraicVector.
   \return A new AlgebraicVector, (obj1/obj2).
 */
    const size_t size = obj2.size();
    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1 / obj2[i];}
    return temp;
}

/***********************************************************************************
*     Math routines
************************************************************************************/
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::absolute() const {
/*! Element-wise application of the absolute value function.
   \return A new AlgebraicVector<T>.
 */
    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = absolute((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::trunc() const {
/*! Element-wise application of the truncation function.
   \return A new AlgebraicVector<T>.
 */
    using std::trunc;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = trunc((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::round() const {
/*! Element-wise application of the round function.
   \return A new AlgebraicVector<T>.
 */
    using std::round;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = round((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::mod() const {
/*! Element-wise application of the mod function.
   \return A new AlgebraicVector<T>.
 */
    using std::trunc;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = mod((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::pow(const int p) const {
/*! Element-wise application of the integer power function.
   \param[in] p power to raise each element to.
   \return A new AlgebraicVector<T>.
 */
    using std::pow;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = pow((*this)[i], p);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::pow(const double p) const {
/*! Element-wise application of the double power function.
   \param[in] p power to raise each element to.
   \return A new AlgebraicVector<T>.
 */
    using std::pow;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = pow((*this)[i], p);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::root(const int p) const {
/*! Element-wise application of the p-th root function.
   \param[in] p root to be computed.
   \return A new AlgebraicVector<T>.
 */
    using DACE::root;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = root((*this)[i],p);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::minv() const {
/*! Element-wise application of the multiplicative inverse function.
   \return A new AlgebraicVector<T>.
 */
    using DACE::minv;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = minv((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::sqr() const {
/*! Element-wise application of the square function.
   \return A new AlgebraicVector<T>.
 */
    using DACE::sqr;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = sqr((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::sqrt() const {
/*! Element-wise application of the square root function.
   \return A new AlgebraicVector<T>.
 */
    using std::sqrt;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = sqrt((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::isrt() const {
/*! Element-wise application of the inverse square root function.
   \return A new AlgebraicVector<T>.
 */
    using DACE::isrt;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = isrt((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::cbrt() const {
/*! Element-wise application of the cube root function.
   \return A new AlgebraicVector<T>.
 */
    using std::cbrt;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = cbrt((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::icbrt() const {
/*! Element-wise application of the inverse cube root function.
   \return A new AlgebraicVector<T>.
 */
    using DACE::icbrt;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = sqrt((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::hypot(const AlgebraicVector<T> &obj) const {
/*! Element-wise application of the hypotenuse function hypot(x,y).
   \param[in] obj AlgebraicVector<T> second argument.
   \return A new AlgebraicVector<T>.
*/
    using std::hypot;

    const size_t size = this->size();
    if(obj.size() != size)
        throw std::runtime_error("DACE::AlgebraicVector<T>::hypot(): Vectors must have the same length.");

    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = hypot((*this)[i], obj[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::exp() const {
/*! Element-wise application of the exponential function.
   \return A new AlgebraicVector<T>.
 */
    using std::exp;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = exp((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::log() const {
/*! Element-wise application of the natural logarithm function.
   \return A new AlgebraicVector<T>.
 */
    using std::log;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = log((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::logb(const double b) const {
/*! Element-wise application of the logarithm function relative to given base.
   \param[in] b base for the logarithm
   \return A new AlgebraicVector<T>.
 */
    using DACE::logb;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = logb((*this)[i], b);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::log10() const {
/*! Element-wise application of the decadic logarithm function.
   \return A new AlgebraicVector<T>.
 */
    using std::log10;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = log10((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::log2() const {
/*! Element-wise application of the binary logarithm function.
   \return A new AlgebraicVector<T>.
 */
    using std::log2;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = log2((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::sin() const {
/*! Element-wise application of the sine function.
   \return A new AlgebraicVector<T>.
 */
    using std::sin;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = sin((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::cos() const {
/*! Element-wise application of the cosine function.
   \return A new AlgebraicVector<T>.
 */
    using std::cos;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = cos((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::tan() const {
/*! Element-wise application of the tangent function.
   \return A new AlgebraicVector<T>.
 */
    using std::tan;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = tan((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::asin() const {
/*! Element-wise application of the arcsine function.
   \return A new AlgebraicVector<T>.
 */
    using std::asin;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = asin((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::acos() const {
/*! Element-wise application of the arccosine function.
   \return A new AlgebraicVector<T>.
 */
    using std::acos;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = acos((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::atan() const {
/*! Element-wise application of the arctangent function.
   \return A new AlgebraicVector<T>.
 */
    using std::atan;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = atan((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::atan2(const AlgebraicVector<T> &obj) const {
/*! Element-wise application of the four-quadrant arctangent of Y/X.
    Y is the current object, whereas X is the argument.
   \param[in] obj AlgebraicVector<T> representing X
   \return A new AlgebraicVector<T> with elements in [-pi, pi].
   \throw std::runtime_error
*/
    using std::atan2;

    const size_t size = this->size();
    if(obj.size() != size)
        throw std::runtime_error("DACE::AlgebraicVector<T>::atan2(): Vectors must have the same length.");

    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = atan2((*this)[i], obj[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::sinh() const {
/*! Element-wise application of the hyperbolic sine function.
   \return A new AlgebraicVector<T>.
 */
    using std::sinh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = sinh((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::cosh() const {
/*! Element-wise application of the hyperbolic cosine function.
   \return A new AlgebraicVector<T>.
 */
    using std::cosh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = cosh((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::tanh() const {
/*! Element-wise application of the hyperbolic tangent function.
   \return A new AlgebraicVector<T>.
 */
    using std::tanh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = tanh((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::asinh() const {
/*! Element-wise application of the hyperbolic arcsine function.
   \return A new AlgebraicVector<T>.
 */
    using std::asinh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = asinh((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::acosh() const {
/*! Element-wise application of the hyperbolic arccosine function.
   \return A new AlgebraicVector<T>.
 */
    using std::acosh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = acosh((*this)[i]);}

    return temp;
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::atanh() const {
/*! Element-wise application of the hyperbolic arctangent function.
   \return A new AlgebraicVector<T>.
 */
    using std::atanh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = atanh((*this)[i]);}

    return temp;
}

/***********************************************************************************
*    Vector routines
************************************************************************************/
template<typename T> template<typename V> typename PromotionTrait<T,V>::returnType AlgebraicVector<T>::dot(const AlgebraicVector<V> &obj) const {
/*! Compute the dot product with another AlgebraicVector.
   \param[in] obj the other AlgebraicVector.
   \return A scalar value,.
   \throw std::runtime_error
 */
    const size_t size = this->size();
    if(size != obj.size())
          throw std::runtime_error("DACE::AlgebraicVector<T>::dot(): Vectors must have the same length.");

    typename PromotionTrait<T,V>::returnType temp = 0.0;
    for(size_t i=0; i<size; i++) {
        temp += (*this)[i] * obj[i];}

    return temp;
}

template<typename T> template<typename V> AlgebraicVector<typename PromotionTrait<T,V>::returnType> AlgebraicVector<T>::cross(const AlgebraicVector<V> &obj) const {
/*! Compute the cross product with another 3D AlgebraicVector.
   \param[in] obj The other AlgebraicVector.
   \return A new AlgebraicVector.
   \throw std::runtime_error
 */
    if((this->size() != 3) || (obj.size() != 3))
        throw std::runtime_error("DACE::AlgebraicVector<T>::cross(): Inputs must be 3 element AlgebraicVectors.");

    AlgebraicVector<typename PromotionTrait<T,V>::returnType> temp(3);

    temp[0] = ((*this)[1] * obj[2]) - ((*this)[2] * obj[1]);
    temp[1] = ((*this)[2] * obj[0]) - ((*this)[0] * obj[2]);
    temp[2] = ((*this)[0] * obj[1]) - ((*this)[1] * obj[0]);

    return temp;
}

template<typename T> T AlgebraicVector<T>::length() const {
/*! Compute the length (Euclidean vector norm).
   \return Length of the vector.
 */
    using std::sqrt; using DACE::sqr;      // Implementational note: these using statements are very subtle and absolutely needed.
                                            // They force the compiler to perform argument dependent lookup (ADL) which then finds
                                            // the correct sqrt() and sqr() functions even if they are not in DACE:: or std::!
    const size_t size = this->size();
    T norm = 0.0;
    for(size_t i=0; i<size; i++) {
        norm = norm + sqr((*this)[i]);}

    return sqrt(norm);
}

template<typename T> AlgebraicVector<T> AlgebraicVector<T>::normalize() const {
/*! Normalize the vector.
   \return An AlgebraicVector<T> of unit length.
 */
    using DACE::minv;

    AlgebraicVector<T> temp(*this);
    temp *= minv(this->length());
    return temp;
}

/***********************************************************************************
*     Polynomial evaluation routines
************************************************************************************/
template<> template<typename V> V AlgebraicVector<DA>::eval(const V &args) const {
/*! Evaluate a vector of polynomials with any vector type V with arguments
    and return a vector of results of the same type V.
   \param[in] args vector (e.g. AlgebraicVector<>) of arguments
   \return A new vector of same type as argument args containing the
    results of the evaluation.
   \note This DA specific function is only available in AlgebraicVector<DA>,
    when called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
 */

    return compiledDA(*this).eval(args);
}

template<> template<typename U> AlgebraicVector<U> AlgebraicVector<DA>::eval(const std::initializer_list<U> l) const {
/*! Evaluate a vector of polynomials with an braced initializer list of type U
    and return an AlgebraicVector of type U with the results.
   \param[in] l Braced initializer list containing the arguments.
   \return A new AlgebraicVector of type U containing the results of the evaluation.
   \note C++ is not able to derive the type of elements of an initializer list automatically.
    That means eval() must be called explicitly as e.g. eval<double>({1.0, 2.0, 3.0}) when
    used with initializer lists.
 */
    return compiledDA(*this).eval<U>(l);
}

template<> template<typename U> AlgebraicVector<U> AlgebraicVector<DA>::evalScalar(const U &arg) const {
/*! Evaluate a vector of polynomials with a single arithmetic type U argument.
   \param[in] arg single variable of arithmetic type T of the first independent DA variable.
   \return The result of the evaluation.
   \note This DA specific function is only available in AlgebraicVector<DA>,
    when called on AlgebraicVectors of other types (e.g. double), a compiler
    error will be the result.
   \note To be used only for single polynomial evaluation. For multiple
    evaluations of the same vector of polynomials use the corresponding method
    in class compiledDA.
   \sa compiledDA
   \sa AlgebraicVector::compile()
 */
    return compiledDA(*this).evalScalar(arg);
}

/***********************************************************************************
*     DA norm routines
************************************************************************************/
template<typename T> AlgebraicVector<double> AlgebraicVector<T>::norm(const unsigned int type) const {
/*! Element-wise application of the norm function.
   \param[in] type type of norm to be computed. See DA::norm.
   \return A new AlgebraicVector<T>.
   \sa DA::norm()
 */
    using DACE::norm;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = norm((*this)[i], type);}

    return temp;
}

/***********************************************************************************
*     Input/Output routines
************************************************************************************/
template<typename U> std::ostream& operator<<(std::ostream &out, const AlgebraicVector<U> &obj) {
/*! Output a vector to a C++ output stream.
   \param[in] out standard output stream.
   \param[in] obj AlgebraicVector<U> to be written to the stream
   \return Reference to output stream out.
 */
    const size_t size = obj.size();

    out << "[[[ " << size << " vector" << std::endl;
    for(size_t i=0; i<size;i++) {
        out << obj[i] << std::endl;}
    out << "]]]" << std::endl;

    return out;
}

template<typename U> std::istream& operator>>(std::istream &in, AlgebraicVector<U> &obj) {
/*! Read a vector from a C++ input stream.
   \param[in] in standard input stream.
   \param[in] obj AlgebraicVector<U> to be read from the stream
   \return Reference to input stream in.
 */
    std::string init_line;
    size_t vec_size;

    // try to read the first line
    getline(in, init_line);
    if(in.good()) {
        // retrieve the size of the vector to be read
        std::size_t found = init_line.find_first_of(' ');
        std::string size_str(init_line,4,found-4);
        if(!(std::istringstream(size_str) >> vec_size)) vec_size = 0;

        // resize the object to meet the size of the vector to be read
        obj.resize(vec_size);

        // fill the AlgebraicVector
        for (size_t i = 0; in.good() && (i < vec_size); i++) {
            in >> obj[i];

            // check the next character
            if (in.peek() == '\n')       // the previous operator>> does not consume the \n character when an AlgebraicVector<T> (with T != DA) is considered
                in.ignore();            // ignore the next character
        }
        // skip the line at the end of a AlgebraicVector (containing ]]])
        getline(in, init_line);
    }else{
        obj.resize(0);}

    return in;
}

template<typename T> std::string AlgebraicVector<T>::toString() const {
/*! Convert to string.
    \return A string.
 */
    std::ostringstream strs;
    strs << *this << std::endl;

    return strs.str();
}

/***********************************************************************************
*     Non-member functions
************************************************************************************/
template<typename T> AlgebraicVector<double> cons(const AlgebraicVector<T> &obj) {
/*! Return the constant parts.
   \return An AlgebraicVector<double>.
   \sa AlgebraicVector<T>::cons
 */
    return obj.cons();
}

template<typename T> AlgebraicVector<T> absolute(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the absolute value function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::absolute
 */
    return obj.absolute();
}

template<typename T> AlgebraicVector<T> trunc(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the truncation function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::trunc
 */
    return obj.trunc();
}

template<typename T> AlgebraicVector<T> round(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the round function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::round
 */
    return obj.round();
}

template<typename T> AlgebraicVector<T> mod(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the mod function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::mod
 */
    return obj.mod();
}

template<typename T> AlgebraicVector<T> pow(const AlgebraicVector<T> &obj, const int p) {
/*! Element-wise application of the integer power function.
   \param[in] obj AlgebraicVector<T>.
   \param[in] p power to raise each element to.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::pow
 */
    return obj.pow(p);
}

template<typename T> AlgebraicVector<T> pow(const AlgebraicVector<T> &obj, const double p) {
/*! Element-wise application of the double power function.
   \param[in] obj AlgebraicVector<T>.
   \param[in] p power to raise each element to.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::pow
 */
    return obj.pow(p);
}

template<typename T> AlgebraicVector<T> root(const AlgebraicVector<T> &obj, const int p) {
/*! Element-wise application of the root function.
   \param[in] obj AlgebraicVector<T>.
   \param[in] p root to be computed.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::root
 */
    return obj.root(p);
}

template<typename T> AlgebraicVector<T> minv(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the multiplicative inverse function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::minv
 */
    return obj.minv();
}

template<typename T> AlgebraicVector<T> sqr(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the square function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::sqr
 */
    return obj.sqr();
}

template<typename T> AlgebraicVector<T> sqrt(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the square root function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::sqrt
 */
    return obj.sqrt();
}

template<typename T> AlgebraicVector<T> isrt(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the inverse square root function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::isrt
 */
    return obj.isrt();
}

template<typename T> AlgebraicVector<T> cbrt(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the cube root function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::cbrt
 */
    return obj.cbrt();
}

template<typename T> AlgebraicVector<T> icbrt(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the inverse cube root function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::icbrt
 */
    return obj.icbrt();
}

template<typename T> AlgebraicVector<T> hypot(const AlgebraicVector<T> &X, const AlgebraicVector<T> &Y) {
/*! Element-wise application of the hypotenuse function.
   \param[in] X AlgebraicVector<T> containing X.
   \param[in] Y AlgebraicVector<T> containing Y.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::hypot
 */
    return X.hypot(Y);
}

template<typename T> AlgebraicVector<T> exp(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the exponential function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::exp
 */
    return obj.exp();
}

template<typename T> AlgebraicVector<T> log(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the natural logarithm function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::log
 */
    return obj.log();
}

template<typename T> AlgebraicVector<T> logb(const AlgebraicVector<T> &obj, const double b) {
/*! Element-wise application of the logarithm function relative to given base.
   \param[in] obj AlgebraicVector<T>.
   \param[in] b base for the logarithm
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::logb
 */
    return obj.logb(b);
}

template<typename T> AlgebraicVector<T> log10(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the decadic logarithm function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::log10
 */
    return obj.log10();
}

template<typename T> AlgebraicVector<T> log2(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the binary logarithm function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::log2
 */
    return obj.log2();
}

template<typename T> AlgebraicVector<T> sin(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the sine function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::sin
 */
    return obj.sin();
}

template<typename T> AlgebraicVector<T> cos(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the cosine function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::cos
 */
    return obj.cos();
}

template<typename T> AlgebraicVector<T> tan(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the tangent function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::tan
 */
    return obj.tan();
}

template<typename T> AlgebraicVector<T> asin(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the arcsine function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::asin
 */
    return obj.asin();
}

template<typename T> AlgebraicVector<T> acos(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the arccosine function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::acos
 */
    return obj.acos();
}

template<typename T> AlgebraicVector<T> atan(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the arctangent function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::atan
 */
    return obj.atan();
}

template<typename T> AlgebraicVector<T> atan2(const AlgebraicVector<T> &Y, const AlgebraicVector<T> &X) {
/*! Element-wise application of the four-quadrant arctangent of Y/X.
   \param[in] Y AlgebraicVector<T> containing Y.
   \param[in] X AlgebraicVector<T> containing X.
   \return A new AlgebraicVector<T> with elements in [-pi, pi].
   \sa AlgebraicVector<T>::atan2
 */
    return Y.atan2(X);
}

template<typename T> AlgebraicVector<T> sinh(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the hyperbolic sine function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::sinh
 */
    return obj.sinh();
}

template<typename T> AlgebraicVector<T> cosh(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the hyperbolic cosine function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::cosh
 */
    return obj.cosh();
}

template<typename T> AlgebraicVector<T> tanh(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the hyperbolic tangent function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::tanh
 */
    return obj.tanh();
}

template<typename T> AlgebraicVector<T> asinh(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the hyperbolic arcsine function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::asinh
 */
    return obj.asinh();
}

template<typename T> AlgebraicVector<T> acosh(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the hyperbolic arccosine function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::acosh
 */
    return obj.acosh();
}

template<typename T> AlgebraicVector<T> atanh(const AlgebraicVector<T> &obj) {
/*! Element-wise application of the hyperbolic arctangent function.
   \param[in] obj AlgebraicVector<T>.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::atanh
 */
    return obj.atanh();
}

template<typename U, typename V> typename PromotionTrait<U,V>::returnType dot(const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2) {
/*! Compute the dot product between two AlgebraicVectors.
   \param[in] obj1 a AlgebraicVector.
   \param[in] obj2 a AlgebraicVector.
   \return A scalar value,.
 */
    return obj1.dot(obj2);
}

template<typename U, typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> cross(const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2) {
/*! Compute the cross product between two 3D AlgebraicVectors.
   \param[in] obj1 a AlgebraicVector.
   \param[in] obj2 a AlgebraicVector.
   \return A new AlgebraicVector,.
 */
    return obj1.cross(obj2);
}

template<typename T> T length(const AlgebraicVector<T> &obj) {
/*! Compute the length (Euclidean vector norm).
   \param[in] obj AlgebraicVector<T>.
   \return Length of the vector.
 */
    return obj.length();
}

template<typename T> AlgebraicVector<T> normalize(const AlgebraicVector<T> &obj) {
/*! Normalize an AlgebraicVector<T>.
   \param[in] obj An AlgebraicVector<T> to normalize.
   \return An AlgebraicVector<T> of unit length.
   \sa AlgebraicVector<T>::normalize
 */
    return obj.normalize();
}

template<typename V> V eval(const AlgebraicVector<DA> &obj, const V &args) {
/*! Evaluate an AlgebraicVector<DA> with a vector type V of arguments
    and return a vector of type V with the results.
   \param[in] obj An AlgebraicVector<DA>.
   \param[in] args Vector type V containing the arguments.
   \return A new vector of type V containing the results of the evaluation process.
   \sa AlgebraicVector<T>::eval()
 */
    return obj.eval(args);
}

template<typename T> AlgebraicVector<T> eval(const AlgebraicVector<DA> &obj, const std::initializer_list<T> l) {
/*! Evaluate an AlgebraicVector<DA> with an braced initializer list of type T
    and return an AlgebraicVector of type T with the results.
   \param[in] obj An AlgebraicVector<DA>.
   \param[in] l Braced initializer list containing the arguments.
   \return A new AlgebraicVector of type T containing the results of the evaluation.
   \note C++ is not able to derive the type of elements of an initializer list automatically.
    That means eval() must be called explicitly as e.g. eval<double>(x, {1.0, 2.0, 3.0}) when
    used with initializer lists.
   \sa AlgebraicVector<T>::eval()
 */
    return obj.eval<T>(l);
}

template<typename U> AlgebraicVector<U> evalScalar(const AlgebraicVector<DA> &obj, const U &arg) {
/*! Evaluate an AlgebraicVector<DA> with a single scalar argument of type U
    and return an AlgebraicVector<U> containing the results.
   \param[in] obj The AlgebraicVector<T> to evaluate.
   \param[in] arg The argument of type U.
   \return A new AlgebraicVector<U> containing the results of the evaluation process.
   \sa AlgebraicVector<T>::evalScalar()
 */
    return obj.evalScalar(arg);
}

template<typename T> AlgebraicVector<double> norm(const AlgebraicVector<T> &obj, const unsigned int type) {
/*! Element-wise application of the norm function.
   \param[in] obj AlgebraicVector<T>.
   \param[in] type type of norm to be computed. See DA::norm.
   \return A new AlgebraicVector<T>.
   \sa AlgebraicVector<T>::norm
   \sa DA::norm
 */
    return obj.norm(type);
}


}
#endif /* DINAMICA_ALGEBRAICVECTOR_T_H_ */
