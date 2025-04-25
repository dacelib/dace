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
 * AlgebraicVector.h
 *
 *  Created on: Sep. 10, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_ALGEBRAICVECTOR_H_
#define DINAMICA_ALGEBRAICVECTOR_H_

// C++ stdlib classes required for interface definition
#include <vector>
#include <initializer_list>

// DACE classes required for interface definition (DA.h needed for DA::getMaxOrder(), DA::getMaxVariables() default arguments)
#include "dace/PromotionTrait.h"
#include "dace/DA.h"

namespace DACE{

// forward declarations
#ifdef WITH_ALGEBRAICMATRIX
template<typename T> class AlgebraicMatrix;
#endif

/*! Generic vector class to handle vectors of algebraic types and their algebraic operations. */
template<typename T> class AlgebraicVector : public std::vector<T>
{
public:
    /***********************************************************************************
    *     Constructors
    ************************************************************************************/
    AlgebraicVector();                                                                        //!< Default constructor
    explicit AlgebraicVector(const size_t size);                                              //!< Constructor with size
    AlgebraicVector(const size_t size, const T &d);                                           //!< Constructor with size and elements value
    AlgebraicVector(const std::vector<T> &v);                                                 //!< Copy constructor
    AlgebraicVector(const std::vector<T> &v, const size_t first, const size_t last);          //!< Extraction constructor
    AlgebraicVector(std::initializer_list<T> l);                                              //!< Constructor from braced initializer list

    /***********************************************************************************
    *     Element and coefficient access / extraction routines
    ************************************************************************************/
    AlgebraicVector<T> extract(const size_t first, const size_t last) const;            //!< Return the subvector containing the elements between first and last, inclusively
    template<typename U> AlgebraicVector<typename PromotionTrait<T,U>::returnType> concat(const std::vector<U> &obj) const;
                                                                                                    //!< Return a new vector containing the elements of this vector followed by those of obj
    AlgebraicVector<double> cons() const;                                                           //!< Return vector containing only the costant parts of each element
#ifdef WITH_ALGEBRAICMATRIX
    AlgebraicMatrix<double> linear() const;                                                         //!< Return the linear parts in the form of a Matrix
#else
    std::vector< std::vector<double> > linear() const;                                              //!< Return the linear parts in the form of a vector of vectors
#endif /* WITH_ALGEBRAICMATRIX */

    /***********************************************************************************
    *     Operator overloads
    ************************************************************************************/
    AlgebraicVector<T> operator-() const;
    template<typename U> AlgebraicVector<T>& operator+=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator+=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator-=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator-=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator*=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator*=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator/=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator/=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator<<(const std::vector<U> &obj);                 //!< Concatenation operator.

    /***********************************************************************************
    *     Math routines
    ************************************************************************************/
    AlgebraicVector<T> absolute() const;                                                            //!< Element-wise absolute value function.
    AlgebraicVector<T> trunc() const;                                                               //!< Element-wise truncation.
    AlgebraicVector<T> round() const;                                                               //!< Element-wise rounding.
    AlgebraicVector<T> mod() const;                                                                 //!< Element-wise modulo.
    AlgebraicVector<T> pow(const int p) const;                                                      //!< Element-wise exponentiation to (integer) power.
    AlgebraicVector<T> pow(const double p) const;                                                   //!< Element-wise exponentiation to (double) power.
    AlgebraicVector<T> root(const int p = 2) const;                                                 //!< Element-wise p-th root.
    AlgebraicVector<T> minv() const;                                                                //!< Element-wise multiplicative inverse.
    AlgebraicVector<T> sqr() const;                                                                 //!< Element-wise square.
    AlgebraicVector<T> sqrt() const;                                                                //!< Element-wise square root.
    AlgebraicVector<T> isrt() const;                                                                //!< Element-wise inverse square root.
    AlgebraicVector<T> cbrt() const;                                                                //!< Element-wise cube root.
    AlgebraicVector<T> icbrt() const;                                                               //!< Element-wise inverse cube root.
    AlgebraicVector<T> hypot(const AlgebraicVector<T> &obj) const;                                  //!< Element-wise hypotenuse.
    AlgebraicVector<T> exp() const;                                                                 //!< Element-wise exponential.
    AlgebraicVector<T> log() const;                                                                 //!< Element-wise natural logarithm.
    AlgebraicVector<T> logb(const double b = 10.0) const;                                           //!< Element-wise logarithm wrt a given base.
    AlgebraicVector<T> log10() const;                                                               //!< Element-wise logarithm to base 10.
    AlgebraicVector<T> log2() const;                                                                //!< Element-wise logarithm to base 2.
    AlgebraicVector<T> sin() const;                                                                 //!< Element-wise sine.
    AlgebraicVector<T> cos() const;                                                                 //!< Element-wise cosine.
    AlgebraicVector<T> tan() const;                                                                 //!< Element-wise tangent.
    AlgebraicVector<T> asin() const;                                                                //!< Element-wise arcsine.
    AlgebraicVector<T> acos() const;                                                                //!< Element-wise arccosine.
    AlgebraicVector<T> atan() const;                                                                //!< Element-wise arctangent.
    AlgebraicVector<T> atan2(const AlgebraicVector<T> &obj) const;                                  //!< Element-wise arctangent in [-pi, pi].
    AlgebraicVector<T> sinh() const;                                                                //!< Element-wise hyperbolic sine.
    AlgebraicVector<T> cosh() const;                                                                //!< Element-wise hyperbolic cosine.
    AlgebraicVector<T> tanh() const;                                                                //!< Element-wise hyperbolic tangent.
    AlgebraicVector<T> asinh() const;                                                               //!< Element-wise hyperbolic arcsine.
    AlgebraicVector<T> acosh() const;                                                               //!< Element-wise hyperbolic arccosine.
    AlgebraicVector<T> atanh() const;                                                               //!< Element-wise hyperbolic arctangent.

    /***********************************************************************************
    *    Vector routines
    ************************************************************************************/
    template<typename V> typename PromotionTrait<T,V>::returnType dot(const AlgebraicVector<V> &obj) const;
                                                                                                    //!< Dot product (scalar product, inner product) of two vectors.
    template<typename V> AlgebraicVector<typename PromotionTrait<T,V>::returnType> cross(const AlgebraicVector<V> &obj) const;
                                                                                                    //!< Cross product of two vectors of length 3.
    T length() const;                                                                               //!< Length of the vector in Euclidean norm.
    AlgebraicVector<T> normalize() const;                                                           //!< Normalized vector of unit length along this vector.
    // XXX: various Jacobians, gradients, curls, etc?

    /***********************************************************************************
    *     Special routines (DA related)
    ************************************************************************************/
    AlgebraicVector<T> deriv(const unsigned int p) const;                                           //!< Derivative of each element with respect to given variable. DA only.
    AlgebraicVector<T> integ(const unsigned int p) const;                                           //!< Integration of each element with respect to given variable. DA only.
    template<typename V> V eval(const V &args) const;                                               //!< Generic evaluation of a AlgebraicVector<DA> with arguments. DA only.
    template<typename U> AlgebraicVector<U> eval(const std::initializer_list<U> l) const;           //!< Generic evaluation of an AlgebraicVector<DA> with braced initializer list. DA only.
    template<typename U> AlgebraicVector<U> evalScalar(const U &arg) const;                         //!< Generic evaluation of a AlgebraicVector<DA> with single argument.  DA only.
    compiledDA compile() const;                                                                     //!< Compile current DA for efficient repeated evaluation. DA only.
    AlgebraicVector<T> plug(const unsigned int var, const double val = 0.0) const;                  //!< Partial evaluation to replace given independent DA variable by value val. DA only.
    AlgebraicVector<T> trim(const unsigned int min, const unsigned int max = DA::getMaxOrder()) const;
                                                                                                    //!< Trim the coefficients of each components to particular orders. DA only.
    AlgebraicVector<T> invert() const;                                                              //!< Inverse function of the AlgebraicVector<DA>. DA only.

    /********************************************************************************
    *     DA norm routines
    *********************************************************************************/
    AlgebraicVector<double> norm(const unsigned int type = 0) const;                                //!< Element-wise DA norm
    /* XXX: define and add the norm estimation routines from DA including convergence radius estimation
    std::vector<double> orderNorm(const unsigned int var = 0, const unsigned int type = 0) const;
                                                                            //!< Different types of norms over coefficients of each order separately
    std::vector<double> estimNorm(const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = DA::getMaxOrder()) const;
                                                                            //!< Estimate of different types of order sorted norms
    Interval bound() const;                                                 //!< Estimate range bound over [-1,1] for each independent variable
    double convRadius(const double eps, const unsigned int type = 1) const; //!< Estimate the convergence radius of the current DA.
    */

    /********************************************************************************
    *     Static factory routines
    *********************************************************************************/
    static AlgebraicVector<DA> identity(const size_t n = DA::getMaxVariables());              //!< Create an AlgebraicVector<DA> containing the identity in n dimensions. DA only.

    /***********************************************************************************
    *     Input/Output routines
    ************************************************************************************/
    std::string toString() const;                                                                   //!< Convert the vector into a human readable string.

private:
#ifndef WITH_ALGEBRAICMATRIX
    static void matrix_inverse(std::vector< std::vector<double> > &A);        // Private helper routine for double precision matrix inversion
#endif /* WITH_ALGEBRAICMATRIX */
};

// operators
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator+( const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator+( const AlgebraicVector<U> &obj1, const V &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator+( const U &obj1, const AlgebraicVector<V> &obj2);

template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator-( const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator-( const AlgebraicVector<U> &obj1, const V &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator-( const U &obj1, const AlgebraicVector<V> &obj2);

template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator*( const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator*( const AlgebraicVector<U> &obj1, const V &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator*( const U &obj1, const AlgebraicVector<V> &obj2);

template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator/( const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator/( const AlgebraicVector<U> &obj1, const V &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator/( const U &obj1, const AlgebraicVector<V> &obj2);

template<typename U> std::ostream& operator<<(std::ostream &out, const AlgebraicVector<U> &obj);    //!< Overload output stream operator.
template<typename U> std::istream& operator>>(std::istream &in, AlgebraicVector<U> &obj);           //!< Overload input stream operator.

// Declaration of external functional style wrappers to access AlgebraicVector functions
template<typename T> AlgebraicVector<double> cons(const AlgebraicVector<T> &obj);
#ifdef WITH_ALGEBRAICMATRIX
template<typename T> AlgebraicMatrix<double> linear(const AlgebraicVector<T> &obj);
#else
template<typename T> std::vector< std::vector<double> > linear(const AlgebraicVector<T> &obj);
#endif /* WITH_ALGEBRAICMATRIX */
template<typename T> AlgebraicVector<T> deriv(const AlgebraicVector<T> &obj, const unsigned int p);
template<typename T> AlgebraicVector<T> integ(const AlgebraicVector<T> &obj, const unsigned int p);
template<typename T> AlgebraicVector<T> pow(const AlgebraicVector<T> &obj, const int p);
template<typename T> AlgebraicVector<T> root(const AlgebraicVector<T> &obj, const int p = 2);
template<typename T> AlgebraicVector<T> minv(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> sqr(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> sqrt(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> isrt(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> exp(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> log(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> logb(const AlgebraicVector<T> &obj, const double b = 10.0);
template<typename T> AlgebraicVector<T> sin(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> cos(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> tan(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> asin(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> acos(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> atan(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> atan2(const AlgebraicVector<T> &obj1, const AlgebraicVector<T> &obj2);
template<typename T> AlgebraicVector<T> sinh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> cosh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> tanh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> asinh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> acosh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> atanh(const AlgebraicVector<T> &obj);
template<typename U, typename V> typename PromotionTrait<U,V>::returnType dot(const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename U, typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> cross(const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename T> T vnorm(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> normalize(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> trim(const AlgebraicVector<T> &obj, unsigned int min, unsigned int max = DA::getMaxOrder());
template<typename T, typename V> V eval(const AlgebraicVector<T> &obj, const V &args);
template<typename T, typename U> AlgebraicVector<U> eval(const AlgebraicVector<T> &obj, const std::initializer_list<U> l);
template<typename T, typename U> AlgebraicVector<U> evalScalar(const AlgebraicVector<T> &obj, const U &arg);
template<typename T> compiledDA compile(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> plug(const AlgebraicVector<T> &obj, const unsigned int var, const double val = 0.0);
template<typename T> AlgebraicVector<double> norm(const AlgebraicVector<T> &obj, const unsigned int type = 0);

// specializations for various DA specific routines implemented and instantiated directly in the library instead of in a template
#ifdef WITH_ALGEBRAICMATRIX
template<> DACE_API AlgebraicMatrix<double> AlgebraicVector<DA>::linear() const;
template<> DACE_API AlgebraicMatrix<double> linear(const AlgebraicVector<DA> &obj);
#else
template<> DACE_API std::vector< std::vector<double> > AlgebraicVector<DA>::linear() const;
template<> DACE_API void AlgebraicVector<DA>::matrix_inverse(std::vector< std::vector<double> > &A);
template<> DACE_API std::vector< std::vector<double> > linear(const AlgebraicVector<DA> &obj);
#endif /* WITH_ALGEBRAICMATRIX */
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::trim(const unsigned int min, const unsigned int max) const;
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::deriv(const unsigned int p) const;
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::integ(const unsigned int p) const;
template<> DACE_API compiledDA AlgebraicVector<DA>::compile() const;
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::plug(const unsigned int var, const double val) const;
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::invert() const;
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::identity(const size_t n);
template<> DACE_API AlgebraicVector<DA> trim(const AlgebraicVector<DA> &obj, unsigned int min, unsigned int max);
template<> DACE_API AlgebraicVector<DA> deriv(const AlgebraicVector<DA> &obj, const unsigned int p);
template<> DACE_API AlgebraicVector<DA> integ(const AlgebraicVector<DA> &obj, const unsigned int p);
template<> DACE_API compiledDA compile(const AlgebraicVector<DA> &obj);
template<> DACE_API AlgebraicVector<DA> plug(const AlgebraicVector<DA> &obj, const unsigned int var, const double val);

// shortcuts for common vector types
typedef AlgebraicVector<DA> vectorDA;       //!< Shorthand notation for AlgebraicVector<DA>.
typedef AlgebraicVector<double> vectordb;   //!< Shorthand notation for AlgebraicVector<double>.
}

#endif /* DINAMICA_ALGEBRAICVECTOR_H_ */
