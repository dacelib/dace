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
 * DA.h
 *
 *  Created on: Feb 24, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_DA_H_
#define DINAMICA_DA_H_

// DACE C++ interface version (must match the version returned by DACEVER)
#define DACE_CPP_MAJOR (2)
#define DACE_CPP_MINOR (1)

// C++ stdlib classes used in this public interface
#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <initializer_list>

#include "dace/dacecore.h"

namespace DACE{

// forward declarations
class compiledDA;
class storedDA;
class DACEException;
class Monomial;
class Interval;
class DA;
template<typename T> class AlgebraicVector;

/*! Basic DA class representing a single polynomial. */
class DACE_API DA
{
    friend class compiledDA;
    friend class storedDA;

private:
    static bool initialized;                                                //!< Indicates if DA::init() was called
    static std::stack<unsigned int> TOstack;                                //!< Truncation order stack
    DACEDA m_index;                                                         //!< Index to the DA vector

public:
    /********************************************************************************
    *     DACE Setup
    *********************************************************************************/
    static void init(const unsigned int ord, const unsigned int nvar);      //!< DACE initialization
    static bool isInitialized();                                            //!< Get DACE initialization status
    static void version(int &maj, int &min, int &patch);                    //!< DACE core routines version
    static void checkVersion();                                             //!< Check DACE C++ interface and core routines version for compatibility
    static double setEps(const double eps);                                 //!< Set truncation epsilon
    static double getEps();                                                 //!< Get truncation epsilon
    static double getEpsMac();                                              //!< Get the machine epsilon
    static unsigned int getMaxOrder();                                      //!< Get the maximum order
    static unsigned int getMaxVariables();                                  //!< Get the maximum number of variables
    static unsigned int getMaxMonomials();                                  //!< Get the maximum monomials
    static unsigned int getTO();                                            //!< Get truncation order
    static unsigned int setTO(const unsigned int ot = DA::getMaxOrder());   //!< Set truncation order
    static void pushTO(const unsigned int ot = DA::getMaxOrder());          //!< Set truncation order and save previous value
    static void popTO();                                                    //!< Restore truncation order to value before last DA::pushTO()

    /********************************************************************************
    *     Constructors & Destructors
    *********************************************************************************/
    DA();                                                                   //!< Default constructor
    DA(const DA &da);                                                       //!< Copy constructor

    DA(DA &&da);                                                            //!< Move constructor

    explicit DA(const int i, const double c = 1.0);                         //!< Constructor for DA identities
    explicit DA(const unsigned int i, const double c = 1.0);                //!< Constructor for DA identities
    DA(const double c);                                                     //!< Constructor for DA constants
    ~DA() throw();                                                          //!< Destructor

    /********************************************************************************
    *     Coefficient access and extraction routines
    *********************************************************************************/
    int isnan() const;
    int isinf() const;
    double cons() const;                                                    //!< Get constant part of a DA
    AlgebraicVector<double> linear() const;                                 //!< Get linear part of a DA
    AlgebraicVector<DA> gradient() const;                                   //!< Gradient vector with respect to all independent DA variables
    double getCoefficient(const std::vector<unsigned int> &jj) const;                //!< Get specific coefficient
    void setCoefficient(const std::vector<unsigned int> &jj, const double coeff);    //!< Set specific coefficient
    Monomial getMonomial(const unsigned int npos) const;                    //!< Get the Monomial at given position
    void getMonomial(const unsigned int npos, Monomial &m) const;            //!< Extract the Monomial at given position
    std::vector<Monomial> getMonomials() const;                             //!< Get std::vector of all non-zero Monomials

    /********************************************************************************
    *     Assignments
    *********************************************************************************/
    DA& operator=(DA &&da);                                                 //!< Move assignment from existing DA

    DA& operator=(const DA &da);                                            //!< Assignment from existing DA
    DA& operator=(const double c);                                          //!< Assignment from existing double

    DA& operator+=(const DA &da);                                           //!< Add another DA to this DA
    DA& operator+=(const double c);                                         //!< Add another double to this DA

    DA& operator-=(const DA &da);                                           //!< Subtract another DA from this DA
    DA& operator-=(const double c);                                         //!< Subtract another double from this DA

    DA& operator*=(const DA &da);                                           //!< Multiply this DA by another DA
    DA& operator*=(const double c);                                         //!< Multiply this DA by another double

    DA& operator/=(const DA &da);                                           //!< Divide this DA by another DA
    DA& operator/=(const double c);                                         //!< Divide this DA by another double

    /********************************************************************************
    *     Algebraic operations
    *********************************************************************************/
    DA operator-() const;                                                   //!< Negate this DA

    friend DA DACE_API operator+(const DA &da1, const DA &da2);             //!< Addition between two DAs
    friend DA DACE_API operator+(const DA &da, const double c);             //!< Addition between a DA and a constant
    friend DA DACE_API operator+(const double c, const DA &da);             //!< Addition between a constant and a DA

    friend DA DACE_API operator-(const DA &da1, const DA &da2);             //!< Subtraction between two DAs
    friend DA DACE_API operator-(const DA &da, const double c);             //!< Subtraction between a DA and a constant
    friend DA DACE_API operator-(const double c, const DA &da);             //!< Subtraction between a constant and a DA

    friend DA DACE_API operator*(const DA &da1, const DA &da2);             //!< Multiplication between two DAs
    friend DA DACE_API operator*(const DA &da, const double c);             //!< Multiplication between a DA and a constant
    friend DA DACE_API operator*(const double c, const DA &da);             //!< Multiplication between a constant and a DA

    friend DA DACE_API operator/(const DA &da1, const DA &da2);             //!< Division between two DAs
    friend DA DACE_API operator/(const DA &da, const double c);             //!< Division between a DA and a constant
    friend DA DACE_API operator/(const double c, const DA &da);             //!< Division between a constant and a DA

    /********************************************************************************
    *     Math routines
    *********************************************************************************/
    DA multiplyMonomials(const DA &da) const;                                //!< Multiply the DA with the argument da monomial by monomial (i.e. coefficient-wise)
    DA divide(const unsigned int var, const unsigned int p = 1) const;      //!< Divide by an independent variable to some power
    DA deriv(const unsigned int i) const;                                   //!< Derivative with respect to given variable
    DA deriv(const std::vector<unsigned int> ind) const;                    //!< Derivative with respect to given variables
    DA integ(const unsigned int i) const;                                   //!< Integral with respect to given variable
    DA integ(const std::vector<unsigned int> ind) const;                    //!< Integral with respect to given variables
    DA trim(const unsigned int min, const unsigned int max = DA::getMaxOrder()) const;
                                                                            //!< Trim the coefficients only include orders between min and max, inclusively
    DA trunc() const;                                                       //!< Truncate the constant part to an integer
    DA round() const;                                                       //!< Round the constant part to an integer
    DA mod(const double p) const;                                           //!< Modulo of the constant part
    DA pow(const int p) const;                                              //!< Exponentiation to given (integer) power
    DA pow(const double p) const;                                           //!< Exponentiation to given double power
    DA root(const int p = 2) const;                                         //!< p-th root
    DA minv() const;                                                        //!< multiplicative inverse
    DA sqr() const;                                                         //!< square
    DA sqrt() const;                                                        //!< square root
    DA isrt() const;                                                        //!< inverse square root
    DA cbrt() const;                                                        //!< cubic root
    DA icrt() const;                                                        //!< inverse cubic root
    DA hypot(const DA &da) const;                                           //!< hypotenuse
    DA exp() const;                                                         //!< exponential
    DA log() const;                                                         //!< natural logarithm
    DA logb(const double b = 10.0) const;                                   //!< logarithm with respect to given base
    DA log10() const;                                                       //!< 10 based logarithm
    DA log2() const;                                                        //!< 2 based logarithm
    DA sin() const;                                                         //!< sine
    DA cos() const;                                                         //!< cosine
    DA tan() const;                                                         //!< tangent
    DA asin() const;                                                        //!< arcsine
    DA acos() const;                                                        //!< arccosine
    DA atan() const;                                                        //!< arctangent
    DA atan2(const DA &da) const;                                           //!< arctangent in [-pi, pi]
    DA sinh() const;                                                        //!< hyperbolic sine
    DA cosh() const;                                                        //!< hyperbolic cosine
    DA tanh() const;                                                        //!< hyperbolic tangent
    DA asinh() const;                                                       //!< hyperbolic arcsine
    DA acosh() const;                                                       //!< hyperbolic arccosine
    DA atanh() const;                                                       //!< hyperbolic arctangent
    DA erf() const;                                                         //!< error function
    DA erfc() const;                                                        //!< complementary error function
    DA BesselJFunction(const int n) const;                                  //!< Bessel function J_n
    DA BesselYFunction(const int n) const;                                  //!< Bessel function Y_n
    DA BesselIFunction(const int n, const bool scaled = false) const;       //!< Bessel function I_n
    DA BesselKFunction(const int n, const bool scaled = false) const;       //!< Bessel function K_n
    DA GammaFunction() const;                                               //!< Gamma function
    DA LogGammaFunction() const;                                            //!< Logartihmic Gamma function
    DA PsiFunction(const unsigned int n) const;                             //!< Psi function

    /********************************************************************************
    *    Norm and estimation routines
    *********************************************************************************/
    unsigned int size() const;                                              //!< Number of non-zero coefficients
    double abs() const;                                                     //!< Maximum absolute value of all coefficients
    double norm(const unsigned int type = 0) const;                         //!< Different types of norms over all coefficients
    std::vector<double> orderNorm(const unsigned int var = 0, const unsigned int type = 0) const;
                                                                            //!< Different types of norms over coefficients of each order separately
    std::vector<double> estimNorm(const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = DA::getMaxOrder()) const;
                                                                            //!< Estimate of different types of order sorted norms
    std::vector<double> estimNorm(std::vector<double> &err, const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = DA::getMaxOrder()) const;
                                                                            //!< Estimate of different types of order sorted norms with error of the estimate
    Interval bound() const;                                                 //!< Estimate range bound over [-1,1] for each independent variable
    double convRadius(const double eps, const unsigned int type = 1) const; //!< Estimate the convergence radius of the current DA

    /********************************************************************************
    *     DACE polynomial evaluation routines
    *********************************************************************************/
    template<class T> T eval(const std::vector<T> &args) const;             //!< Evaluation with a vector of arguments (not efficient for repeated evaluation!)
    template<class T> T eval(const T args[], const unsigned int length) const; //!< Evaluation with an array of arguments (not efficient for repeated evaluation!)
    template<class T> T evalScalar(const T &arg) const;                     //!< Evaluation with a single arithmetic type T argument (not efficient for repeated evaluation!)
    compiledDA compile() const;                                             //!< Compile current DA for efficient repeated evaluation
    DA plug(const unsigned int var, const double val = 0.0) const;          //!< Partial evaluation to replace given independent DA variable by value val
    double evalMonomials(const DA &values) const;                           //!< evaluate the DA providing the value of every monomial in da
    DA replaceVariable(const unsigned int from = 0, const unsigned int to = 0, const double val = 1.0) const;
                                                                            //!< Replace variable number from by val times variable number to
    DA scaleVariable(const unsigned int var = 0, const double val = 1.0) const;
                                                                            //!< Scale variable number var by val
    DA translateVariable(const unsigned int var = 0, const double a = 1.0, const double c = 0.0) const;
                                                                            //!< Translate variable number var to a*var + c

    /********************************************************************************
    *     DACE input/output routines
    *********************************************************************************/
    std::string toString() const;                                           //!< Convert to string representation
    void write(std::ostream &os) const;                                     //!< Write binary representation of DA to stream
    friend DACE_API std::ostream& operator<< (std::ostream &out, const DA &da);      //!< Output to C++ stream in text form
    friend DACE_API std::istream& operator>> (std::istream &in, DA &da);             //!< Input from C++ stream in text form

    /********************************************************************************
    *     Static factory routines
    *********************************************************************************/
    static DA random(const double cm);                                      //!< Create new DA filled with random coefficients
    static DA identity(const unsigned int var);                             //!< Create the DA identity for independent DA variable var
    static DA fromString(const std::string &str);                           //!< Create new DA from string
    static DA fromString(const std::vector<std::string> &str);              //!< Create new DA from vector of strings
    static DA read(std::istream &is);                                       //!< Create new DA from blob read from binary file

    /********************************************************************************
    *     DACE various routines
    *********************************************************************************/
    static void memdump();                                                  //!< Dump memory status to screen (for debugging purposes only)
};

/********************************************************************************
*     DACE non-member functions
*********************************************************************************/
DACE_API int isnan(const DA &da);
DACE_API int isinf(const DA &da);
DACE_API double cons(const DA &da);
DACE_API AlgebraicVector<double> linear(const DA &da);
DACE_API AlgebraicVector<DA> gradient(const DA &da);

DACE_API DA divide(const DA &da, const unsigned int var, const unsigned int p = 1);
DACE_API DA deriv(const DA &da, const unsigned int i);
DACE_API DA deriv(const DA &da, const std::vector<unsigned int> ind);
DACE_API DA integ(const DA &da, const unsigned int i);
DACE_API DA integ(const DA &da, const std::vector<unsigned int> ind);
DACE_API DA trim(const DA &da, const unsigned int min, const unsigned int max = DA::getMaxOrder());
DACE_API DA trunc(const DA &da);
DACE_API DA round(const DA &da);
DACE_API DA mod(const DA &da, const double p);
DACE_API DA pow(const DA &da, const int p);
DACE_API DA pow(const DA &da, const double p);
DACE_API DA root(const DA &da, const int p = 2);
DACE_API DA minv(const DA &da);
DACE_API DA sqr(const DA &da);
DACE_API DA sqrt(const DA &da);
DACE_API DA isrt(const DA &da);
DACE_API DA cbrt(const DA &da);
DACE_API DA icrt(const DA &da);
DACE_API DA hypot(const DA &da1, const DA &da2);
DACE_API DA exp(const DA &da);
DACE_API DA log(const DA &da);
DACE_API DA logb(const DA &da, const double b = 10.0);
DACE_API DA log10(const DA &da);
DACE_API DA log2(const DA &da);
DACE_API DA sin(const DA &da);
DACE_API DA cos(const DA &da);
DACE_API DA tan(const DA &da);
DACE_API DA asin(const DA &da);
DACE_API DA acos(const DA &da);
DACE_API DA atan(const DA &da);
DACE_API DA atan2(const DA &da1, const DA &da2);
DACE_API DA sinh(const DA &da);
DACE_API DA cosh(const DA &da);
DACE_API DA tanh(const DA &da);
DACE_API DA asinh(const DA &da);
DACE_API DA acosh(const DA &da);
DACE_API DA atanh(const DA &da);
DACE_API DA erf(const DA &da);
DACE_API DA erfc(const DA &da);
DACE_API DA jn(const int n, const DA &da);
DACE_API DA yn(const int n, const DA &da);
DACE_API DA BesselJFunction(const int n, const DA &da);
DACE_API DA BesselYFunction(const int n, const DA &da);
DACE_API DA BesselIFunction(const int n, const DA &da, const bool scaled = false);
DACE_API DA BesselKFunction(const int n, const DA &da, const bool scaled = false);
DACE_API DA tgamma(const DA &da);
DACE_API DA lgamma(const DA &da);
DACE_API DA GammaFunction(const DA &da);
DACE_API DA LogGammaFunction(const DA &da);
DACE_API DA PsiFunction(const unsigned int n, const DA &da);

DACE_API unsigned int size(const DA &da);
DACE_API double abs(const DA &da);
DACE_API double norm(const DA &da, unsigned int type = 0);
DACE_API std::vector<double> orderNorm(const DA &da, unsigned int var = 0, unsigned int type = 0);
DACE_API std::vector<double> estimNorm(const DA &da, unsigned int var = 0, unsigned int type = 0, unsigned int nc = DA::getMaxOrder());
DACE_API std::vector<double> estimNorm(const DA &da, std::vector<double> &err, unsigned int var = 0, unsigned int type = 0, unsigned int nc = DA::getMaxOrder());
DACE_API Interval bound(const DA &da);
DACE_API double convRadius(const DA &da, const double eps, const unsigned int type = 1);

template<class T> T eval(const DA &da, const std::vector<T> &args);
template<class T> T eval(const DA &da, const T args[], const unsigned int length);
template<class T> T evalScalar(const DA &da, const T &arg);
DACE_API compiledDA compile(const DA &da);
DACE_API DA plug(const DA &da, const unsigned int var, const double val = 0.0);
DACE_API DA replaceVariable(const DA &da, const unsigned int from = 0, const unsigned int to = 0, const double val = 1.0);
DACE_API DA scaleVariable(const DA &da, const unsigned int var = 0, const double val = 1.0);
DACE_API DA translateVariable(const DA &da, const unsigned int var = 0, const double a = 1.0, const double c = 0.0);

DACE_API std::string toString(const DA &da);
DACE_API void write(const DA &da, std::ostream &os);



/*! Stored DA class representing a DA vector in a binary, setup independent format. */
class DACE_API storedDA : std::vector<char>
{
private:
    static const unsigned int headerSize;

public:
    storedDA(const DA &da);                            //!< Constructor from a DA.
    storedDA(const std::vector<char> &data);        //!< Constructor from binary data.
    storedDA(std::istream &is);                     //!< Constructor from stream.

    bool isValid() const;                           //!< Is the data a valid DACE blob

    operator DA() const;                            //!< Cast to DA
    operator std::string() const;                   //!< Cast to std::string

    friend DACE_API std::ostream& operator<<(std::ostream &out, const storedDA &sda);      //!< Output to C++ stream in binary form
};

}

#endif /* DINAMICA_DA_H_ */
