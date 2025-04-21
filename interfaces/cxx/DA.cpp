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
 * DA.cpp
 *
 *  Created on: Feb 24, 2014
 *      Author: Dinamica Srl
 */

// C++ stdlib classes used only internally in the implementation
#include <algorithm>
#include <exception>
#include <sstream>
#include <cmath>

// DACE classes
#include "dace/config.h"
#include "dace/compiledDA.h"
#include "dace/DACEException.h"
#include "dace/Monomial.h"
#include "dace/Interval.h"
#include "dace/DA.h"
#include "dace/AlgebraicVector.h"
#include "dace/AlgebraicVector_t.h"

namespace DACE{

/********************************************************************************
*     DACE static class variables
*********************************************************************************/
bool DA::initialized = false;           // not initialized yet
std::stack<unsigned int> DA::TOstack;   // truncation order stack, initially empty

/********************************************************************************
*     DACE Setup
*********************************************************************************/
void DA::init(const unsigned int ord, const unsigned int nvar) {
/*! Initialize the DACE control arrays and set the maximum order and the
    maximum number of variables.
   \note MUST BE CALLED BEFORE ANY OTHER DA ROUTINE CAN BE USED.
   \note This routine performs a mandatory version check to compare the version
    of the C++ interface used to compile the program to the version of the
    DACE library that is linked dynamically at runtime.
   \param[in] ord order of the Taylor polynomials;
   \param[in] nvar number of variables considered.
   \sa DA::checkVersion
 */
    try {
        checkVersion();
    } catch(DACEException const& ex) {
        std::cerr << ex << std::endl;
        std::terminate(); }
    daceInitialize(ord,nvar);
    if(daceGetError()) DACEException();
    initialized = true;
}

bool DA::isInitialized() {
/*! Returns the inizialisation status of the DACE.
   \return Returns true if the DACE has previously been initialized by a call
    to DA::init, or false otherwise.
   \sa DA::init
 */
    return initialized;
}

void DA::version(int &maj, int &min, int &patch) {
/*! Return the major and minor version number of the DACE
    along with the patch level of the library.
   \param[out] maj major DACE version number;
   \param[out] min minor DACE version number;
   \param[out] patch patch level of DACE version.
   \sa DA::checkVersion
 */
    daceGetVersion(maj, min, patch);
    if(daceGetError()) DACEException();
}

void DA::checkVersion() {
/*! Check the DACE core library version linked to this C++ interface
    against the interface version and throw an exception if the versions don't
    match.
   \throw DACE::DACEException
   \note This routine is called automatically by DA::init
    to ensure compatibility with the current runtime environment.
 */
    int maj, min, patch;
    version(maj, min, patch);
    if((maj!=DACE_CPP_MAJOR)||(min!=DACE_CPP_MINOR)||(maj!=DACE_MAJOR_VERSION)||(min!=DACE_MINOR_VERSION)) DACEException(20,99);
}

unsigned int DA::getMaxOrder() {
/*! Return the maximum order currently set for the computations.
   \return The maximum order, or zero if undefined.
   \throw DACE::DACEException
 */

    const unsigned int ord = daceGetMaxOrder();
    if(daceGetError()) DACEException();

    return ord;
}

double DA::setEps(const double eps) {
/*! Set the cutoff value eps to a new value and return the
    previous value.
   \param[in] eps new cutoff value.
   \return The previous cutoff value of eps, or zero if undefined.
   \throw DACE::DACEException
 */
    const double old_eps = daceSetEpsilon(eps);
    if(daceGetError()) DACEException();

    return old_eps;
}

double DA::getEps() {
/*! Return the cutoff value eps currently set for the computations.
   \return The cutoff value of eps, or zero if undefined.
   \throw DACE::DACEException
 */
    const double eps = daceGetEpsilon();
    if(daceGetError()) DACEException();

    return eps;
}

double DA::getEpsMac() {
/*! Return the machine epsilon (pessimistic estimate).
   \return The machine epsilon, or zero if undefined.
   \throw DACE::DACEException
 */
    const double epsmac = daceGetMachineEpsilon();
    if(daceGetError()) DACEException();

    return epsmac;
}

unsigned int DA::getMaxVariables() {
/*! Return the maximum number of variables set for the computations.
   \return The maximum number of variables, or zero if undefined.
   \throw DACE::DACEException
 */
    const unsigned int nvar = daceGetMaxVariables();
    if(daceGetError()) DACEException();

    return nvar;
}

unsigned int DA::getMaxMonomials() {
/*! Return the maximum number of monomials available with the
    order and number of variables specified.
   \return The maximum number of monomials, or zero if undefined.
   \throw DACE::DACEException
 */
    const unsigned int nmmax = daceGetMaxMonomials();
    if(daceGetError()) DACEException();

    return nmmax;
}

unsigned int DA::setTO(const unsigned int ot) {
/*! Set the truncation order ot to a new value and return the
    previous value. All terms larger than the truncation order are discarded
    in subsequent operations.
   \param[in] ot new truncation order.
   \return previous truncation order, or zero if undefined.
   \throw DACE::DACEException
   \sa DA::getTO
   \sa DA::pushTO
   \sa DA::popTO
 */
    const unsigned int old_no = daceSetTruncationOrder(ot);
    if(daceGetError()) DACEException();

    return old_no;
}

unsigned int DA::getTO() {
/*! Return the truncation order currently set for the computations. All terms
    larger than the truncation order are discarded in subsequent operations.
   \return current truncation order, or zero if undefined.
   \throw DACE::DACEException
   \sa DA::setTO
   \sa DA::pushTO
   \sa DA::popTO
 */
    const unsigned int no = daceGetTruncationOrder();
    if(daceGetError()) DACEException();

    return no;
}

void DA::pushTO(const unsigned int ot) {
/*! Set a new truncation order, saving the previous one on the truncation
    order stack. All terms larger than the truncation order are discarded
    in subsequent operations.
   \param[in] ot new truncation order.
   \throw DACE::DACEException
   \sa DA::getTO
   \sa DA::setTO
   \sa DA::popTO
 */
    TOstack.push(daceSetTruncationOrder(ot));
    if(daceGetError()) DACEException();
}

void DA::popTO() {
/*! Restore the previous truncation order from the truncation order stack.
    All terms larger than the truncation order are discarded in subsequent
    operations.
   \throw DACE::DACEException
   \sa DA::getTO
   \sa DA::setTO
   \sa DA::pushTO
 */
    if(!TOstack.empty()){
        daceSetTruncationOrder(TOstack.top());
        TOstack.pop();
        if(daceGetError()) DACEException();}
}

/********************************************************************************
*     Constructors & Destructors
*********************************************************************************/
DA::DA(){
/*! Create an empty DA object representing the constant zero function.
   \throw DACE::DACEException
 */
    daceAllocateDA(m_index, 0);
    if(daceGetError()) DACEException();
}

DA::DA(const DA &da){
/*! Create a copy of a DA object.
   \param[in] da object to be copied.
   \throw DACE::DACEException
 */
    daceAllocateDA(m_index, 0);
    daceCopy(da.m_index, m_index);
    if(daceGetError()) DACEException();
}

DA::DA(DA &&da){
/*! Move all resources from a DA object to us without copying.
   \param[in] da object to be moved.
   \throw DACE::DACEException
 */
    m_index = da.m_index;
    daceInvalidateDA(da.m_index);   // prevent other DA from freeing our memory
    if(daceGetError()) DACEException();
}

DA::DA(const double c){
/*! Create a DA object with the constant part equal to c.
   \note This routine MUST be called with a floating point type as the first argument, e.g. DA(1.0). Expressions involving integer data types such as DA(1) will be interpreted as the first independent DA variable instead of the constant DA object of the given value.
   \param[in] c A double value set as constant part of the DA object.
   \throw DACE::DACEException
 */
    daceAllocateDA(m_index, 0);
    daceCreateConstant(m_index, c);
    if(daceGetError()) DACEException();
}

DA::DA(const unsigned int i, const double c){
/*! Create a DA object as c times the independent variable number i.
   \note When used in its one argument form (with the default argument 1.0 for c), this routine MUST be called with an integer type as the first argument, e.g. DA(1). Expressions involving floating point data types such as DA(1.0) will be interpreted as the constant DA of the given value instead of the first independent variable.
   \param[in] i independent variable number (i=0 means the constant part).
   \param[in] c coefficient corresponding to the given independent variable. By default, this value is assumed to be 1.0.
   \throw DACE::DACEException
 */
    daceAllocateDA(m_index, 0);
    daceCreateVariable(m_index, i, c);
    if(daceGetError()) DACEException();
}

DA::DA(const int i, const double c){
/*! Create a DA object as c times the independent variable number i.
   \note When used in its one argument form (with the default argument 1.0 for c), this routine MUST be called with an integer type as the first argument, e.g. DA(1). Expressions involving floating point data types such as DA(1.0) will be interpreted as the constant DA object of the given value instead of the first independent DA variable.
   \param[in] i independent variable number (i=0 means the constant part).
   \param[in] c coefficient corresponding to the given independent variable. By default, this value is assumed to be 1.0.
   \throw DACE::DACEException
 */
    daceAllocateDA(m_index, 0);
    daceCreateVariable(m_index,(unsigned int) i, c);
    if(daceGetError()) DACEException();
}

DA::~DA() throw(){
/*! Destroy a DA object and free the associated object in the DACE core.
 */
    daceFreeDA(m_index);

    // Never throw an exception from a destructor. Instead, we clear the error and ignore it. There is not much the user could do in this case anyway.
    if(daceGetError()) daceClearError();
}

/********************************************************************************
*     Coefficient access routines
*********************************************************************************/
int DA::isnan() const{
/*! Check if a DA object has any NAN coefficients.
   \return True is any coefficients of the DA object are NAN.
   \throw DACE::DACEException
*/
    const int temp = daceIsNan(m_index);
    if(daceGetError()) DACEException();

    return temp;
}

int DA::isinf() const{
/*! Check if a DA object has any INF coefficients.
   \return True is any coefficients of the DA object are INF.
   \throw DACE::DACEException
*/
    const int temp = daceIsInf(m_index);
    if(daceGetError()) DACEException();

    return temp;
}

double DA::cons() const{
/*! Return the constant part of a DA object.
   \return A double corresponding to the constant part of the DA object.
   \throw DACE::DACEException
 */
    const double temp = daceGetConstant(m_index);
    if(daceGetError()) DACEException();

    return temp;
}

AlgebraicVector<double> DA::linear() const{
/*! Return the linear part of a DA object.
   \return An AlgebraicVector<dobule> containing the linear coefficients of
    each independent DA variable in the DA object.
   \throw DACE::DACEException
 */
     AlgebraicVector<double> v(daceGetMaxVariables());

    daceGetLinear(m_index, v.data()); // Note: v.data() is C++11
    if(daceGetError()) DACEException();

    return v;
}

AlgebraicVector<DA> DA::gradient() const {
/*! Compute the gradient of the DA object.
   \return An AlgebraicVector<DA> containing the derivatives
    of the DA object with respect to all independent DA variables.
   \throw DACE::DACEException
 */

    const unsigned int nvar = daceGetMaxVariables();
    AlgebraicVector<DA> temp(nvar);
    for(unsigned int i = 0; i < nvar; i++){
        temp[i] = deriv(i+1);}
    return temp;
}

double DA::getCoefficient(const std::vector<unsigned int> &jj) const {
/*! Return a specific coefficient of a DA object.
   \param[in] jj vector of the exponents of the coefficient to retrieve.
   \return The coefficient of DA object corresponding to the given vector of exponents.
   \throw DACE::DACEException
 */

    double coeff;
    const unsigned int nvar = daceGetMaxVariables();
    if(jj.size() >= nvar)
    {
        coeff = daceGetCoefficient(m_index, jj.data());
    }
    else
    {
        std::vector<unsigned int> temp(jj);
        temp.resize(nvar, 0);
        coeff = daceGetCoefficient(m_index, temp.data());
    }

    if(daceGetError()) DACEException();

    return coeff;
}

void DA::setCoefficient(const std::vector<unsigned int> &jj, const double coeff) {
/*! Set a specific coefficient into a DA object.
   \param[in] jj vector of the exponents of the coefficient to be set.
   \param[in] coeff value to be set as coefficient.
   \throw DACE::DACEException
 */
    // check arguments
    const unsigned int nvar = daceGetMaxVariables();

    if(jj.size() >= nvar)
    {
        daceSetCoefficient(m_index, jj.data(), coeff);
    }
    else
    {
        std::vector<unsigned int> temp(jj);
        temp.resize(nvar, 0);
        daceSetCoefficient(m_index, temp.data(), coeff);
    }

    if(daceGetError()) DACEException();
}

Monomial DA::getMonomial(const unsigned int npos) const{
/*! Return the Monomial corresponding to the non-zero coefficient at the given
    position in the DA object (monomials use one based indexing!).
   \param[in] npos position within the DA object. The ordering of the Monomials
    within a DA object is implementation dependent and does not correspond
    to the order in which Monomials are listed in the ASCII output routines.
   \return A Monomial object containing both the coefficient and
    the vector of exponents corresponding to the given position.
    If the requested monomial is not present in the DA object, a Monomial
    with coefficient set to 0.0 is returned.
   \throw DACE::DACEException
   \sa Monomial
   \sa DA::getMonomial
   \sa DA::getMonomials
*/
    Monomial m;
    getMonomial(npos, m);
    return m;
}


void DA::getMonomial(const unsigned int npos, Monomial &m) const {
/*! Return the Monomial corresponding to the non-zero coefficient at the given
    position in the DA object (monomials use one based indexing!).
   \param[in] npos position within the DA object. The ordering of the Monomials
    within a DA object is implementation dependent and does not correspond
    to the order in which Monomials are listed in the ASCII output routines.
   \param[out] m the monomial object in which to store the corresponding monomial.
   \throw DACE::DACEException
   \sa Monomial
   \sa DA::getMonomial
   \sa DA::getMonomials
*/

    daceGetCoefficientAt(m_index, (int)npos, m.m_jj.data(), m.m_coeff); // Note: m_jj.data() is C++11
    if(daceGetError()) DACEException();
}

std::vector<Monomial> DA::getMonomials() const{
/*! Return a vector of all Monomials in the DA object (differently from
    getMonomial() where only a single Monomial, corresponding to a specified
    position in the DA object, is returned).
   \return A vector of Monomial objects containing both the coefficient and
    the exponents corresponding to each monomial in the DA object. The monomials
    are returned in the same order as in the DACE ASCII output (that is, they
    are sorted by order).
   \throw DACE::DACEException
   \sa Monomial
   \sa DA::getMonomial
 */
    const unsigned int nord = daceGetMaxOrder();
    const unsigned int s = size();
    std::vector<Monomial> res(s), out(s);

    for(unsigned int i = 0; i < s; i++)
    daceGetCoefficientAt(m_index, i+1, res[i].m_jj.data(), res[i].m_coeff); // Note: m_jj.data() is C++11
    if(daceGetError()) DACEException();

    // compute the order of each monomial
    std::vector<unsigned int> sum(s);
    for(unsigned int i = 0; i < s; i++){
        sum[i] = res[i].order();
    }

    // sort monomials by order
    unsigned int k = 0;
    for(unsigned int ord = 0; ord <= nord; ord++){
        for(unsigned int i = 0; i < s; i++){
            if( sum[i] == ord ){
                out[k] = res[i];
                k++;
            }
        }
    }

    return out;
}

/********************************************************************************
*     Assignments
*********************************************************************************/
DA& DA::operator=(DA &&da){
/*! Move all resources from a DA object to us without copying.
   \param[in] da object to be moved.
   \throw DACE::DACEException
 */
    // do a switch-a-roo: we get theirs, they get ours, and then they can even free it for us!
    std::swap(m_index, da.m_index);

    return *this;
}

DA& DA::operator=(const DA &da){
/*! Copy the content of a given DA object da into the DA object.
   \param[in] da DA object to be copied.
   \return The DA object with the same content of the given DA object.
   \throw DACE::DACEException
 */
    if(this != &da){
        daceCopy(da.m_index, m_index);
        if(daceGetError()) DACEException();
    }

    return *this;
}

DA& DA::operator=(const double c){
/*! Copy a constant polynomial of value c into the DA.
   \param[in] c Constant value to be copied.
   \return The DA object representing the constant function with value c.
   \throw DACE::DACEException
 */
    daceCreateConstant(m_index, c);

    if(daceGetError()) DACEException();

    return *this;
}

DA& DA::operator+=(const DA &da){
/*! Compute the sum between the DA object and the given one.
    The result is directly copied into the current DA object.
   \param[in] da DA object to be added.
   \return The current DA object with modified contents.
   \throw DACE::DACEException
 */
    daceAdd(m_index, da.m_index, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

DA& DA::operator+=(const double c){
/*! Compute the sum between the current DA object and a given constant.
    The result is directly copied into the current DA object.
   \param[in] c constant value to be added.
   \return The current DA object with modified contents.
   \throw DACE::DACEException
 */
    daceAddDouble(m_index, c, m_index);

    if(daceGetError()) DACEException();

    return *this;
}

DA& DA::operator-=(const DA &da){
/*! Compute the difference between the current DA object and the given one.
    The result is directly copied into the current DA object.
   \param[in] da DA object to be subtracted.
   \return The current DA object with modified contents.
   \throw DACE::DACEException
 */
    daceSubtract(m_index, da.m_index, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

DA& DA::operator-=(const double c){
/*! Compute the difference between the current DA object and a given constant.
    The result is directly copied into the current DA object.
   \param[in] c constant value to be subtracted.
   \return The current DA object with modified contents.
   \throw DACE::DACEException
 */
    daceSubtractDouble(m_index, c, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

DA& DA::operator*=(const DA &da){
/*! Compute the product between the current DA object and the given one.
    The result is directly copied into the current DA object.
   \param[in] da DA object to be multiplied.
   \return The current DA object with modified contents.
   \throw DACE::DACEException
 */
    daceMultiply(m_index, da.m_index, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

DA& DA::operator*=(const double c){
/*! Compute the product between the current DA object and a given constant.
    The result is directly copied into the current DA object.
   \param[in] c constant value to be multiplied.
   \return The current DA object with modified contents.
   \throw DACE::DACEException
 */
    daceMultiplyDouble(m_index, c, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

DA& DA::operator/=(const DA &da){
/*! Compute the division between the current DA object and the given one.
    The result is directly copied into the current DA object.
   \param[in] da DA object through which the current DA is divided by.
   \return The current DA object with modified contents.
   \throw DACE::DACEException
 */
    daceDivide(m_index, da.m_index, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

DA& DA::operator/=(const double c){
/*! Compute the division between the current DA object and a given constant.
    The result is directly copied into the current DA object.
   \param[in] c constant value through which the current DA is divided by.
   \return The current DA object with modified contents.
   \throw DACE::DACEException
 */
    daceDivideDouble(m_index, c, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

/********************************************************************************
*     Algebraic operations
*********************************************************************************/
DA DA::operator-() const {
/*! Compute the additive inverse of the given DA object.
    The result is copied in a new DA vector.
   \return A new DA object with the opposite sign.
   \throw DACE::DACEException
*/
    DA temp;
    daceMultiplyDouble(m_index,-1.0,temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator+(const DA &da1, const DA &da2){
/*! Compute the addition between two DA objects.
    The result is copied in a new DA object.
   \param[in] da1 first DA object.
   \param[in] da2 second DA object.
   \return A new DA object containing the result of the operation (da1+da2).
   \throw DACE::DACEException
 */
    DA temp;
    daceAdd(da1.m_index, da2.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator + (const DA &da, const double c){
/*! Compute the addition between a DA object and a given constant.
    The result is copied in a new DA object.
   \param[in] da DA object.
   \param[in] c given constant.
   \return A new DA object containing the result of the operation (da+c).
   \throw DACE::DACEException
 */
    DA temp;
    daceAddDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator + (const double c, const DA &da){
/*! Compute the addition between a given constant and a DA object.
    The result is copied in a new DA object.
   \param[in] c given constant.
   \param[in] da DA object.
   \return A new DA object containing the result of the operation (c+da).
   \throw DACE::DACEException
 */
    DA temp;
    daceAddDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator-(const DA &da1, const DA &da2){
/*! Compute the subtraction between two DA objects.
    The result is copied in a new DA object.
   \param[in] da1 first DA object.
   \param[in] da2 second DA object.
   \return A new DA object containing the result of the operation (da1-da2).
   \throw DACE::DACEException
 */
    DA temp;
    daceSubtract(da1.m_index, da2.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator-(const DA &da, const double c){
/*! Compute the subtraction between a DA object and a given constant.
    The result is copied in a new DA object.
   \param[in] da DA object.
   \param[in] c given constant.
   \return A new DA object containing the result of the operation (da-c).
   \throw DACE::DACEException
 */
    DA temp;
    daceSubtractDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator-(const double c, const DA &da){
/*! Compute the subtraction between a given constant and a DA object.
    The result is copied in a new DA object.
   \param[in] c given constant.
   \param[in] da DA object.
   \return A new DA object containing the result of the operation (c-da).
   \throw DACE::DACEException
 */
    DA temp;
    daceDoubleSubtract(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator*(const DA &da1, const DA &da2){
/*! Compute the multiplication between two DA objects.
    The result is copied in a new DA object.
   \param[in] da1 first DA object.
   \param[in] da2 second DA object.
   \return A new DA object containing the result of the operation (da1*da2).
   \throw DACE::DACEException
 */
    DA temp;
    daceMultiply(da1.m_index, da2.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator*(const DA &da, const double c){
/*! Compute the multiplication between a DA object and a given constant.
    The result is copied in a new DA object.
   \param[in] da DA object.
   \param[in] c given constant.
   \return A new DA object containing the result of the operation (da*c).
   \throw DACE::DACEException
 */
    DA temp;
    daceMultiplyDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator*(const double c, const DA &da){
/*! Compute the multiplication between a given constant and a DA object.
    The result is copied in a new DA object.
   \param[in] c given constant.
   \param[in] da DA object.
   \return A new DA object containing the result of the operation (c*da).
   \throw DACE::DACEException
 */
    DA temp;
    daceMultiplyDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator/(const DA &da1, const DA &da2){
/*! Compute the division between two DA objects.
    The result is copied in a new DA object.
   \param[in] da1 first DA object.
   \param[in] da2 second DA object.
   \return A new DA object containing the result of the operation (da1/da2).
   \throw DACE::DACEException
 */
    DA temp;
    daceDivide(da1.m_index, da2.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator/(const DA &da, const double c){
/*! Compute the division between a DA object and a given constant.
    The result is copied in a new DA object.
   \param[in] da DA object.
   \param[in] c given constant.
   \return A new DA object containing the result of the operation (da/c).
   \throw DACE::DACEException
 */
    DA temp;
    daceDivideDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA operator/(const double c, const DA &da){
/*! Compute the division between a given constant and a DA object.
    The result is copied in a new DA object.
   \param[in] c given constant.
   \param[in] da DA object.
   \return A new DA object containing the result of the operation (c/da).
   \throw DACE::DACEException
 */
    DA temp;
    daceDoubleDivide(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/********************************************************************************
*     Math routines
*********************************************************************************/

DA DA::multiplyMonomials(const DA &da) const {
/*! Multiply the DA vector with another DA vector da monomial by monomial.
    This is the equivalent of coefficient-wise multiplication (like in DA addition).
   \param[in] da DA vector to multiply with coefficient-wise
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::evalMonomial
*/
    DA temp;
    daceMultiplyMonomials(m_index, da.m_index, temp.m_index);
    if (daceGetError()) DACEException();

    return temp;
}

DA DA::divide(const unsigned int var, const unsigned int p) const{
/*! Divide by independent variable var raised to power p.
    The result is copied in a new DA object.
   \param[in] var independente variable number to divide by.
   \param[in] p power of the independent variable.
   \return A new DA object containing the result.
   \throw DACE::DACEException
 */
    DA temp;
    daceDivideByVariable(m_index, var, p, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::deriv(const unsigned int i) const{
/*! Compute the derivative of a DA object with respect to variable i.
    The result is copied in a new DA object.
   \param[in] i variable with respect to which the derivative is calculated.
   \return A new DA object containing the result of the derivation.
   \throw DACE::DACEException
 */
    DA temp;
    daceDifferentiate(i, m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::deriv(const std::vector<unsigned int> ind) const{
/*! Compute the derivative of a DA object with respect to variables ind.
    The result is copied in a new DA object.
   \param[in] ind vector containing the number of derivatives to take for each
    independent variable. If ind has fewer entries than there are independent
    variables, the missing entries are assumed to be zero. If ind has more
    entries than there are independent variables, extra values are ignored.
   \return A new DA object containing the result of the derivation.
   \throw DACE::DACEException
 */
    DA temp(*this);
    const unsigned int size = std::min((unsigned int)ind.size(),daceGetMaxVariables());

    for(unsigned int i=0; i<size; i++)
        for(unsigned int j=0; j<ind[i]; j++)
            daceDifferentiate((i+1), temp.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::integ(const unsigned int i) const{
/*! Compute the integral of a DA object with respect to variable i.
    The result is copied in a new DA object.
   \param[in] i variable with respect to which the integral is calculated.
   \return A new DA object containing the result of the integration.
   \throw DACE::DACEException
 */
    DA temp;
    daceIntegrate(i, m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::integ(const std::vector<unsigned int> ind) const{
/*! Compute the integral of a DA object with respect to variables i.
    The result is copied in a new DA object.
   \param[in] ind vector containing the number of integrals to take for each
    independent variable. If ind has fewer entries than there are independent
    variables, the missing entries are assumed to be zero. If ind has more
    than there are independent variables, extra values are ignored.
   \return A new DA object containing the result of the integration.
   \throw DACE::DACEException
 */
    DA temp(*this);
    const unsigned int size = std::min((unsigned int)ind.size(),daceGetMaxVariables());

    for(unsigned int i=0; i<size; i++)
        for(unsigned int j=0; j<ind[i]; j++)
            daceIntegrate((i+1), temp.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::trim(const unsigned int min, const unsigned int max) const{
/*! Returns a DA object with all monomials of order less than min and greater
    than max removed.
    The result is copied in a new DA object.
   \param[in] min The minimum order to keep in the DA object.
   \param[in] max The maximum order to keep in the DA object.
   \return A new DA object containing the result of the trimming.
   \throw DACE::DACEException
 */
    DA temp;
    daceTrim(m_index, min, max, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::trunc() const{
/*! Truncate the constant part of a DA object to an integer.
    The result is copied in a new DA object.
   \return A new DA object with a truncated constant part.
   \throw DACE::DACEException
 */
    DA temp;
    daceTruncate(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::round() const{
/*! Round the constant part of a DA object to an integer.
    The result is copied in a new DA object.
   \return A new DA object with a rounded constant part.
   \throw DACE::DACEException
 */
    DA temp;
    daceRound(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::mod(const double p) const{
/*! Compute the floating-point remainder of c/p (c modulo p),
    where c is the constant part of the current DA object.
    The result is copied in a new DA object.
   \param[in] p costant with respect to which the modulo function is computed.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceModulo(m_index, p, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::pow(const int p) const{
/*! Elevate a DA object to a given integer power.
    The result is copied in a new DA object.
   \param[in] p power to which the DA object is raised.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    dacePower(m_index, p, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::pow(const double p) const{
/*! Elevate a DA object to a given real power. The constant part must be positive.
    The result is copied in a new DA object.
   \param[in] p power to which the DA object is raised.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    dacePowerDouble(m_index, p, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::root(const int p) const{
/*! Compute the p-th root of a DA object.
    The result is copied in a new DA object.
   \param[in] p root to be computed.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceRoot(m_index, p, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::minv() const{
/*! Compute the multiplicative inverse of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceMultiplicativeInverse(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::sqr() const{
/*! Compute the square of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceSquare(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::sqrt() const{
/*! Compute the square root of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceSquareRoot(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::isrt() const{
/*! Compute the inverse square root of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceInverseSquareRoot(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::cbrt() const{
/*! Compute the cubic root of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceCubicRoot(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::icrt() const{
/*! Compute the inverse cubic root of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceInverseCubicRoot(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::hypot(const DA &da) const{
/*! Compute the hypotenuse (sqrt(a*a+b*b)) of a DA object and the given DA argument.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceHypotenuse(m_index, da.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::exp() const{
/*! Compute the exponential of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceExponential(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::log() const{
/*! Compute the natural logarithm of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceLogarithm(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::logb(const double b) const{
/*! Compute the logarithm of a DA object with respect to a given base.
    The result is copied in a new DA object.
   \param[in] b base with respect to which the logarithm is computed (base 10 set as default base).
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceLogarithmBase(m_index, b, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::log10() const{
/*! Compute the 10 based logarithm of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceLogarithm10(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::log2() const{
/*! Compute the 2 based logarithm of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceLogarithm2(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::sin() const{
/*! Compute the sine of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceSine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::cos() const{
/*! Compute the cosine of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceCosine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::tan() const{
/*! Compute the tangent of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceTangent(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::asin() const{
/*! Compute the arcsine of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceArcSine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::acos() const{
/*! Compute the arccosine of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceArcCosine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::atan() const{
/*! Compute the arctangent of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceArcTangent(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::atan2(const DA &da) const{
/*! Compute the four-quadrant arctangent of Y/X. Y is the current DA object,
    whereas X is the given da.
    The result is copied in a new DA object.
   \param[in] da DA object
   \return A new DA object containing the result of the operation Y/X in [-pi, pi].
   \throw DACE::DACEException
 */
    DA temp;
    daceArcTangent2(m_index, da.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::sinh() const{
/*! Compute the hyperbolic sine of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceHyperbolicSine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::cosh() const{
/*! Compute the hyperbolic cosine of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceHyperbolicCosine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::tanh() const{
/*! Compute the hyperbolic tangent of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceHyperbolicTangent(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::asinh() const{
/*! Compute the hyperbolic arcsine of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceHyperbolicArcSine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::acosh() const{
/*! Compute the hyperbolic arccosine of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceHyperbolicArcCosine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::atanh() const{
/*! Compute the hyperbolic arctangent of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceHyperbolicArcTangent(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}


DA DA::erf() const{
/*! Compute the error function of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceErrorFunction(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::erfc() const{
/*! Compute the complementary error function of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceComplementaryErrorFunction(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::BesselJFunction(const int n) const{
/*! Compute the n-th Bessel function of first type J_n of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Bessel function
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \note The DA must have non-negative constant part while the order is allowed to be negative.
   \note This function fails if the result is too large to be represented in double precision.
 */
    DA temp;
    daceBesselJFunction(m_index, n, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::BesselYFunction(const int n) const{
/*! Compute the n-th Bessel function of second type Y_n of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Bessel function
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \note The DA must have non-negative constant part while the order is allowed to be negative.
   \note This function fails if the result is too large to be represented in double precision.
 */
    DA temp;
    daceBesselYFunction(m_index, n, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::BesselIFunction(const int n, const bool scaled) const{
/*! Compute the n-th modified Bessel function of first type I_n of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Bessel function
   \param[in] scaled if true, the modified Bessel function is scaled
    by a factor exp(-x), i.e. exp(-x)I_n(x) is returned.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \note The DA must have non-negative constant part while the order is allowed to be negative.
   \note This function fails if the result is too large to be represented in double precision.
 */
    DA temp;
    daceBesselIFunction(m_index, n, scaled, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::BesselKFunction(const int n, const bool scaled) const{
/*! Compute the n-th modified Bessel function of second type K_n of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Bessel function
   \param[in] scaled if true, the modified Bessel function is scaled
    by a factor exp(x), i.e. exp(x)K_n(x) is returned.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \note The DA must have non-negative constant part while the order is allowed to be negative.
   \note This function fails if the result is too large to be represented in double precision.
 */
    DA temp;
    daceBesselKFunction(m_index, n, scaled, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::GammaFunction() const{
/*! Compute the Gamma function of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceGammaFunction(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::LogGammaFunction() const{
/*! Compute the Logarithmic Gamma function (i.e. the natural logarithm of Gamma) of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    daceLogGammaFunction(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::PsiFunction(const unsigned int n) const{
/*! Compute the n-th order Psi function, i.e. the (n+1)st derivative of the Logarithmic Gamma function, of a DA object.
    The result is copied in a new DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
 */
    DA temp;
    dacePsiFunction(m_index, n, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/********************************************************************************
*    Norm and estimation routines
*********************************************************************************/
unsigned int DA::size() const{
/*! Return the number of non-zero coefficients of a DA object.
   \return The number of non-zero coefficients stored in the DA object.
   \throw DACE::DACEException
 */
    unsigned int res;
    res=daceGetLength(m_index);
    if(daceGetError()) DACEException();

    return res;
}

double DA::abs() const{
/*! Compute the max norm of a DA object.
   \return A double corresponding to the result of the operation.
   \throw DACE::DACEException
 */
    double c;
    c=daceAbsoluteValue(m_index);
    if(daceGetError()) DACEException();

    return c;
}

double DA::norm(const unsigned int type) const{
/*! Compute different types of norms for a DA object.
   \param[in] type type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
   \return A double corresponding to the result of the operation.
   \throw DACE::DACEException
 */
    double c;
    c=daceNorm(m_index, type);
    if(daceGetError()) DACEException();

    return c;
}

std::vector<double> DA::orderNorm(const unsigned int var, const unsigned int type) const{
/*! Extract different types of order sorted norms from a DA object.
   \param[in] var order\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable var
   \param[in] type type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
   \return A double corresponding to the result of the operation.
   \throw DACE::DACEException
 */
    std::vector<double> v(daceGetMaxOrder()+1);
    daceOrderedNorm(m_index, var, type, v.data()); // Note: v.data() is C++11
    if(daceGetError()) DACEException();

    return v;
}

std::vector<double> DA::estimNorm(const unsigned int var, const unsigned int type, const unsigned int nc) const{
/*! Estimate different types of order sorted norms for terms of a DA object
    up to a specified order.
   \param[in] var order\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable var
   \param[in] type type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
   \param[in] nc maximum order to be estimated (Default order = Max order)
   \return A double corresponding to the result of the operation.
   \throw DACE::DACEException
   \note If estimation is not possible, zero is returned for all requested orders.
 */
    std::vector<double> v(nc+1);
    daceEstimate(m_index, var, type, v.data(), NULL ,nc); // Note: v.data() is C++11
    if(daceGetError()) DACEException();

    return v;
}

std::vector<double> DA::estimNorm(std::vector<double> &err, const unsigned int var, const unsigned int type, const unsigned int nc) const{
/*! Estimate different types of order sorted norms for terms of a DA object
    up to a specified order with error estimates.
   \param[out] err returns the amount by which the estimate underestimates the actual ordered norm of the terms in the polynomial up to the minimum of nc or the maximum computation order.
   \param[in] var order\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable var
   \param[in] type type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
   \param[in] nc maximum order to be estimated (Default order = Max order)
   \return A double corresponding to the result of the operation.
   \throw DACE::DACEException
   \note If estimation is not possible, zero is returned for all requested orders.
 */
    std::vector<double> v(nc+1);
    err.resize(std::min(nc, daceGetMaxOrder())+1);
    daceEstimate(m_index, var, type, v.data(), err.data(), nc);
    if(daceGetError()) DACEException();

    return v;
}

Interval DA::bound() const{
/*! Compute lower and upper bounds of a DA object.
   \return An Interval object containing both the lower and the upper bound
    of the DA object.
   \throw DACE::DACEException
   \sa Interval
 */
    Interval i;
    daceGetBounds(m_index, i.m_lb, i.m_ub);
    if(daceGetError()) DACEException();

    return i;
}

double DA::convRadius(const double eps, const unsigned int type) const{
/*! Estimate the convergence radius of the DA object.
   \param[in] eps requested tolerance.
   \param[in] type type of norm (sum norm is used as default)
   \return A double corresponding to the estimated convergence radius.
   \throw DACE::DACEException
 */
    const unsigned int ord = daceGetTruncationOrder();

    std::vector<double> res = estimNorm(0, type, ord+1);
    return std::pow(eps/res[ord+1],1.0/(ord+1));
}

/********************************************************************************
*     DACE polynomial evaluation routines
*********************************************************************************/
compiledDA DA::compile() const{
/*! Compile current DA object and create a compiledDA object.
   \return The compiled DA object.
   \throw DACE::DACEException
 */
    return compiledDA(*this);
}

DA DA::plug(const unsigned int var, const double val) const{
/*! Partial evaluation of a DA object. In the DA object, variable var is
    replaced by the value val. The resulting DA object is returned.
   \param[in] var variable number to be replaced
   \param[in] val value by which to replace the variable
   \return A new DA object containing the resulting DA object.
   \throw DACE::DACEException
 */
    DA temp;
    daceEvalVariable(m_index,var,val,temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

double DA::evalMonomials(const DA &values) const {
/*! Evaluates the DA vector using the coefficients in argument values as the values for each monomial.
    This is equivalent to a monomial-wise dot product of two DA vectors.
   \param[in] values DA vector containing the values of each monomial
   \return The result of the evaluation.
   \throw DACE::DACEException
   \sa DA::multiplyMonomial
*/
    const double res = daceEvalMonomials(m_index, values.m_index);
    if (daceGetError()) DACEException();

    return res;
}

DA DA::replaceVariable(const unsigned int from, const unsigned int to, const double val) const{
/*! Partial evaluation of a DA object. In the DA object, variable from is
    replaced by the value val times variable to. The resulting DA object is returned.
   \param[in] from variable number to be replaced
   \param[in] to variable number to be inserted instead
   \param[in] val value by which to scale the inserted variable
   \return A new DA object containing the resulting DA object.
   \throw DACE::DACEException
   \sa DA::replaceVariable
 */
    DA temp;
    daceReplaceVariable(m_index, from, to, val, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::scaleVariable(const unsigned int var, const double val) const{
/*! Scaling of an independent variable. In the DA object, variable var is
    replaced by the value val times var. The resulting DA object is returned.
   \param[in] var variable number to be scaled
   \param[in] val value by which to scale the variable
   \return A new DA object containing the resulting DA object.
   \throw DACE::DACEException
   \sa DA::scaleVariable
 */
    DA temp;
    daceScaleVariable(m_index, var, val, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::translateVariable(const unsigned int var, const double a, const double c) const{
/*! Affine translation of an independent variable. In the DA object, variable var is
    replaced by a*var + c. The resulting DA object is returned.
   \param[in] var variable number to be translated
   \param[in] a value by which to scale the variable
   \param[in] c value by which to shift the variable
   \return A new DA object containing the resulting DA object.
   \throw DACE::DACEException
   \sa DA::translateVariable
 */
    DA temp;
    daceTranslateVariable(m_index, var, a, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/********************************************************************************
*     DACE input/output routines
*********************************************************************************/
std::string DA::toString() const{
/*!  Convert DA object to string.
    \return A string.
   \throw DACE::DACEException
 */
    // initialize 2D char array
    unsigned int nstr = daceGetMaxMonomials() + 2;
    char *ss = new char[nstr*DACE_STRLEN];

    // call dacewrite
    daceWrite( m_index, ss, nstr );

    // copy from char array to string
    std::string s;
    for(unsigned int i=0; i<(unsigned int)nstr; i++){
        ss[(i+1)*DACE_STRLEN-1] = '\0'; // should already be done by the Fortran code but just to be sure terminate string
        s.append(&ss[i*DACE_STRLEN]);
        s.append(1,'\n');
    }

    // delete 2D char array
    delete[] ss;

    if(daceGetError()) DACEException();

    return s;
}

std::ostream& operator<<(std::ostream &out, const DA &da){
/*! Overload of std::operator<< in iostream.
   \param[in] out standard output stream.
   \param[in] da DA object to be printed in the stream
   \return Standard output stream.
   \throw DACE::DACEException
   \sa DA::toString
 */
    out << toString(da);

    return out;
}


std::istream& operator>>(std::istream &in, DA &da){
/*! Overload of std::operator>> in iostream. Reads both string and binary
    DA representations from a file.
   \param[in] in standard input stream.
   \param[in] da DA object to be created from the stream
   \return Standard input stream.
   \throw DACE::DACEException
   \note When using binary IO operations, make sure the stream is opened in ios_base::binary mode!
    Some C++ libraries are known to mangle the input otherwise which will break the ability to read binary DA objects.
    Setting the binary flag for all IO (also text based) does not affect the output and is recommended.
   \sa DA::fromString
 */
    storedDA sda(in);

    if(sda.isValid())
    {
        da = sda;                                       // automatically cast sDA to DA
        return in;
    }
    else
    {
        const std::string endstr = "------------------------------------------------";  // from daceio.c
        std::string line = sda;                         // automatically cast sDA to string
        std::vector<std::string> strs;
        strs.reserve(5000);                             // estimate that 5000 lines will be an OK guess for many use cases

        if(line.length() > 0)
        {
            // parse the content of line into string array
            const std::string::iterator end = line.end();
            std::string::iterator p = line.begin();
            while(p != end)
            {
                const std::string::iterator p0 = p;
                while(p != end && *p != '\n') p++;
                strs.emplace_back(p0, p);
                if(p != end) p++;
            }

            // complete the last line, in case it was not a full line that was read
            if(*(end-1) != '\n')
            {
                getline(in, line);
                strs[strs.size()-1] += line;
            }
        }

        if (!strs.empty())
        {
            if(!strs.back().empty())
            {
                // check that last line is not the terminator line and in case remove it
                if(strs.back().compare(4, 31, endstr, 0, 31) == 0)
                {
                    strs.pop_back();

                }
                else
                {
                    // read the istream until end string is found and put each line in the string vector (end condition taken from daceio.c)
                    for(getline(in, line); in.good() && (line.compare(4, 31, endstr, 0, 31) != 0); getline(in, line))
                    strs.push_back(line);
                }
            }
            // convert string vector to DA
            da = DA::fromString(strs);
        }
    }

    return in;
}

void DA::write(std::ostream &os) const{
/*!  Write a binary representation of the DA as a blob to os.
    \throw DACE::DACEException
 */
    os << storedDA(*this);
}


/********************************************************************************
*     DACE static factory routines
*********************************************************************************/
DA DA::random(const double cm){
/*! Create a DA object and fill with random entries.
   \param[in] cm filling factor (for cm < 0, the DA object is filled with random numbers;
    for cm > 0, the DA object is filled with weighted decaying numbers).
   \throw DACE::DACEException
 */
    DA temp;
    daceCreateRandom(temp.m_index, cm);
    if(daceGetError()) DACEException();

    return temp;
}

DA DA::identity(const unsigned int var){
/*! Create a DA object representing the identity function in independent DA
    variable number var. This is just an alias for DA::DA(var).
   \param[in] var The independent DA variable number
   \throw DACE::DACEException
   \sa DA::DA
 */
    return DA((int)var);
}

DA DA::fromString(const std::string &str){
/*! Convert a string to DA object.
   \param[in] str string.
   \return A DA object.
   \throw DACE::DACEException
 */
    std::istringstream ssin(str);
    DA temp;
    ssin >> temp;

    return temp;
}

DA DA::fromString(const std::vector<std::string> &str){
/*! Convert a vector of strings to DA object.
   \param[in] str vector of strings, each representing one line of the input.
   \return A DA object.
   \throw DACE::DACEException
 */
    // create 2D char array
    const unsigned int nstr = (unsigned int)str.size();
    char *ss = new char[nstr*DACE_STRLEN];

    // fill the array with blanks
    for(unsigned int i=0; i<nstr*DACE_STRLEN; i++)
        ss[i] = ' ';

    // fill char array with rows of the string vector
    for(unsigned int i=0; i<nstr; i++){
        str[i].copy(&ss[i*DACE_STRLEN], DACE_STRLEN);}

    // call to daceread
    DA da;
    daceRead(da.m_index, ss, nstr);

    // delete 2D char array
    delete[] ss;

    if(daceGetError()) DACEException();

    return da;
}

DA DA::read(std::istream &is){
/*!  Read a binary representation of a DA from is.
    \throw DACE::DACEException
    \note When using binary IO operations, make sure the stream is opened in ios_base::binary mode!
     Some C++ libraries are known to mangle the input otherwise which will break the ability to read binary DA objects.
     Setting the binary flag for all IO (also text based) does not affect the output and is recommended.
 */
    storedDA sda(is);
    return sda;         // automatically cast to DA
}

/********************************************************************************
*     DACE various routines
*********************************************************************************/
void DA::memdump(){
    daceMemoryDump();
}

/********************************************************************************
*     DACE non-member functions
*********************************************************************************/
int isnan(const DA &da) {
/*! Check if a DA object has any NAN coefficients.
   \param[in] da a given DA object.
   \return True if any coefficients of the DA object are NAN.
   \throw DACE::DACEException
*/
    return da.isnan();}

int isinf(const DA &da) {
/*! Check if a DA object has any INF coefficients.
   \param[in] da a given DA object.
   \return True if any coefficients of the DA object are INF.
   \throw DACE::DACEException
*/
    return da.isinf();}

double cons(const DA &da) {
/*! Return the constant part of a DA object.
   \param[in] da a given DA object.
   \return A double corresponding to the constant part of the DA object.
   \throw DACE::DACEException
 */
    return da.cons();}

AlgebraicVector<double> linear(const DA &da) {
/*! Return the linear part of a DA object.
   \param[in] da a given DA object.
   \return An AlgebraicVector<dobule> containing the linear coefficients of
    each independent DA variable in the given DA object.
   \throw DACE::DACEException
 */
    return da.linear();}

AlgebraicVector<DA> gradient(const DA &da) {
/*! Compute the gradient of a DA object.
   \param[in] da the given DA object.
   \return A AlgebraicVector<DA> containing the derivatives
    of the DA object with respect to all independent DA variables.
   \throw DACE::DACEException
 */

    return da.gradient();}

DA divide(const DA &da, const unsigned int var, const unsigned int p){
/*! Divide by independent variable var raised to power p.
    The result is copied in a new DA object.
   \param[in] da DA object.
   \param[in] var variable number to divide by.
   \param[in] p power of variable var to divide by.
   \return A new DA object containing the result of the division.
   \throw DACE::DACEException
   \sa DA::divide
 */
    return da.divide(var, p);}


DA deriv(const DA &da, const unsigned int i){
/*! Compute the derivative of a DA object with respect to variable i.
    The result is copied in a new DA object.
   \param[in] da DA object.
   \param[in] i variable with respect to which the derivative is calculated.
   \return A new DA object containing the result of the derivation.
   \throw DACE::DACEException
   \sa DA::deriv
 */
    return da.deriv(i);}

DA deriv(const DA &da, const std::vector<unsigned int> ind){
/*! Compute the derivative of a DA object with respect to variables ind.
    The result is copied in a new DA object.
   \param[in] da DA object.
   \param[in] ind vector containing the number of derivatives to take for each
    independent variable. If ind has fewer entries than there are independent
    variables, the missing entries are assumed to be zero. If ind has more
    entries than there are independent variables, extra values are ignored.
   \return A new DA object containing the result of the derivation.
   \throw DACE::DACEException
   \sa DA::deriv
 */
    return da.deriv(ind);}

DA integ(const DA &da, const unsigned int i){
/*! Compute the integral of a DA object with respect to variable i.
    The result is copied in a new DA object.
   \param[in] da DA object.
   \param[in] i variable with respect to which the integral is calculated.
   \return A new DA object containing the result of the integration.
   \throw DACE::DACEException
   \sa DA::integ
 */
    return da.integ(i);}

DA integ(const DA &da, const std::vector<unsigned int> ind){
/*! Compute the integral of a DA object with respect to variable i.
    The result is copied in a new DA object.
   \param[in] da DA object.
   \param[in] ind vector containing the number of derivatives to take for each
    independent variable. If ind has fewer entries than there are independent
    variables, the missing entries are assumed to be zero. If ind has more
    entries than there are independent variables, extra values are ignored.
   \return A new DA object containing the result of the integration.
   \throw DACE::DACEException
   \sa DA::integ
 */
    return da.integ(ind);}

DA trim(const DA &da, const unsigned int min, const unsigned int max){
/*! Returns a DA object with all monomials of order less than min and greater
    than max removed (trimmed).
    The result is copied in a new DA object.
   \param[in] da DA object.
   \param[in] min The minimum order to keep in the DA object.
   \param[in] max The maximum order to keep in the DA object.
   \return A new DA object containing the result of the trimming.
   \throw DACE::DACEException
   \sa DA::trim
 */
    return da.trim(min,max);}

DA trunc(const DA &da){
/*! Truncate the constant part of a DA object to an integer.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object with a truncated constant part.
   \throw DACE::DACEException
   \sa DA::trunc
 */
    return da.trunc();}

DA round(const DA &da){
/*! Round the constant part of a DA object to an integer.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object with a rounded constant part.
   \throw DACE::DACEException
   \sa DA::round
 */
    return da.round();}

DA mod(const DA &da, double p){
/*! Compute the floating-point remainder of c/p (c modulo p),
    where c is the constant part of the given DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \param[in] p costant with respect to which the modulo function is computed.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::mod
 */
    return da.mod(p);}

DA pow(const DA &da, int p){
/*! Raise a DA object to a given integer power.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \param[in] p power at which the DA object is elevated.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::pow
 */
    return da.pow(p);}


DA pow(const DA &da, double p){
/*! Raise a DA object to a given integer power. The constant part must be positive.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \param[in] p power at which the DA object is elevated.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::pow
 */
    return da.pow(p);}


DA root(const DA &da, int p){
/*! Compute the p-th root of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \param[in] p root to be computed.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::root
 */
    return da.root(p);}

DA minv(const DA &da){
/*! Compute the multiplicative inverse of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::minv
 */
    return da.minv();}

DA sqr(const DA &da){
/*! Compute the square of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::sqr
 */
    return da.sqr();}

DA sqrt(const DA &da){
/*! Compute the square root of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::sqrt
 */
    return da.sqrt();}

DA isrt(const DA &da){
/*! Compute the inverse square root of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::isrt
 */
    return da.isrt();}


DA cbrt(const DA &da){
/*! Compute the cubic root of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::cbrt
 */
    return da.cbrt();}

DA icrt(const DA &da){
/*! Compute the inverse cubic root of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::icrt
 */
    return da.icrt();}

DA hypot(const DA &da1, const DA &da2){
/*! Compute the hypotenuse (sqrt(da1*da1 + da2*da2)) of two DA objects.
    The result is copied in a new DA object.
   \param[in] da1 first DA object.
   \param[in] da2 second DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::hypot
 */
    return da1.hypot(da2);}


DA exp(const DA &da){
/*! Compute the exponential of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::exp
 */
    return da.exp();}

DA log(const DA &da){
/*! Compute the natural logarithm of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::log
 */
    return da.log();}

DA logb(const DA &da, const double b){
/*! Compute the logarithm of a DA object with respect to a given base.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \param[in] b base with respect to which the logarithm is computed (base 10 set as default base).
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::logb
 */
    return da.logb(b);}


DA log10(const DA &da){
/*! Compute the 10 based logarithm of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::log10
 */
    return da.log10();}

DA log2(const DA &da){
/*! Compute the 2 based logarithm of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::log2
 */
    return da.log2();}


DA sin(const DA &da){
/*! Compute the sine of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::sin
 */
    return da.sin();}

DA cos(const DA &da){
/*! Compute the cosine of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::cos
 */
    return da.cos();}

DA tan(const DA &da){
/*! Compute the tangent of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::tan
 */
    return da.tan();}

DA asin(const DA &da){
/*! Compute the arcsine of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::asin
 */
    return da.asin();}

DA acos(const DA &da){
/*! Compute the arccosine of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::acos
 */
    return da.acos();}

DA atan(const DA &da){
/*! Compute the arctangent of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::atan
 */
    return da.atan();}

DA atan2(const DA &da1, const DA &da2){
/*! Compute the four-quadrant arctangent of da1/da2.
    The result is copied in a new DA object.
   \param[in] da1 DA object
   \param[in] da2 DA object
   \return A new DA object containing the result of the operation in [-pi, pi].
   \throw DACE::DACEException
   \sa DA::atan2
 */
    return da1.atan2(da2);}

DA sinh(const DA &da){
/*! Compute the hyperbolic sine of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::sinh
 */
    return da.sinh();}

DA cosh(const DA &da){
/*! Compute the hyperbolic cosine of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::cosh
 */
    return da.cosh();}

DA tanh(const DA &da){
/*! Compute the hyperbolic tangent of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::tanh
 */
    return da.tanh();}

DA asinh(const DA &da){
/*! Compute the hyperbolic arcsine of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::asinh
 */
    return da.asinh();}

DA acosh(const DA &da){
/*! Compute the hyperbolic arccosine of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::acosh
 */
    return da.acosh();}

DA atanh(const DA &da){
/*! Compute the hyperbolic arctangent of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::atanh
 */
    return da.atanh();}


DA erf(const DA &da){
/*! Compute the error function of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::erf
 */
    return da.erf();}

DA erfc(const DA &da){
/*! Compute the complementary error function of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::erfc
 */
    return da.erfc();}

DA jn(const int n, const DA &da){
/*! Compute the n-th Bessel function of first type J_n of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Bessel function.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \note The DA must have non-negative constant part while the order is allowed to be negative.
   \note Alias of BesselJFunction for C compatible naming.
   \sa DA::BesselJFunction
 */
    return da.BesselJFunction(n);}

DA yn(const int n, const DA &da){
/*! Compute the n-th Bessel function of second type Y_n of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Bessel function.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \note The DA must have non-negative constant part while the order is allowed to be negative.
   \note Alias of BesselYFunction for C compatible naming.
   \sa DA::BesselYFunction
 */
    return da.BesselYFunction(n);}

DA BesselJFunction(const int n, const DA &da){
/*! Compute the n-th Bessel function of first type J_n of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Bessel function
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \note The DA must have non-negative constant part while the order is allowed to be negative.
    \sa DA::BesselJFunction
*/
    return da.BesselJFunction(n);}

DA BesselYFunction(const int n, const DA &da){
/*! Compute the n-th Bessel function of second type Y_n of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Bessel function
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \note The DA must have non-negative constant part while the order is allowed to be negative.
   \sa DA::BesselYFunction
 */
    return da.BesselYFunction(n);}

DA BesselIFunction(const int n, const DA &da, const bool scaled){
/*! Compute the n-th modified Bessel function of first type I_n of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Bessel function
   \param[in] da a given DA object.
   \param[in] scaled if true, the modified Bessel function is scaled
    by a factor exp(-x), i.e. exp(-x)I_n(x) is returned.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \note The DA must have non-negative constant part while the order is allowed to be negative.
   \sa DA::BesselIFunction
 */
    return da.BesselIFunction(n, scaled);}

DA BesselKFunction(const int n, const DA &da, const bool scaled){
/*! Compute the n-th modified Bessel function of second type K_n of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Bessel function
   \param[in] da a given DA object.
   \param[in] scaled if true, the modified Bessel function is scaled
    by a factor exp(-x), i.e. exp(-x)K_n(x) is returned.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \note The DA must have non-negative constant part while the order is allowed to be negative.
   \sa DA::BesselKFunction
 */
    return da.BesselKFunction(n, scaled);}

DA tgamma(const DA &da){
/*! Compute the Gamma function of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::GammaFunction
   \note Alias of GammaFunction() for C99 compatible naming.
 */
    return da.GammaFunction();}

DA lgamma(const DA &da){
/*! Compute the Log Gamma function of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::LogGammaFunction
   \note Alias of LogGammaFunction() for C99 compatible naming.
 */
    return da.LogGammaFunction();}

DA GammaFunction(const DA &da){
/*! Compute the Gamma function of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::GammaFunction
 */
    return da.GammaFunction();}

DA LogGammaFunction(const DA &da){
/*! Compute the Log Gamma function of a DA object.
    The result is copied in a new DA object.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::LogGammaFunction
 */
    return da.LogGammaFunction();}

DA PsiFunction(const unsigned int n, const DA &da){
/*! Compute the n-th Psi function of a DA object.
    The result is copied in a new DA object.
   \param[in] n order of the Psi function to compute.
   \param[in] da a given DA object.
   \return A new DA object containing the result of the operation.
   \throw DACE::DACEException
   \sa DA::PsiFunction
 */
    return da.PsiFunction(n);}

unsigned int size(const DA &da){
/*! Return the number of non-zero coefficients of a DA object.
   \param[in] da a given DA object.
   \return The number of non-zero coefficients of the DA object.
   \throw DACE::DACEException
   \sa DA::size
 */
    return da.size();}

double abs(const DA &da){
/*! Compute the max norm of a DA object.
   \param[in] da a given DA object.
   \return A double corresponding to the result of the operation.
   \throw DACE::DACEException
   \sa DA::abs
 */
    return da.abs();}

double norm(const DA &da, unsigned int type){
/*! Compute different types of norms for a DA object.
   \param[in] da a given DA object.
   \param[in] type type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
   \return A double corresponding to the result of the operation.
   \throw DACE::DACEException
   \sa DA::norm
 */
    return da.norm(type);}

std::vector<double> orderNorm(const DA &da, const unsigned int var, const unsigned int type){
/*! Compute different types of order sorted norms for terms of a DA object.
   \param[in] da a given DA object.
   \param[in] var order\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable var
   \param[in] type type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
   \return A double corresponding to the result of the operation.
   \throw DACE::DACEException
   \sa DA::onorm
 */
    return da.orderNorm(var, type);}

std::vector<double> estimNorm(const DA &da, const unsigned int var, const unsigned int type, const unsigned int nc){
/*! Estimate different types of order sorted norms for terms of a DA object
    up to a specified order.
   \param[in] da a given DA object.
   \param[in] var order\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable var
   \param[in] type type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
   \param[in] nc maximum order (Default order = Max order)
   \return A double corresponding to the result of the operation.
   \throw DACE::DACEException
   \note If estimation is not possible, zero is returned for all requested orders.
   \sa DA::estim
 */
    return da.estimNorm(var, type, nc);}


std::vector<double> estimNorm(const DA &da, std::vector<double> &err, const unsigned int var, const unsigned int type, const unsigned int nc){
/*! Estimate different types of order sorted norms for terms of a DA object
    up to a specified order with error estimates.
   \param[in] da a given DA object.
   \param[out] err returns the amount by which the estimate underestimates the actual ordered norm of the terms in the polynomial up to the minimum of nc or the maximum computation order.
   \param[in] var order\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable var
   \param[in] type type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
   \param[in] nc maximum order (Default order = Max order)
   \return A double corresponding to the result of the operation.
   \throw DACE::DACEException
   \note If estimation is not possible, zero is returned for all requested orders.
   \sa DA::estim
 */
    return da.estimNorm(err, var, type, nc);}

Interval bound(const DA &da){
/*! Compute lower and upper bounds of a DA object.
   \param[in] da a given DA object.
   \return An Interval object containing both the lower and the upper bound
    of the DA object.
   \throw DACE::DACEException
   \sa DA::bound
 */
    return da.bound();}

double convRadius(const DA &da, const double eps, const unsigned int type){
/*! Estimate the convergence radius of the given DA.
   \param[in] da the given DA object.
   \param[in] eps requested tolerance.
   \param[in] type type of norm (sum norm is used as default)
   \return A double corresponding to the estimated convergence radius.
   \throw DACE::DACEException
   \sa DA::conv_radius
 */
    return da.convRadius(eps, type);}

DA plug(const DA &da, const unsigned int var, const double val){
/*! Partial evaluation of a DA object. In the DA object, variable var is
    replaced by the value val. The resulting DA object is returned.
   \param[in] da a given DA object.
   \param[in] var variable number to be replaced
   \param[in] val value by which to replace the variable
   \return A new DA object containing the resulting DA object.
   \throw DACE::DACEException
   \sa DA::plug
 */
    return da.plug(var, val);}


DA replaceVariable(const DA &da, const unsigned int from, const unsigned int to, const double val){
/*! Partial evaluation of a DA object. In the DA object, variable from is
    replaced by the value val times variable to. The resulting DA object is returned.
   \param[in] da a given DA object.
   \param[in] from variable number to be replaced
   \param[in] to variable number to be inserted instead
   \param[in] val value by which to scale the inserted variable
   \return A new DA object containing the resulting DA object.
   \throw DACE::DACEException
   \sa DA::replaceVariable
 */
    return da.replaceVariable(from, to, val);}

DA scaleVariable(const DA &da, const unsigned int var, const double val){
/*! Scaling of an independent variable. In the DA object, variable var is
    replaced by the value val times var. The resulting DA object is returned.
   \param[in] da a given DA object.
   \param[in] var variable number to be scaled
   \param[in] val value by which to scale the variable
   \return A new DA object containing the resulting DA object.
   \throw DACE::DACEException
   \sa DA::scaleVariable
 */
    return da.scaleVariable(var, val);}

DA translateVariable(const DA &da, const unsigned int var, const double a, const double c){
/*! Affine translation of an independent variable. In the DA object, variable var is
    replaced by a*var + c. The resulting DA object is returned.
   \param[in] da a given DA object.
   \param[in] var variable number to be translated
   \param[in] a value by which to scale the variable
   \param[in] c value by which to shift the variable
   \return A new DA object containing the resulting DA object.
   \throw DACE::DACEException
   \sa DA::translateVariable
 */
    return da.translateVariable(var, a, c);}

compiledDA compile(const DA &da){
/*! Compile a given DA object and create a compiledDA object.
   \return A compiled DA object.
   \throw DACE::DACEException
   \sa DA::compile
 */
    return da.compile();}

std::string toString(const DA &da){
/*! Convert DA object to a string.
   \param[in] da a given DA object.
   \return A string.
   \throw DACE::DACEException
   \sa DA::toString
 */
    return da.toString();}


void write(const DA &da, std::ostream &os){
/*! Write binary representation of DA object to given stream os.
   \param[in] da a given DA object.
   \param[in] os the output stream to write to.
   \throw DACE::DACEException
   \sa DA::write
 */
    return da.write(os);}


// static class variables
const unsigned int storedDA::headerSize = daceBlobSize(NULL);

// create new storedDA from an existing DA
storedDA::storedDA(const DA &da){
    unsigned int len;

    daceExportBlob(da.m_index, NULL, len);
    resize(len);
    daceExportBlob(da.m_index, data(), len);
    if(daceGetError()) DACEException();
}

// create new storedDA by copying from a buffer
storedDA::storedDA(const std::vector<char> &data) : std::vector<char>(data){
}

// create a new storedDA by reading from a stream
storedDA::storedDA(std::istream &is) : std::vector<char>(storedDA::headerSize){
    // read blob header
    is.read(data(), headerSize);
    if(is.gcount() != headerSize)
    {
        resize((size_t)is.gcount());
        return;
    }

    // check validity of blob header and read the rest
    const unsigned int len = daceBlobSize(data());
    if(len == 0)
    {
        return;
    }
    else if(len > headerSize)
    {
        resize(len);
        is.read(data()+headerSize, len-headerSize);
        if(is.gcount() != (len-headerSize))
        {
            resize((size_t)(headerSize+is.gcount()));
            return;
        }
    }
}

// return if this storedDA data appears to be valid
bool storedDA::isValid() const{
    const size_t s1 = size();
    // check first that we have the minimum amount of data for the DACE blob header
    if(s1 < headerSize)
        return false;

    const unsigned int s2 = daceBlobSize(data());
    // is the blob valid
    if(s2 == 0)
        return false;

    // do we have the amount of data claimed in the header
    return s1 >= s2;
}

// convert data to string
storedDA::operator std::string() const{
    return std::string(data(), size());
}

// convert to DA
storedDA::operator DA() const{
    DA da;

    if(isValid())
    {
        daceImportBlob(data(), da.m_index);
        if(daceGetError()) DACEException();
    }
    else
    {
        DACEException(15, 111);     // XXX: number
    }

    return da;
}

// write binary data to ostream
std::ostream& operator<<(std::ostream &out, const storedDA &sda) {
    out.write(sda.data(), sda.size());

    return out;
}

}
