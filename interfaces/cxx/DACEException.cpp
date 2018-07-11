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
 * DACEException.cpp
 *
 *  Created on: Mar 11, 2014
 *      Author: Dinamica Srl
 */

// C++ stdlib classes used only internally in the implementation
#include <sstream>
#include <iostream>

// DACE classes
#include "DA/DACEException.h"
#include <dacecore.h>

namespace DACE{

    int DACEException::severity = 6;
    bool DACEException::warning = true;

    /********************************************************************************
    *     Constructors & Destructors
    *********************************************************************************/
    DACEException::DACEException(){
    /*! Create a DACEException object from an existing severity and error codes.
     */
        //m_x = FC_GLOBAL(dacegetxerr,DACEGETXERR)();
        //m_yy = FC_GLOBAL(dacegetyyerr,DACEGETYYERR)();
        //FC_GLOBAL(daceclrerr,DACECLRERR)();
        m_x = daceGetErrorX();
        m_yy = daceGetErrorYY();
        updateMessage();
        daceClearError();
        execute();
    }

    DACEException::DACEException(const int exc_sv, const int exc_id){
    /*! Create a DACEException object from given severity and ID codes.
       \param exc_sv severity code of the error.
       \param exc_id ID code of the error.
     */
        m_x = exc_sv;
        m_yy = exc_id;
        updateMessage();
        execute();
    }

    DACEException::~DACEException() throw(){
    /*! Destroy the DACEException object.
     */
        // nothing to do, just overloading the virtual destructor of the parent class
    }

    /********************************************************************************
    *     Private member functions
    *********************************************************************************/
    void DACEException::updateMessage(){
    /*! Update the error message of this exception based on its ID.
     */
        struct errstrings{
            int ID;
            const char* msg;};

        static const errstrings DACEerr[] = {
            { 000, "DACE: Unknown DACE error. Contact Dinamica SRL for filing a bug report."},
/*
            { 101, "DACEDAL: Attempt to deallocate protected or unallocated variable, ignored"},
            { 102, "DACESETNOT: Truncation order set to value larger then maximum order, setting to maximum order"},
            { 103, "DACEEST: Not enough non-zero monomials found, returned estimate may be inaccurate"},
            { 104, "DACEREAD: Line numbering out of order, ignoring line numbers"},
            { 105, "DACEREAD: DA vector contains more variables than current setup, ignoring coefficient"},
            { 106, "DACEREAD: Duplicate monomial in DA vector"},
            { 107, "DACEINI: Requested order has been increased to the minimum required order of 1"},
            { 108, "DACEINI: Requested number of variables has been increased to the minimum required number of variables 1"},

            { 601, "DACEVAR: Requested independent variable is out of bounds of current setup, returning zero DA"},
            { 602, "DACEPOK: Not enough storage to insert monomials, truncating"},
            { 603, "DACEVAR: Not enough storage to set up variable, returning zero DA"},
            { 604, "DACECOEF: Not enough storage to set coefficient, truncating"},
            { 605, "DACEDIVC: Divide by zero, returning zero DA"},
            { 606, "DACEINT: Requested independent variable out of bounds of current setup, returning zero DA"},
            { 607, "DACEDER: Requested independent variable out of bounds of current setup, returning zero DA"},
            { 608, "DACEMINV: Divide by zero, returning zero DA"},
            { 609, "DACESQRT: Negative constant part in square root, returning zero DA"},
            { 610, "DACEISQRT: Negative constant part in inverse square root, returning zero DA"},
            { 611, "DACELOG: Negative constant part in logarithm, returning zero DA"},
            { 612, "DACETAN: Cosine is zero in tangent, returning zero DA"},
            { 613, "DACEASIN: Constant part is out of domain [-1,1] in arcsine, returning zero DA"},
            { 614, "DACEACOS: Constant part is out of domain [-1,1] in arccosine, returning zero DA"},
            { 617, "DACEACOSH: Constant part is out of domain [1,infinity) in hyperbolic arccosine, returning zero DA"},
            { 618, "DASEATANH: Constant part is out of domain [-1,1], returning zero DA"},
            { 619, "DACEONORM: Requested independent variable out of bounds of current setup, returning zero"},
            { 620, "DACEEST: Maximum order must be at least 2 in order to use DA estimation"},
            { 622, "DACEPLUG: Requested independent variable out of bounds of current setup, returning zero DA"},
            { 623, "DACELOGB: Logarithm base must be positive, returning zero DA"},
            { 624, "DACEREAD: Not enough lines provided to read a DA vector, returning zero DA"},
            { 625, "DACEREAD: Unrecognized DA input format, returning zero DA"},
            { 627, "DACEENC: Invalid exponents with order larger than maximum order, returning zero"},
            { 628, "DACEDEC: Invalid DA codes provided, returning zero"},
            { 629, "DACEROOT: Zero-th root does not exists, returning zero DA"},
            { 630, "DACEROOT: Negative or zero constant part in even root, returning zero DA"},
            { 631, "DACEROOT: Zero constant part in odd root, returning zero DA"},
            { 632, "DACEPAC: Not enough storage in the target object, truncating"},
            { 633, "DACECMUL: Not enough storage in the target object, truncating"},
            { 634, "DACELIN: Not enough storage in the target object, truncating"},
            { 635, "DACEINT: Not enough storage in the target object, truncating"},
            { 636, "DACEDER: Not enough storage in the target object, truncating"},
            { 637, "DACECOP: Not enough storage in the target object, truncating"},
            { 638, "DACEPUSH: Not enough space to store the provided data, truncating"},
            { 639, "DACEPULL: Not enough space to store the requested data, truncating"},
            { 640, "DACETRIM: Not enough space to store the requested data, truncating"},
            { 641, "DACENORM: Unknown norm type, resorting to max norm"},
            { 642, "DACEONORM: Unknown norm type, resorting to max norm"},
            { 643, "DACETRIM: Not enough storage in the target object, truncating"},
            { 644, "DACEWRITE: DA vector contains more monomials than expected, truncating"},

            { 901, "DACEINI: Internal array size exceeds available addressing space"},

            {1098, "DACEINI: Internal error, size of generated DA coding arrays does not match number of monomials"},
            {1097, "DACEINI: Internal error, generated DA coding arrays are faulty"},
            {1096, "DACEINI: Internal error, unable to correctly allocate scratch variables"},
            {1094, "DACEINI: Internal error, memory allocation failure"},
            {1093, "DACEREALL: Internal error, memory allocation failure"},
            {1092, "DACEALL: DACE has not been initialized"},
            {1090, "DACEINFO: DA Object not allocated"},
*/
            {1101, "DA::getCoeff: Not enough exponents, missing exponents treated as zero"},
            {1102, "DA::setCoeff: Not enough exponents, missing exponents treated as zero"},
            {1103, "DA::getCoeff: More exponents than variables, ignoring extra exponents"},
            {1104, "DA::setCoeff: More exponents than variables, ignoring extra exponents"},

            {1604, "compiledDA::compiledDA(): Dimension lower than 1"},
            {1605, "compiledDA::eval: Argument size lower than the number of variables in the polynomial"},

            {2099, "DA::checkVersion: DACE C++ interface header file and DACE core library version mismatch"}
        };
        static const int length = sizeof(DACEerr)/sizeof(errstrings);

        const int id = m_x*100+m_yy;
        int i = 0;
        std::stringstream s;
        
        if(m_x > 10)
        {
            for(i=length-1; (i>0)&&(DACEerr[i].ID != id); i--);

            
            s << DACEerr[i].msg << " (ID: " << DACEerr[i].ID << ")" ;
            
        }
        else
        {
            s << daceGetErrorMSG() << " (ID: " << id << ")";
        }
        msg = s.str();
    }

    void DACEException::execute() const{
    /*! Execute this exception, i.e. throw or print warning based on current user settings.
       \throw DACE::DACEException
     */
        const int sev = m_x%11; // modulo 11 to handle both DACE and C++ interface severity codes

        if(sev>=severity){
            throw *this;}
        else if(warning){
            std::cerr << "Warning: " << msg << std::endl;}
    }

    /********************************************************************************
    *     Public member functions
    *********************************************************************************/
    const char* DACEException::what() const throw(){
    /*! Return a human readable error string representing this exception.
       \return A C string containing the error message.
     */
        return msg.c_str();
    }

    /********************************************************************************
    *     Static member functions
    *********************************************************************************/
    void DACEException::setSeverity(const int n){
    /*! Set a value for the severity level. Errors with severity code greater
        or equal to this value will generate an exception object.
        Severity levels are:\n
        0 = Warning:   Informative, no action required\n
        1 = Warning:   Serius, possible wrong implementation\n
        6 = Error:     Recoverable, assumptions have been made\n
        9 = Error:     Unrecoverable, new call to init is required to
                       reinitialize DACE, interface objects are no longer valid\n
        10 = Critical: Crash in the DACE, just printing as much as possible
                       and dying.
       \param n severity value
     */
        severity = n;
    }

    void DACEException::setWarning(const bool w){
    /*! Set the current status for visualization of warnings.
       \param w status of warnings\n
        true: turn on printing of warnings\n
        false: turn off printing of warnings
     */
        warning = w;
    }

    /********************************************************************************
    *     Friend functions
    *********************************************************************************/
	DACE_API std::ostream& operator<< (std::ostream &out, const DACEException &ex){
    /*! Overload of std::operator<< in iostream.
       \param[in] out standard output stream.
       \param[in] ex Exception to be printed in the stream
       \return Output stream out.
     */
        return out << ex.msg << std::endl;
    }

}
