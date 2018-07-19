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
 * DACEException.h
 *
 *  Created on: Mar 11, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_DACEEXCEPTION_H_
#define DINAMICA_DACEEXCEPTION_H_

// C++ stdlib classes used in this public interface
#include <exception>
#include <string>
#include <ostream>

#include "dace/Def.h"

namespace DACE{

/*! DACEException class containing methods for error handling within the DACE C++ interface. */
class DACE_API DACEException : public std::exception
{
private:
    int m_x;                    //!< Severity code
    int m_yy;                   //!< Error code
	std::string msg;            //!< Error message
    static int severity;        //!< Default severity code
    static bool warning;        //!< Default warning status
    DACE_LOCAL void execute() const;       //!< Execute the exception
    DACE_LOCAL void updateMessage();       //!< Update the error message

public:
    DACEException();                                    //!< Default constructor
    DACEException(const int exc_sv, const int exc_id);  //!< Constructor
    ~DACEException() throw();                           //!< Destructor

    const char* what() const throw();                   //!< Convert exception to string
    static void setSeverity(const int n);               //!< Select the desired severity code
    static void setWarning(const bool w);               //!< Select the warning status

	friend DACE_API std::ostream& operator<< (std::ostream &out, const DACEException &ex); //!< Overload output stream operator
};

}
#endif /* DINAMICA_DACEEXCEPTION_H_ */
