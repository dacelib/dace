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
 *  daceerror.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/** \addtogroup DACE Core 
 *  @{
 */

#include <stdio.h>

#include "dacebase.h"
#include "daceaux.h"
#include "daceerror.h"


/********************************************************************************
 *     DACE error state routine
 ********************************************************************************/

/*! Return the current error XYY code.
   \return The error code in XYY form
   \sa daceerror.h
*/
unsigned int daceGetError()
{
    return DACEDbg.ierr;
}

/*! Return the current error X code.
   \return The X error code X
   \sa daceerror.h
*/
unsigned int daceGetErrorX()
{
    return DACEDbg.ixerr;
}

/*! Return the current error YY code.
   \return The YY error code
   \sa daceerror.h
*/
unsigned int daceGetErrorYY()
{
    return DACEDbg.iyyerr;
}

/*! Return the function name of current generated error.
   \return The function name originating the error.
*/
const char* daceGetErrorFunName()
{
    return DACEDbg.name;
}

/*! Return the current error message.
   \return The current error message string
   \sa daceerror.h
*/
const char* daceGetErrorMSG()
{
    return DACEDbg.msg;
}

/*     ***********************
 *
 *     THIS SUBROUTINE CLEARS THE ERROR STATE.
 *     Should be called by the interface once the error has been corrected
 *
 */
/*! Clear the error code.
*/
void daceClearError()
{
    DACEDbg.ierr = 0;
    DACEDbg.ixerr = 0;
    DACEDbg.iyyerr = 0;
    strncpy(DACEDbg.name,"",ERROR_FUN_SIZE);
    strncpy(DACEDbg.msg,"",ERROR_MSG_SIZE);
}

/*!   This subroutine serves as an error handler for errors within the dace.      
      It is intended mostly for development and debugging. More descriptive error messages should be
      displayed by the user interface.
 
      The error codes are defined as XYY with X indicating the severity and
      YY corresponding to the actual error code
 
      Severity Levels X

      1  = Info:  Informative, no action required
      
      3  = Warning:  Serious, possibly incorrect use of DACE routines
      
      6  = Error:    Recoverable, result may not be correct or assumptions have
                     been made
      
      9  = Error:    Unrecoverable, new call to DACEINI is required to
                     reinitialize DACE, DACE objects are no longer valid
      
      10 = Critical: Crash in the DACE, just printing as much as possible
                     and dying
 
      Currently used error codes XYY are defined in daceerror.h
      
      \param[in] c an error string representing the error
      \param[in] ix is the error severity code
      \param[in] iyy is the error code
      \sa daceerror.h
 */

void daceSetError(const char *c, const unsigned int ix, const unsigned int iyy)
{
    // check if it is a critical error
    if(ix == DACE_PANIC)
    {
        fprintf(stderr, "DACE critical error %u in %s:\n%s\nbye bye!\n", DACEerr[iyy].ID, c, DACEerr[iyy].msg);
        exit(1);
    }
    else
    {
        if(ix>DACEDbg.ixerr)
        {
            DACEDbg.ierr = ix*100+iyy;      // set error code xyy
            DACEDbg.ixerr = ix;
            DACEDbg.iyyerr = iyy;
            strncpy(DACEDbg.name,c,ERROR_FUN_SIZE);
            strncpy(DACEDbg.msg,c,ERROR_FUN_SIZE);
            strncat(DACEDbg.msg,": ",ERROR_MSG_SIZE-strlen(DACEDbg.name));
            strncat(DACEDbg.msg,DACEerr[DACEDbg.iyyerr].msg,ERROR_MSG_SIZE-strlen(DACEDbg.name)-2);
        }
    }
}
/** @}*/