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
 * Def.h
 *
 *  Created on: Aug 20, 2015
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_DEF_H_
#define DINAMICA_DEF_H_

// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define DACE_HELPER_DLL_IMPORT __declspec(dllimport)
  #define DACE_HELPER_DLL_EXPORT __declspec(dllexport)
  #define DACE_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define DACE_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define DACE_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define DACE_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define DACE_HELPER_DLL_IMPORT
    #define DACE_HELPER_DLL_EXPORT
    #define DACE_HELPER_DLL_LOCAL
  #endif
#endif

// Now we use the generic helper definitions above to define DACE_API and DACE_LOCAL.
// DACE_API is used for the public API symbols. It either DLL imports or DLL exports (or does nothing for static build)
// DACE_LOCAL is used for non-api symbols.

#ifdef DACE_DLL // defined if DACE is compiled as a DLL
  #ifdef DACE_DLL_EXPORTS // defined if we are building the DACE DLL (instead of using it)
    #define DACE_API DACE_HELPER_DLL_EXPORT
  #else
    #define DACE_API DACE_HELPER_DLL_IMPORT
  #endif // DACE_DLL_EXPORTS
  #define DACE_LOCAL DACE_HELPER_DLL_LOCAL
#else // DACE_DLL is not defined: this means DACE is a static lib.
  #define DACE_API
  #define DACE_LOCAL
#endif // DACE_DLL

#endif //DINAMICA_DEF_H_
