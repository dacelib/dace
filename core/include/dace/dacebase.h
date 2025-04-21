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
 *  dacebase.h
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/*
    This file contains all routines in the public interface to the DACE core.
    DACE core users or high level interfaces should only include this file.

    Legacy users of DACE 1 can additionally include dacecompat.h to obtain a
    source compatible mapping to the old DACE 1 function names and semantics.
*/
/** \addtogroup DACE Core
 *  @{
 */

#ifndef DINAMICA_DACEBASE_H_
#define DINAMICA_DACEBASE_H_

#include <stdbool.h>        // for bool type

#include "dace/config.h"

// Maximum line length of the DACE string I/O interface
#define DACE_STRLEN (140)

#ifdef __cplusplus
    // in C++ pass values by reference (which is ABI compatible with C pointers in all supported C++ compilers)
    #define REF(x) &(x)
    #define ARG(x) (x)
    extern "C" {
#else /* __cplusplus */
    // in C, and hence when compiling the DACE core library, pass values by pointer
    #define REF(x) *(x)
    #define ARG(x) &(x)
#endif  /* __cplusplus */

// Error handling symbolic constants
#define DACE_INFO       1
#define DACE_WARNING    3
#define DACE_ERROR      6
#define DACE_SEVERE     9
#define DACE_PANIC     10

#if DACE_MEMORY_MODEL == DACE_MEMORY_HYBRID || DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    // Type of a DACE DA object
    typedef int DACEDA;
#elif DACE_MEMORY_MODEL == DACE_MEMORY_DYNAMIC
    // A DACE variable
    typedef struct dvariable {
        unsigned int len, max;
        struct dmonomial *mem;
    } variable;

    // Type of a DACE DA object
    typedef variable DACEDA;
#else
    #error Invalid DACE memory model selected!
#endif

/// @cond
/********************************************************************************
*     DACE initialization and state related routines
*********************************************************************************/
DACE_API void daceInitialize(unsigned int no, unsigned int nv);
DACE_API void daceInitializeThread();
DACE_API void daceCleanupThread();
DACE_API void daceGetVersion(int REF(imaj), int REF(imin), int REF(ipat));
DACE_API double daceSetEpsilon(const double deps);
DACE_API double daceGetEpsilon();
DACE_API double daceGetMachineEpsilon();
DACE_API unsigned int daceGetMaxOrder();
DACE_API unsigned int daceGetMaxVariables();
DACE_API unsigned int daceGetMaxMonomials();
DACE_API unsigned int daceGetTruncationOrder();
DACE_API unsigned int daceSetTruncationOrder(const unsigned int fnot);

/********************************************************************************
*     DACE error state routine
*********************************************************************************/
DACE_API unsigned int daceGetError();
DACE_API unsigned int daceGetErrorX();
DACE_API unsigned int daceGetErrorYY();
DACE_API const char* daceGetErrorFunName();
DACE_API const char* daceGetErrorMSG();
DACE_API void daceClearError();

/********************************************************************************
*     DACE memory handling routines
*********************************************************************************/
DACE_API void daceAllocateDA(DACEDA REF(inc), const unsigned int len);
DACE_API void daceFreeDA(DACEDA REF(inc));
DACE_API void daceInvalidateDA(DACEDA REF(inc));
DACE_API void daceMemoryDump();     // really just an internal debugging routine

/********************************************************************************
*     DACE variable creation routines
*********************************************************************************/
DACE_API void daceCreateVariable(DACEDA REF(ina), const unsigned int i, const double ckon);
DACE_API void daceCreateMonomial(DACEDA REF(ina), const unsigned int jj[], const double ckon);
DACE_API void daceCreateConstant(DACEDA REF(ina), const double ckon);
DACE_API void daceCreateFilled(DACEDA REF(ina), const double ckon);
DACE_API void daceCreateRandom(DACEDA REF(ina), const double cm);

/********************************************************************************
*     DACE coefficient access routines
*********************************************************************************/
DACE_API double daceGetConstant(const DACEDA REF(ina));
DACE_API void daceGetLinear(const DACEDA REF(ina), double c[]);
DACE_API double daceGetCoefficient(const DACEDA REF(ina), const unsigned int jj[]);
DACE_API double daceGetCoefficient0(const DACEDA REF(ina), const unsigned int ic);
DACE_API void daceSetCoefficient(DACEDA REF(ina), const unsigned int jj[], const double cjj);
DACE_API void daceSetCoefficient0(DACEDA REF(ina), const unsigned int ic, const double cjj);
DACE_API void daceGetCoefficientAt(const DACEDA REF(ina), const unsigned int npos, unsigned int jj[], double REF(cjj));
DACE_API unsigned int daceGetLength(const DACEDA REF(ina));

/********************************************************************************
*     DACE DA copying and filtering
*********************************************************************************/
DACE_API void daceCopy(const DACEDA REF(ina), DACEDA REF(inb));
DACE_API void daceCopyFiltering(const DACEDA REF(ina), DACEDA REF(inb));
DACE_API void daceFilter(const DACEDA REF(ina), DACEDA REF(inb), const DACEDA REF(inc));
DACE_API void daceTrim(const DACEDA REF(ina), const unsigned int imin, const unsigned int imax, DACEDA REF(inc));
DACE_API unsigned int daceIsNan(const DACEDA REF(ina));
DACE_API unsigned int daceIsInf(const DACEDA REF(ina));

/********************************************************************************
*     Basic DACE arithmetic operations
*********************************************************************************/
DACE_API void daceWeightedSum(const DACEDA REF(ina), const double afac, const DACEDA REF(inb), const double bfac, DACEDA REF(inc));
DACE_API void daceAdd(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
DACE_API void daceSubtract(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
DACE_API void daceMultiply(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
DACE_API void daceMultiplyMonomials(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
DACE_API void daceDivide(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
DACE_API void daceSquare(const DACEDA REF(ina), DACEDA REF(inb));
DACE_API void daceAddDouble(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
DACE_API void daceDoubleSubtract(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
DACE_API void daceSubtractDouble(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
DACE_API void daceMultiplyDouble(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
DACE_API void daceDivideDouble(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
DACE_API void daceDoubleDivide(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
DACE_API void daceDivideByVariable(const DACEDA REF(ina), const unsigned int var, const unsigned int p, DACEDA REF(inc));
DACE_API void daceDifferentiate(const unsigned int idif, const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceIntegrate(const unsigned int iint, const DACEDA REF(ina), DACEDA REF(inc));

/********************************************************************************
*     DACE intrinsic function routines
*********************************************************************************/
DACE_API void daceTruncate(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceRound(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceModulo(const DACEDA REF(ina), const double p, DACEDA REF(inc));
DACE_API void dacePowerDouble(const DACEDA REF(ina), const double p, DACEDA REF(inc));
DACE_API void dacePower(const DACEDA REF(ina), const int np, DACEDA REF(inc));
DACE_API void daceRoot(const DACEDA REF(ina), const int np, DACEDA REF(inc));
DACE_API void daceMultiplicativeInverse(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceSquareRoot(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceInverseSquareRoot(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceCubicRoot(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceInverseCubicRoot(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceHypotenuse(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
DACE_API void daceExponential(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceLogarithm(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceLogarithmBase(const DACEDA REF(ina), const double b, DACEDA REF(inc));
DACE_API void daceLogarithm10(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceLogarithm2(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceSine(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceCosine(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceTangent(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceArcSine(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceArcCosine(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceArcTangent(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceArcTangent2(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
DACE_API void daceHyperbolicSine(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceHyperbolicCosine(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceHyperbolicTangent(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceHyperbolicArcSine(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceHyperbolicArcCosine(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceHyperbolicArcTangent(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceErrorFunction(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceComplementaryErrorFunction(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceBesselIFunction(const DACEDA REF(ina), const int n, const bool scaled, DACEDA REF(inc));
DACE_API void daceBesselJFunction(const DACEDA REF(ina), const int n, DACEDA REF(inc));
DACE_API void daceBesselKFunction(const DACEDA REF(ina), const int n, const bool scaled, DACEDA REF(inc));
DACE_API void daceBesselYFunction(const DACEDA REF(ina), const int n, DACEDA REF(inc));
DACE_API void daceLogGammaFunction(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void daceGammaFunction(const DACEDA REF(ina), DACEDA REF(inc));
DACE_API void dacePsiFunction(const DACEDA REF(ina), const unsigned int n, DACEDA REF(inc));

/********************************************************************************
*     DACE norm and norm estimation routines
*********************************************************************************/
DACE_API double daceAbsoluteValue(const DACEDA REF(ina));
DACE_API double daceNorm(const DACEDA REF(ina), const unsigned int ityp);
DACE_API void daceOrderedNorm(const DACEDA REF(ina), const unsigned int ivar, const unsigned int ityp, double onorm[]);
DACE_API void daceEstimate(const DACEDA REF(ina), const unsigned int ivar, const unsigned int ityp, double c[], double err[], const unsigned int nc);
DACE_API void daceGetBounds(const DACEDA REF(ina), double REF(alo), double REF(aup));

/********************************************************************************
*     DACE polynomial evaluation routines
*********************************************************************************/
DACE_API double daceEvalMonomials(const DACEDA REF(ina), const DACEDA REF(inb));
DACE_API void daceReplaceVariable(const DACEDA REF(ina), const unsigned int from, const unsigned int to, const double val, DACEDA REF(inc));
DACE_API void daceEvalVariable(const DACEDA REF(ina), const unsigned int nvar, const double val, DACEDA REF(inc));
DACE_API void daceScaleVariable(const DACEDA REF(ina), const unsigned int nvar, const double val, DACEDA REF(inc));
DACE_API void daceTranslateVariable(const DACEDA REF(ina), const unsigned int nvar, const double a, const double c, DACEDA REF(inc));
DACE_API void daceEvalTree(const DACEDA *das[], const unsigned int count, double ac[], unsigned int REF(nterm), unsigned int REF(nvar), unsigned int REF(nord));

/********************************************************************************
*     DACE input/output routines
*********************************************************************************/
DACE_API void daceWrite(const DACEDA REF(ina), char *strs, unsigned int REF(nstrs));
DACE_API void daceRead(DACEDA REF(ina), char *strs, unsigned int nstrs);
DACE_API void dacePrint(const DACEDA REF(ina));
DACE_API unsigned int daceExportBlob(const DACEDA REF(ina), void *blob, unsigned int REF(size));
DACE_API unsigned int daceBlobSize(const void *blob);
DACE_API void daceImportBlob(const void *blob, DACEDA REF(inc));

/********************************************************************************
*     DACE miscellaneous routines
*********************************************************************************/
DACE_API double daceRandom();
/// @endcond
#ifdef __cplusplus
    }
#endif /* _cplusplus */
/** @}*/
#endif /* DINAMICA_DACEBASE_H_ */
