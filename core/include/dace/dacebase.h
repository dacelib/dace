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

// DACE version
#define DACE_MAJOR_VERSION 2
#define DACE_MINOR_VERSION 0
#define DACE_COSY          0

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
void daceInitialize(unsigned int no, unsigned int nv);              
void daceInitializeThread();                                        
void daceGetVersion(int REF(imaj), int REF(imin), int REF(icos));   
double daceSetEpsilon(const double deps);                           
double daceGetEpsilon();                                            
double daceGetMachineEpsilon();                                     
unsigned int daceGetMaxOrder();                                    
unsigned int daceGetMaxVariables();                                 
unsigned int daceGetMaxMonomials();                                 
unsigned int daceGetTruncationOrder();                              
unsigned int daceSetTruncationOrder(const unsigned int fnot);       

/********************************************************************************
*     DACE error state routine
*********************************************************************************/
unsigned int daceGetError();                                       
unsigned int daceGetErrorX();                                       
unsigned int daceGetErrorYY();                                      
const char* daceGetErrorFunName();                                  
const char* daceGetErrorMSG();                                      
void daceClearError();                                              

/********************************************************************************
*     DACE memory handling routines
*********************************************************************************/
void daceAllocateDA(DACEDA REF(inc), const unsigned int len);       
void daceFreeDA(DACEDA REF(inc));                                   
void daceInvalidateDA(DACEDA REF(inc));                             
void daceMemoryDump();     // really just an internal debugging routine

/********************************************************************************
*     DACE variable creation routines
*********************************************************************************/
void daceCreateVariable(DACEDA REF(ina), const unsigned int i, const double ckon);          
void daceCreateMonomial(DACEDA REF(ina), const unsigned int jj[], const double ckon);       
void daceCreateConstant(DACEDA REF(ina), const double ckon);                                
void daceCreateFilled(DACEDA REF(ina), const double ckon);                                  
void daceCreateRandom(DACEDA REF(ina), const double cm);                                    

/********************************************************************************
*     DACE coefficient access routines
*********************************************************************************/
double daceGetConstant(const DACEDA REF(ina));                                              
void daceGetLinear(const DACEDA REF(ina), double c[]);                                      
double daceGetCoefficient(const DACEDA REF(ina), const unsigned int jj[]);                  
double daceGetCoefficient0(const DACEDA REF(ina), const unsigned int ic);                   
void daceSetCoefficient(DACEDA REF(ina), const unsigned int jj[], const double cjj);       
void daceSetCoefficient0(DACEDA REF(ina), const unsigned int ic, const double cjj);         
void daceGetCoefficientAt(const DACEDA REF(ina), const unsigned int npos, unsigned int jj[], double REF(cjj));  
unsigned int daceGetLength(const DACEDA REF(ina));                                          

/********************************************************************************
*     DACE DA copying and filtering
*********************************************************************************/  
void daceCopy(const DACEDA REF(ina), DACEDA REF(inb));                                      
void daceCopyFiltering(const DACEDA REF(ina), DACEDA REF(inb));                            
void daceFilter(const DACEDA REF(ina), DACEDA REF(inb), const DACEDA REF(inc));             
void daceTrim(const DACEDA REF(ina), const unsigned int imin, const unsigned int imax, DACEDA REF(inc));

/********************************************************************************
*     Basic DACE arithmetic operations
*********************************************************************************/
void daceWeightedSum(const DACEDA REF(ina), const double afac, const DACEDA REF(inb), const double bfac, DACEDA REF(inc));
void daceAdd(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
void daceSubtract(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
void daceMultiply(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
void daceMultiplyMonomials(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
void daceDivide(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
void daceSquare(const DACEDA REF(ina), DACEDA REF(inb));
void daceAddDouble(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
void daceDoubleSubtract(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
void daceSubtractDouble(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
void daceMultiplyDouble(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
void daceDivideDouble(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
void daceDoubleDivide(const DACEDA REF(ina), const double ckon, DACEDA REF(inb));
void daceDivideByVariable(const DACEDA REF(ina), const unsigned int var, const unsigned int p, DACEDA REF(inc));
void daceDifferentiate(const unsigned int idif, const DACEDA REF(ina), DACEDA REF(inc));
void daceIntegrate(const unsigned int iint, const DACEDA REF(ina), DACEDA REF(inc));

/********************************************************************************
*     DACE intrinsic function routines
*********************************************************************************/
void daceTruncate(const DACEDA REF(ina), DACEDA REF(inc));
void daceRound(const DACEDA REF(ina), DACEDA REF(inc));
void daceModulo(const DACEDA REF(ina), const double p, DACEDA REF(inc));
void dacePowerDouble(const DACEDA REF(ina), const double p, DACEDA REF(inc));
void dacePower(const DACEDA REF(ina), const int np, DACEDA REF(inc));
void daceRoot(const DACEDA REF(ina), const int np, DACEDA REF(inc));
void daceMultiplicativeInverse(const DACEDA REF(ina), DACEDA REF(inc));
void daceSquareRoot(const DACEDA REF(ina), DACEDA REF(inc));
void daceInverseSquareRoot(const DACEDA REF(ina), DACEDA REF(inc));
void daceCubicRoot(const DACEDA REF(ina), DACEDA REF(inc));
void daceInverseCubicRoot(const DACEDA REF(ina), DACEDA REF(inc));
void daceHypotenuse(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
void daceExponential(const DACEDA REF(ina), DACEDA REF(inc));
void daceLogarithm(const DACEDA REF(ina), DACEDA REF(inc));
void daceLogarithmBase(const DACEDA REF(ina), const double b, DACEDA REF(inc));
void daceLogarithm10(const DACEDA REF(ina), DACEDA REF(inc));
void daceLogarithm2(const DACEDA REF(ina), DACEDA REF(inc));
void daceSine(const DACEDA REF(ina), DACEDA REF(inc));
void daceCosine(const DACEDA REF(ina), DACEDA REF(inc));
void daceTangent(const DACEDA REF(ina), DACEDA REF(inc));
void daceArcSine(const DACEDA REF(ina), DACEDA REF(inc));
void daceArcCosine(const DACEDA REF(ina), DACEDA REF(inc));
void daceArcTangent(const DACEDA REF(ina), DACEDA REF(inc));
void daceArcTangent2(const DACEDA REF(ina), const DACEDA REF(inb), DACEDA REF(inc));
void daceHyperbolicSine(const DACEDA REF(ina), DACEDA REF(inc));
void daceHyperbolicCosine(const DACEDA REF(ina), DACEDA REF(inc));
void daceHyperbolicTangent(const DACEDA REF(ina), DACEDA REF(inc));
void daceHyperbolicArcSine(const DACEDA REF(ina), DACEDA REF(inc));
void daceHyperbolicArcCosine(const DACEDA REF(ina), DACEDA REF(inc));
void daceHyperbolicArcTangent(const DACEDA REF(ina), DACEDA REF(inc));
void daceErrorFunction(const DACEDA REF(ina), DACEDA REF(inc));
void daceComplementaryErrorFunction(const DACEDA REF(ina), DACEDA REF(inc));
void daceBesselJFunction(const DACEDA REF(ina), const int n, DACEDA REF(inc));
void daceBesselYFunction(const DACEDA REF(ina), const int n, DACEDA REF(inc));

/********************************************************************************
*     DACE norm and norm estimation routines
*********************************************************************************/
double daceAbsoluteValue(const DACEDA REF(ina));
double daceNorm(const DACEDA REF(ina), const unsigned int ityp);
void daceOrderedNorm(const DACEDA REF(ina), const unsigned int ivar, const unsigned int ityp, double onorm[]);
void daceEstimate(const DACEDA REF(ina), const unsigned int ivar, const unsigned int ityp, double c[], double err[], const unsigned int nc);
void daceGetBounds(const DACEDA REF(ina), double REF(alo), double REF(aup));

/********************************************************************************
*     DACE polynomial evaluation routines
*********************************************************************************/
double daceEvalMonomials(const DACEDA REF(ina), const DACEDA REF(inb));
void daceReplaceVariable(const DACEDA REF(ina), const unsigned int from, const unsigned int to, const double val, DACEDA REF(inc));
void daceEvalVariable(const DACEDA REF(ina), const unsigned int nvar, const double val, DACEDA REF(inc));
void daceScaleVariable(const DACEDA REF(ina), const unsigned int nvar, const double val, DACEDA REF(inc));
void daceTranslateVariable(const DACEDA REF(ina), const unsigned int nvar, const double a, const double c, DACEDA REF(inc));
void daceEvalTree(const DACEDA *das[], const unsigned int count, double ac[], unsigned int REF(nterm), unsigned int REF(nvar), unsigned int REF(nord));

/********************************************************************************
*     DACE input/output routines
*********************************************************************************/
void daceWrite(const DACEDA REF(ina), char *strs, unsigned int REF(nstrs));
void daceRead(DACEDA REF(ina), char *strs, unsigned int nstrs);
void dacePrint(const DACEDA REF(ina));
unsigned int daceExportBlob(const DACEDA REF(ina), void *blob, unsigned int REF(size));
unsigned int daceBlobSize(const void *blob);
void daceImportBlob(const void *blob, DACEDA REF(inc));

/********************************************************************************
*     DACE miscellaneous routines
*********************************************************************************/
double daceRandom();
/// @endcond
#ifdef __cplusplus
    }
#endif /* _cplusplus */
/** @}*/
#endif /* DINAMICA_DACEBASE_H_ */
