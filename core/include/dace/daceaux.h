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
 *  daceaux.h
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/*
    This file contains all internal DACE auxiliary functions used by the DACE core.
    It is not meant to be included publicly by DACE users or high level interfaces.
*/
/** \addtogroup DACE Core
 *  @{
 */

#ifndef DINAMICA_DACEAUX_H_
#define DINAMICA_DACEAUX_H_

#include <stdlib.h>     // for size_t
#include <stdbool.h>    // for bool
#include <string.h>

#include "dace/config.h"

// macro to help tell supported compiler which branch is most likely
#if __GNUC__ || __clang__
    #define UNLIKELY(expr) __builtin_expect(!!(expr), 0)
    #define LIKELY(expr) __builtin_expect(!!(expr), 1)
#else
    #define UNLIKELY(expr) expr
    #define LIKELY(expr) expr
#endif

// DACE internal data structure
typedef struct dcom {
    unsigned int *ie1, *ie2, *ieo, *ia1, *ia2;
    unsigned int nomax, nvmax, nv1, nv2, nmmax;
    double epsmac;
} dacecom;

// DACE thread local data structure
typedef struct dcom_t {
    unsigned int nocut;
    double eps;
#ifdef DACE_FILTERING
    unsigned int *ifi;
    unsigned int lfi;
#endif
} dacecom_t;

#define ERROR_FUN_SIZE 64
#define ERROR_MSG_SIZE 256
// DACE error management
typedef struct ddbg {
    unsigned int ierr, ixerr, iyyerr;
    char name[ERROR_FUN_SIZE];
    char msg[ERROR_MSG_SIZE];
} dacedbg;

// Basic memory structure of a monomial
typedef struct dmonomial {
    double cc;
    unsigned int ii;
} monomial;

// Basic memory structure of an extended monomial (these are written and read from binary files, so enforce no padding!)
#pragma pack(push,1)
typedef struct dextendedmonomial {
    unsigned int i1, i2;
    double cc;
} extended_monomial;
#pragma pack(pop)

// 32 bit magic marker to identify beginning of binary DACE blobs (in little endian the letters DA0 and the record separator 1E)
#define DACE_BINARY_MAGIC (0x1E304144)

// Memory structure of a DACE blob (not public!)
#pragma pack(push,1)
struct daceblob {
    unsigned int magic;
    unsigned int no, nv1, nv2, len;
    extended_monomial monomials[1];
};
#pragma pack(pop)

// external global symbols for common DACE structures (actually allocated in dacememory.c)
extern dacecom DACECom;
extern DACE_THREAD_LOCAL dacecom_t DACECom_t;
extern DACE_THREAD_LOCAL dacedbg DACEDbg;

/// @cond
// math utility routines
unsigned int umin(const unsigned int a, const unsigned int b);
unsigned int umax(const unsigned int a, const unsigned int b);
double pown(double a, unsigned int b);
int npown(int a, unsigned int b);

#if DACE_MEMORY_MODEL == DACE_MEMORY_HYBRID || DACE_MEMORY_MODEL == DACE_MEMORY_DYNAMIC
    // dynamic memory allocation wrappers
    void* dacecalloc(size_t count, size_t size);
    void* dacemalloc(size_t size);
    void* dacemalloc0(size_t size);
    void* dacerealloc(void* ptr, size_t size);
    void dacefree(void* ptr);
#endif

// internal memory related routines
void daceFreeMemory();
void daceVariableInformation(const DACEDA *inc, monomial **ipoc, unsigned int *ilmc, unsigned int *illc);
void daceSetLength(DACEDA *inc, const size_t len);
bool daceIsSameObject(const DACEDA *ina, const DACEDA *inb);

// basic monomial coding integer related routines
unsigned int daceEncodeExponents(const unsigned int p[], const unsigned int no, const unsigned int nv);
unsigned int daceEncode(const unsigned int p[]);
unsigned int daceDecodeExponents(unsigned int ic, const unsigned int no, const unsigned int nv, unsigned int p[]);
void daceDecode(const unsigned int ic, unsigned int p[]);
unsigned int daceCountMonomials(unsigned int no, unsigned int nv);
unsigned int daceNextMonomial(unsigned int p[], const unsigned int no, const unsigned int nv);
unsigned int daceNextOrderedMonomial(unsigned int p[], const unsigned int no, const unsigned int nv);

// internal routines
void daceInitializeThread0();
void daceSetError(const char *c, const unsigned int ix, const unsigned int iyy);
void dacePack(double *restrict cc, DACEDA *restrict inc);
void daceMultiplicativeInverse0(const DACEDA *ina, DACEDA *inc, const double a0);
int BesselWrapper(const double x, const int n0, const int n1, const int type, double *bz);
int ModifiedBesselWrapper(const double x, const int n0, const int n1, const int type, double *bz);
void daceEvaluateBesselFunction(const DACEDA *ina, const double bz[], const double type, const double ktype, DACEDA *inc);
void daceEvaluateScaledModifiedBesselFunction(const DACEDA *ina, const double bz[], const double type, DACEDA *inc);
void daceLogGammaFunction0(const DACEDA *ina, const double a0, DACEDA *inc);
void daceEvaluateSeries(const DACEDA *ina, const double xf[], DACEDA *inc);
/// @endcond
/** @}*/
#endif /* DINAMICA_DACEAUX_H_ */
