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
 *  daceinit.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/** \addtogroup DACE Core
 *  @{
 */

#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "dace/config.h"
#include "dace/dacebase.h"
#include "dace/daceaux.h"

#ifdef WITH_PTHREAD
    #include <pthread.h>
#endif

/*! Set up the ordering and addressing arrays in the common data structure
    and initialize DA memory.
   \note MUST BE CALLED BEFORE ANY OTHER DA ROUTINE CAN BE USED.
   \note Also initializes the truncation order to the maximum computation order
   and disables the DA epsilon cutoff by setting it to 0.0.
   \param[in] no order of the Taylor polynomials;
   \param[in] nv number of variables considered.
   \sa daceTruncationOrder
   \sa daceSetEpsilon
 */
void daceInitialize(unsigned int no, unsigned int nv)
{
   // reset error state
   daceClearError();

    // check inputs
    if(no < 1) {
        daceSetError(__func__, DACE_INFO, 67);
        no = 1;
    }

    if(nv < 1) {
        daceSetError(__func__, DACE_INFO, 68);
        nv = 1;
    }

    // Find machine epsilon
    DACECom.epsmac = 1.0;
    while(1.0+DACECom.epsmac > 1.0)
    {
         DACECom.epsmac = DACECom.epsmac/2.0;
    }
    DACECom.epsmac = DACECom.epsmac*2.0;

    // Reset memory, purging all previous DA objects, if any.
    daceFreeMemory();

    // Compute and check required internal array sizes to be compatible with 32 bit address space
    const double clia = pown(no+1, (nv+1)/2);  // length of reverse lookup array

    if(clia >= pown(2.0, 32))
    {
         daceSetError(__func__, DACE_SEVERE, 11);
         return;
    }

    const unsigned int lia = (unsigned int)clia;  // length of reverse lookup array
    const unsigned int lea = daceCountMonomials(no, nv);     // number of monomials = length of forward lookup array

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    if(no > DACE_STATIC_NOMAX || nv > DACE_STATIC_NVMAX || lea > DACE_STATIC_NMMAX || lia > DACE_STATIC_LIAMAX)
    {
         daceSetError(__func__, DACE_SEVERE, 11);
         return;
    }
    // use static local memory for addressing arrays
    static unsigned int ie1[DACE_STATIC_NMMAX], ie2[DACE_STATIC_NMMAX], ieo[DACE_STATIC_NMMAX], ia1[DACE_STATIC_LIAMAX+1], ia2[DACE_STATIC_LIAMAX+1];
    DACECom.ie1 = ie1;
    DACECom.ie2 = ie2;
    DACECom.ieo = ieo;
    DACECom.ia1 = ia1;
    DACECom.ia2 = ia2;
#else
    // (re)allocate addressing arrays
    dacefree(DACECom.ie1);
    dacefree(DACECom.ie2);
    dacefree(DACECom.ieo);
    dacefree(DACECom.ia1);
    dacefree(DACECom.ia2);

    DACECom.ie1 = dacecalloc(lea, sizeof(unsigned int));
    DACECom.ie2 = dacecalloc(lea, sizeof(unsigned int));
    DACECom.ieo = dacecalloc(lea, sizeof(unsigned int));
    DACECom.ia1 = dacecalloc(lia+1, sizeof(unsigned int));
    DACECom.ia2 = dacecalloc(lia+1, sizeof(unsigned int));
#endif

    // Fill the addressing arrays.
    // This is a complete rewrite from scratch of the original algorithm by Berz.
    // It is probably slightly less efficient but much more readable.
    const unsigned int nv1 = (nv+1)/2, nv2 = nv-nv1;          // number of variables in first and second coding integer
    unsigned int i = 0;                                       // total number of monomials
    unsigned int no1 = 0, no2 = 0;

    // allocate and set exponents to zero
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
	unsigned int p1[DACE_STATIC_NVMAX] = {0}, p2[DACE_STATIC_NVMAX] = {0}; // array of powers of first and second part of monomial
#else
    unsigned int *p1, *p2;                                    // array of powers of first and second part of monomial
    p1 = dacecalloc(nv1, sizeof(unsigned int));
    p2 = dacecalloc(nv2, sizeof(unsigned int));
#endif

    // enumerate all monomials of nv1 variables up to order no in p1
    do {
        unsigned int exp1 = daceEncodeExponents(p1, no, nv1);  // exponents of the first nv1 variables
        unsigned int i0 = i;
        DACECom.ia1[exp1] = i0;                             // store base index of monomial in reverse lookup
        // enumerate all monomials of nv2 variables up to order no-no1 in p2
        do {
            DACECom.ie1[i] = exp1;                          // exponents of first nv1 variables
            DACECom.ie2[i] = daceEncodeExponents(p2, no, nv2); // exponents of last nv2 variables
            DACECom.ieo[i] = no1+no2;                       // order of this monomial
            DACECom.ia2[DACECom.ie2[i]] = i-i0;             // store offset of monomial in reverse lookup
            i++;                                            // increase total monomial count
        } while((no2 = daceNextOrderedMonomial(p2, no-no1, nv2)) > 0);
    } while((no1 = daceNextOrderedMonomial(p1, no, nv1)) > 0);

    // free memory
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(p1);
    dacefree(p2);
#endif

    // Crosscheck number of generated monomials. Should never fail.
    if(i != lea) {
        daceSetError(__func__, DACE_PANIC, 5);
        exit(1);
    }

    // check generated adressing arrays ie1,ie2,ia1,ia2. Should never fail.
    for(i = 0; i < lea; i++)
    {
        const unsigned int nn = DACECom.ia1[DACECom.ie1[i]] + DACECom.ia2[DACECom.ie2[i]];
        if(nn != i)
        {
            daceSetError(__func__, DACE_PANIC, 6);
            exit(1);
        }
    }

    // Store various constants in the global data block
    DACECom.nomax = no;
    DACECom.nvmax = nv;
    DACECom.nv1 = nv1;
    DACECom.nv2 = nv2;
    DACECom.nmmax = lea;

    // Set up thread local settings
    daceInitializeThread0();
}

/*! Set up thread local data structures at the beginning of a new thread.
   \note The main thread must call daceInitialize once before spawning new threads.
   All spawned threads must then call daceInitializeThread to initialize the
   thread before performing any other operations.
   \note daceInitialize MUST NOT be called again by any thread while other threads
   are active.
   \note Also initializes the truncation order to the maximum computation order
   and disables the DA epsilon cutoff by setting it to 0.0.
   \sa daceInitialize
   \sa daceTruncationOrder
   \sa daceSetEpsilon
 */
void daceInitializeThread()
{
    daceClearError();
    daceInitializeThread0();
}

/*! Clean up thread local data structures at the end of thread's life time.
 \note Each spawned thread (except for the main thread) should call daceCleanupThread
 before exitting to ensure any thread local memory is properly release. No DACE operations
 must be performed after calling daceCleanupThread.
 \sa daceInitializeThread
 */
void daceCleanupThread()
{
#ifdef DACE_FILTERING
    #if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
        dacefree(DACECom_t.ifi);
    #endif
#endif
}

/*! Set up thread local data structures without resetting error.
    Also initializes the truncation order to the maximum computation order
    and disables the DA epsilon cutoff by setting it to 0.0.
*/
void daceInitializeThread0()
{
    // Set initial cutoff to zero, disabling the optimization
    DACECom_t.eps = 0.0;

    // Set truncation order to full order
    DACECom_t.nocut = DACECom.nomax;

#ifdef DACE_FILTERING
    #if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
        DACE_THREAD_LOCAL static unsigned int ifi[DACE_STATIC_NMMAX];
        DACECom_t.ifi = ifi;
    #else
        dacefree(DACECom_t.ifi);
        DACECom_t.ifi = dacecalloc(DACECom.nmmax, sizeof(unsigned int));
    #endif
    DACECom_t.lfi = 0;
#endif
}

/*! This subroutine returns the major and minor version number of the DACE.
    These values can be checked by the interface to ensure compatibility.
   \param[out] imaj Major version number
   \param[out] imin Minor version number
   \param[out] ipat Patch version number
 */
void daceGetVersion(int *imaj, int *imin, int *ipat)
{
    *imaj = DACE_MAJOR_VERSION;
    *imin = DACE_MINOR_VERSION;
    *ipat = DACE_PATCH_VERSION;
}

/*! Set cutoff value to eps and return the previous value.
   \param[in] eps New cutoff value at or below which coefficients can be flushed to
    zero for efficiency purposes.
   \return The previous value of the cutoff
   \note This feature can have severe unintended consequences if used incorrectly!
    Flushing occurs for any intermediate result also within the DACE, and can result
    in wrong solutions whenever DA coefficients become very small relative to epsilon.
    For example, division by a large DA divisor can cause the (internally calculated)
    multiplicative inverse to be entirely flushed to zero, resulting in a zero DA
    quotient independently of the size of the dividend.
   \sa daceGetEpsilon
 */
double daceSetEpsilon(const double eps)
{
    const double old_eps = DACECom_t.eps;
    DACECom_t.eps = fabs(eps);
    return old_eps;
}

/*! Get the cutoff value eps.
   \return The current value of the cutoff
   \sa daceSetEpsilon
 */
double daceGetEpsilon()
{
    return DACECom_t.eps;
}

/*! Get machine epsilon value.
   \return The experimentally determined machine epsilon
 */
double daceGetMachineEpsilon()
{
    return DACECom.epsmac;
}

/*! Get the maximum computation order set in the initialization routine.
   \return The maximum computation order
   \sa daceInitialize
 */
unsigned int daceGetMaxOrder()
{
    return DACECom.nomax;
}

/*! Get the maximum number of variables set in the initialization routine.
   \return The maximum number of variables
   \sa daceInitialize
 */
unsigned int daceGetMaxVariables()
{
    return DACECom.nvmax;
}

/*! Get the maximum number of monomials for the current setup.
   \return The maximum number of monomials
   \sa daceInitialize
 */
unsigned int daceGetMaxMonomials()
{
    return DACECom.nmmax;
}

/*! Get the current truncation order set for computations.
   \return The current truncation order
   \sa daceSetTruncationOrder
 */
unsigned int daceGetTruncationOrder()
{
    return DACECom_t.nocut;
}

/*! Set the current truncation order for future computations.
   \param[in] fnot The new truncation order
   \return The previous truncation order
   \sa daceGetTruncationOrder
 */
unsigned int daceSetTruncationOrder(const unsigned int fnot)
{
    if(fnot > DACECom.nomax)
    {
        daceSetError(__func__, DACE_INFO, 62);
    }
    const unsigned int old_nocut = DACECom_t.nocut;
    DACECom_t.nocut = umax(umin(fnot, DACECom.nomax), 1);
    return old_nocut;
}
/** @}*/
