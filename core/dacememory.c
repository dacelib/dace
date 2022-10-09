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
 *  dacememory.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/** \addtogroup DACE Core
 *  @{
 */

/********************************************************************************
 *     DACE memory handling routines
 *********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "dace/config.h"
#include "dace/dacebase.h"
#include "dace/daceaux.h"


// Global DACE data structure allocation
dacecom DACECom = { 0 };                        // !< DACE common data block for all threads
DACE_THREAD_LOCAL dacecom_t DACECom_t = { 0 };  // !< DACE common block with local data for each thread
DACE_THREAD_LOCAL dacedbg DACEDbg = { 0 };      // !< DACE common block for error handling


// The DACE has three types of memory allocators
//   * DACE_MEMORY_HYBRID: malloc/free a block of DAs (default)
//   * DACE_MEMORY_STATIC: compile time allocation of a fixed size block of DAs
//   * DACE_MEMORY_DYNAMIC: malloc/free each individual DA

#if DACE_MEMORY_MODEL == DACE_MEMORY_HYBRID || DACE_MEMORY_MODEL == DACE_MEMORY_STATIC

#ifdef WITH_PTHREAD
    #include <pthread.h>

    // mutexes to protect memory management in case of multihreaded execution with pthreads
    static pthread_mutex_t dace_memory_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

// A DACE variable
typedef struct dvariable {
    int len, max;
    unsigned int mem;
} variable;

// Global memory array
typedef struct dmem {
    monomial *mem;
    variable *var;
    unsigned int nda, mda, nst, lda, lst;
} dacemem;

// Global DACE memory data structure and static memory blocks (only local to this file, not visible outside)
static dacemem DACEMem = { 0 };

/*! Reallocate DACE internal memory.
   \param[in] nvar Minimum number of variables the new memory allocation should support
   \param[in] nmem Minimum amount of monomials the new memory allocation should support
*/
void daceReallocateMemory(const unsigned int nvar, const unsigned int nmem)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_HYBRID
    // Determine new lda and lst (special handling for first call, i.e. lda=0)
    unsigned int nlda;
    if(DACEMem.lda == 0)
    {
         nlda = umax(2000, nvar);
    }
    else
    {
         nlda = umax(umin(2000, DACEMem.lda/3), nvar);
    }
    unsigned int nlst = umax(nlda*DACECom.nmmax, nmem);

    DACEMem.lst = DACEMem.lst + nlst;
    DACEMem.lda = DACEMem.lda + nlda;

    // Reallocate new arrays
    DACEMem.mem = dacerealloc(DACEMem.mem, DACEMem.lst*sizeof(monomial));
    DACEMem.var = dacerealloc(DACEMem.var, DACEMem.lda*sizeof(variable));
#elif DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    // can only "allocate" static memory block once, error hard if called again (we're out of memory!)
    if(DACEMem.lst > 0)
    {
        daceSetError(__func__, DACE_PANIC, 2);
        exit(1);
    }

    // static memory blocks
    static monomial dace_static_mem[DACE_STATIC_MEM_SIZE];
    static variable dace_static_var[DACE_STATIC_VAR_SIZE];

    // assign static memory blocks and sizes
    DACEMem.lst = DACE_STATIC_MEM_SIZE;
    DACEMem.lda = DACE_STATIC_VAR_SIZE;
    DACEMem.mem = dace_static_mem;
    DACEMem.var = dace_static_var;
#endif
}

/*! Allocate storage for a DA vector with memory length len.
   \param[out] inc Index of the newly created variable
   \param[in] len Length of the variable to allocate. If len = 0 the length is
    automatically determined to be large enough for any DA vector (i.e. len=nmmax).
 */
void daceAllocateDA(DACEDA *inc, const unsigned int len)
{
    // sanity check to see if daceini has been called before
    if(UNLIKELY(DACECom.nmmax == 0))
    {
        daceSetError(__func__, DACE_PANIC, 3);
        exit(1);
    }

    // check for sane length or autosize
    unsigned int ilen = len;
    if(LIKELY(ilen == 0))
        ilen = DACECom.nmmax;

#ifdef WITH_PTHREAD
    pthread_mutex_lock(&dace_memory_mutex);
#endif
    // Update mda to point to the lowest free variable
    while((DACEMem.mda < DACEMem.nda) && (DACEMem.var[DACEMem.mda].max > 0))
         DACEMem.mda++;

     // Find the next free variable of appropriate size
    for( unsigned int i = DACEMem.mda; i < DACEMem.nda; i++ )
    {
        if(DACEMem.var[i].max <= -((int)ilen))
        {
            *inc = i;
            DACEMem.var[i].max *= -1;
            DACEMem.var[i].len = 0;
            if(i == DACEMem.mda) DACEMem.mda++;
#ifdef WITH_PTHREAD
            pthread_mutex_unlock(&dace_memory_mutex);
#endif
            return;
        }
    }

    // No appropriate variable found, allocate new one

    // Reallocate more memory if we ran out of variables or stack
    if(UNLIKELY(((DACEMem.nda+1) > DACEMem.lda) || ((DACEMem.nst+ilen) > DACEMem.lst))) daceReallocateMemory(1, ilen);

    // Allocate new variable
    *inc = DACEMem.nda;
    DACEMem.var[DACEMem.nda].mem = DACEMem.nst;
    DACEMem.var[DACEMem.nda].max = ilen;
    DACEMem.var[DACEMem.nda].len = 0;
    DACEMem.nst += ilen;
    if(DACEMem.mda == DACEMem.nda) DACEMem.mda++;
    DACEMem.nda++;
#ifdef WITH_PTHREAD
    pthread_mutex_unlock(&dace_memory_mutex);
#endif
}

/*! Deallocate DA vector inc.
   \param[in] inc Index of the DA variable to free
 */
void daceFreeDA(DACEDA *inc)
{
    // sanity check to see if daceini has been called before
    if(UNLIKELY(DACECom.nmmax == 0))
    {
        daceSetError(__func__, DACE_PANIC, 3);
        exit(1);
    }
    else if(UNLIKELY(*inc == -1))
    {
        // already free or otherwise reallocated
        return;
    }

#ifdef WITH_PTHREAD
    pthread_mutex_lock(&dace_memory_mutex);
#endif
    // Check for sane argument
    if(UNLIKELY((*inc < 0) || (*inc >= (int)DACEMem.nda) || (DACEMem.var[*inc].max <= 0)))
    {
        daceSetError(__func__, DACE_INFO, 61);
#ifdef WITH_PTHREAD
        pthread_mutex_unlock(&dace_memory_mutex);
#endif
        return;
    }

    // Mark variable as free and update minimum free variable
    DACEMem.var[*inc].max *= -1;
    if(*inc < (int)DACEMem.mda)
        DACEMem.mda = *inc;

    // fully release variable(s) on top of stack
    if(*inc == DACEMem.nda-1)
    {
         DACEMem.nda--;
         while((DACEMem.nda > 0) && (DACEMem.var[DACEMem.nda-1].max <= 0))
            DACEMem.nda--;
         if(UNLIKELY(DACEMem.nda == 0))
         {
            DACEMem.nst = 0;
         }
         else
         {
            DACEMem.nst = DACEMem.var[DACEMem.nda-1].mem + DACEMem.var[DACEMem.nda-1].max;
         }
    }

    *inc = -1;
#ifdef WITH_PTHREAD
    pthread_mutex_unlock(&dace_memory_mutex);
#endif
}

/*! Invalidate DA vector inc without deallocating associated memory.
   \param[in] inc Index of the DA variable to invalidate
 */
void daceInvalidateDA(DACEDA *inc)
{
    *inc = -1;
}

/*! Dump information about the current memory management status to stdout.
 */
void daceMemoryDump()
{
    printf("DACE core memory model: " DACE_MEMORY_MODEL_STR "\n");
    printf("Allocated system memory:\n");
    double m = DACEMem.lda*(4.0*3.0)/1024.0/1024.0;
    printf("\t# of DAs = %d \t(%d MB)\n", DACEMem.lda, (unsigned int)m);
    m = DACEMem.lst*(4.0*2.0+8.0)/1024.0/1024.0;
    printf("\tstack    = %d \t(%d MB)\n", DACEMem.lst, (unsigned int)m);

    printf("Allocated internal memory:\n");
    m = (double)DACEMem.nda/DACEMem.lda*100.0;
    printf("\t# of DAs = %d \t(%d %%)\n", DACEMem.nda, (unsigned int)m);
    m = (double)DACEMem.nst/DACEMem.lst*100.0;
    printf("\tstack    = %d \t(%d %%)\n", DACEMem.nst, (unsigned int)m);

    if(DACEMem.lda > 0)
    {
         printf("Objects:\n");
         unsigned int j = 0, l = 0;
         for( unsigned int i = 0; i < DACEMem.nda; i++ )
         {
            if(DACEMem.var[i].max > 0)
            {
               j++;
               l += DACEMem.var[i].max;
            }
         }
         m = (double)j/umax(1, DACEMem.nda)*100.0;
         printf("\tused objects     = %d \t(%d %%)\n", j, (unsigned int)m);
         m = l*(4.0*2.0+8.0)/1024.0/1024.0;
         printf("\tused stack       = %d \t(%d MB)\n", l, (unsigned int)m);
         m = (double)(DACEMem.nda-j)/umax(1, DACEMem.nda)*100.0;
         printf("\tunused objects   = %d \t(%d %%)\n", DACEMem.nda-j, (unsigned int)m);
    }
}

/*! Extract internal information about a DA object.
   \param[in] inc Pointer to the DA object to extract information from
   \param[out] ipoc Pointer to an array of monomials allocated for this variable
   \param[out] ilmc Pointer where to store the maximum number of monomials allocated in this DA object
   \param[out] illc Pointer where to store the currently used length of this DA object
*/
void daceVariableInformation(const DACEDA *inc, monomial **ipoc, unsigned int *ilmc, unsigned int *illc)
{
    if(UNLIKELY(*inc >= 0 && *inc < (int)DACEMem.nda))
    {
        *ipoc = DACEMem.mem+DACEMem.var[*inc].mem;
        *ilmc = DACEMem.var[*inc].max;
        *illc = DACEMem.var[*inc].len;
        if(*ilmc > 0) return;
    }

    *ipoc = NULL;
    *ilmc = 0;
    *illc = 0;
    daceSetError(__func__, DACE_PANIC, 4);
    exit(1);
}

/*! Set the length of a DACE DA object.
   \param[in] inc The DACE DA object to operate on
   \param[in] len The new length of the object
*/
void daceSetLength(DACEDA *inc, const size_t len)
{
    if(UNLIKELY(DACEMem.var[*inc].max < len))
    {
        daceSetError(__func__, DACE_PANIC, 7);
        exit(1);
        // catastrophic error because we may have written past the end of the variable and contaminated memory there
    }
    DACEMem.var[*inc].len = (int)len;
}

/*! Compare if two DACE DA objects refer to the same underlying memory (i.e. are the same object).
   \param[in] ina The first DACE DA object
   \param[in] inb The second DACE DA object
   \return returns true if the two objects are the same, false otherwise
*/
bool daceIsSameObject(const DACEDA *ina, const DACEDA *inb)
{
    return *ina == *inb;    // compare the variable numbers
}

/*! Free the entire DACE memory and purge all DA objects that may have been
    allocated before.
*/
void daceFreeMemory()
{
    DACEMem.nda = 0;
    DACEMem.mda = 1;
    DACEMem.nst = 0;
    DACEMem.lda = 0;
    DACEMem.lst = 0;

#if DACE_MEMORY_MODEL == DACE_MEMORY_HYBRID
    // free memory (will be reallocated in next call to daceall)
    dacefree(DACEMem.mem);
    DACEMem.mem = NULL;
    dacefree(DACEMem.var);
    DACEMem.var = NULL;
#endif      // DACE_MEMORY_MODEL == DACE_MEMORY_DYNAMIC
}


#elif DACE_MEMORY_MODEL == DACE_MEMORY_DYNAMIC


/*! Allocate storage for a DA vector with memory length len.
   \param[out] inc Index of the newly created variable
   \param[in] len Length of the variable to allocate. If len = 0 the length is
    automatically determined to be large enough for any DA vector (i.e. len=nmmax).
 */
void daceAllocateDA(DACEDA *inc, const unsigned int len)
{
    // sanity check to see if daceini has been called before
    if(UNLIKELY(DACECom.nmmax == 0))
    {
        daceSetError(__func__, DACE_PANIC, 3);
        exit(1);
    }

    // check for sane length or autosize
    unsigned int ilen = len;
    if(LIKELY(ilen == 0))
        ilen = DACECom.nmmax;

    // just allocate memory dynamically
    inc->len = 0;
    inc->max = ilen;
    inc->mem = (monomial*)dacemalloc(ilen*sizeof(monomial));
}

/*! Deallocate DA vector inc.
   \param[in] inc Index of the DA variable to free
 */
void daceFreeDA(DACEDA *inc)
{
    dacefree(inc->mem);
    inc->len = 0;
    inc->mem = NULL;
}

/*! Invalidate DA vector inc without deallocating associated memory.
   \param[in] inc Index of the DA variable to invalidate
 */
void daceInvalidateDA(DACEDA *inc)
{
    inc->len = 0;
    inc->mem = NULL;
}

/*! Dump information about the current memory management status to stdout.
 */
void daceMemoryDump()
{
    printf("DACE core memory model: " DACE_MEMORY_MODEL_STR "\n");
    return;
}

/*! Extract internal information about a DA object.
   \param[in] inc Pointer to the DA object to extract information from
   \param[out] ipoc Pointer to an array of monomials allocated for this variable
   \param[out] ilmc Pointer where to store the maximum number of monomials allocated in this DA object
   \param[out] illc Pointer where to store the currently used length of this DA object
*/
void daceVariableInformation(const DACEDA *inc, monomial **ipoc, unsigned int *ilmc, unsigned int *illc)
{
    *ipoc = inc->mem;
    *ilmc = inc->max;
    *illc = inc->len;

    if(UNLIKELY(inc->mem == NULL))
    {
        daceSetError(__func__, DACE_PANIC, 4);
        exit(1);
    }
}

/*! Set the length of a DACE DA object.
   \param[in] inc The DACE DA object to operate on
   \param[in] len The new length of the object
*/
void daceSetLength(DACEDA *inc, const size_t len)
{
    if(UNLIKELY(inc->max < len))
    {
        daceSetError(__func__, DACE_PANIC, 7);
        exit(1);
        // catastrophic error because we may have written past the end of the variable and contaminated memory there
    }
    inc->len = (int)len;
}

/*! Compare if two DACE DA objects refer to the same underlying memory (i.e. are the same object).
   \param[in] ina The first DACE DA object
   \param[in] inb The second DACE DA object
   \return returns true if the two objects are the same, false otherwise
*/
bool daceIsSameObject(const DACEDA *ina, const DACEDA *inb)
{
    return ina->mem == inb->mem;        // compare memory addresses
}

/*! Free the entire DACE memory and purge all DA objects that may have been
    allocated before.
*/
void daceFreeMemory()
{
    return;
}

#else
#error Invalid DACE memory model selected!
#endif // DACE_MEMORY_MODEL
/** @}*/
