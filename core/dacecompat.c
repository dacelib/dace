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
 *  dacecompat.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/** \addtogroup DACE Core
 *  @{
 */

#include "dace/config.h"
#include "dace/dacebase.h"
#include "dace/dacecompat.h"
#include "dace/daceaux.h"

/********************************************************************************
 *     Compatibility shims to bridge old and new interface
 *********************************************************************************/

/*! Return the number of non-zero monomials in a DA object.
   \param[in] ina Pointer to DA object to get length of
   \param[out] size Number of non-zero monomials
   \deprecated Has been renamed to daceGetLength()
   \sa daceGetLength()
*/
void dacesize(const DACEDA *ina, unsigned int *size)
{
    *size = daceGetLength(ina);
}

/*! Extract coefficient of a monomial in a DA object.
   \param[in] ina Pointer to DA object to extract monomial coefficient from
   \param[in] jj C array of nvmax exponents identifying the monomial
   \param[out] cjj Pointer where to store the value of the coefficient
   \deprecated Has been replaced by daceGetCoefficient()
   \sa daceGetCoefficient()
*/
void dacepek(const DACEDA *ina, const unsigned int jj[], double *cjj)
{
    *cjj = daceGetCoefficient0(ina, daceEncode(jj));
}

/*! Compute absolute value of a DA object.
    Same as daceNorm(ina, 0).
   \param[in] ina Pointer to DA object to take absolute value of
   \param[out] anorm Pointer where to store the absolute value
   \deprecated Has been replaced by daceNorm() and daceAbsolute()
   \sa daceNorm()
   \sa daceAbsolute()
*/
void daceabs(const DACEDA *ina, double *anorm)
{
    *anorm = daceNorm(ina, 0);
}

/*! Compute absolute value of a DA object.
   \param[in] ina Pointer to DA object to take norm of
   \param[in] ityp Type of norm to take (see daceNorm())
   \param[out] anorm Pointer where to store the norm
   \deprecated Has been replaced by daceNorm()
   \sa daceNorm()
*/
void dacenorm(const DACEDA *ina, const unsigned int ityp, double *anorm)
{
    *anorm = daceNorm(ina, ityp);
}

/*! Compute an evaluation tree for efficient evaluation of DA objects.
   \param[in] das Pointers to DA objects to evaluate
   \param[in] count Number of DA objects in das[]
   \param[out] ac C array containing compiled evaluation tree data
   \param[out] nterm resulting number of terms in evaluation tree
   \param[out] nvar number of variables appearing in evaluation tree
   \param[out] nord maximum order appearing in evaluation tree
   \deprecated Has been replaced by daceEvalTree()
   \sa daceEvalTree()
*/
void dacetree(const DACEDA das[], const unsigned int count, double ac[],  unsigned int *nterm, unsigned int *nvar, unsigned int *nord)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    if(count > 10*DACE_STATIC_NVMAX)
    {
        daceSetError("dacetree: number of DAs too high", DACE_ERROR, 111);
        return;
    }
    const DACEDA *temp[10*DACE_STATIC_NVMAX];
#else
    const DACEDA **temp = dacecalloc(count, sizeof(const DACEDA *));
#endif
    for(unsigned int i = 0; i < count; i++) temp[i] = &das[i];
    daceEvalTree(temp, count, ac, nterm, nvar, nord);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree((void*)temp);
#endif
}
/** @}*/
