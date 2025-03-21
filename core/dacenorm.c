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
 *  dacenorm.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/** \addtogroup DACE Core
 *  @{
 */

#include <math.h>
#include <stdbool.h>

#include "dace/config.h"
#include "dace/dacebase.h"
#include "dace/daceaux.h"


/********************************************************************************
 *     DACE norm and norm estimation routines
 *********************************************************************************/

/*! Compute the absolute value (maximum coefficient norm) of a DA object.
   \param[in] ina Pointer to the DA object to take absolute value of
   \return The absolute value of ina
*/
double daceAbsoluteValue(const DACEDA *ina)
{
    return daceNorm(ina, 0);
}

/*! Compute a norm of a DA object.
   \param[in] ina Pointer to the DA object to take norm of
   \param[in] ityp Type of norm to compute.
     0 = max norm
     1 = sum norm
    >1 = corresponding vector norm
   \return The norm of ina
*/
double daceNorm(const DACEDA *ina, const unsigned int ityp)
{
    monomial *ipoa; unsigned int ilma, illa;
    double anorm = 0.0;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    if(ityp == 0)
    {
        // max norm
        for(monomial *i = ipoa; i < ipoa+illa; i++)
            anorm = fmax(anorm, fabs(i->cc));
    }
    else if(ityp == 1)
    {
        // sum norm
        for(monomial *i = ipoa; i < ipoa+illa; i++)
            anorm += fabs(i->cc);
    }
    else
    {
        // ityp vector norm
        for(monomial *i = ipoa; i < ipoa+illa; i++)
            anorm += pown(fabs(i->cc), ityp);
        anorm = pow(anorm, 1.0/ityp);
    }
    return anorm;
}

/*! Compute an order sorted norm of a DA object.
   \param[in] ina Pointer to the DA object to take norm of
   \param[in] ivar Independent variable with respect to which to group.
     0 = group by monomial order
    >1 = group by given independent variable
   \param[in] ityp Type of norm to compute.
     0 = max norm
     1 = sum norm
    >1 = corresponding vector norm
   \param[out] onorm C array of length nomax+1 containing the grouped estimates
*/
void daceOrderedNorm(const DACEDA *ina, const unsigned int ivar, const unsigned int ityp, double onorm[])
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    for(unsigned int i = 0; i <= DACECom.nomax; i++) onorm[i] = 0.0;

    if(ivar > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        return;
    }

    if(ivar == 0)
    {
        if(ityp == 0)
        {
            // max norm
            for(monomial *i = ipoa; i < ipoa+illa; i++ )
            {
                const unsigned int io = DACECom.ieo[i->ii];
                onorm[io] = fmax(onorm[io], fabs(i->cc));
            }
        }
        else if(ityp == 1)
        {
            // sum norm
            for(monomial *i = ipoa; i < ipoa+illa; i++ )
            {
                const unsigned int io = DACECom.ieo[i->ii];
                onorm[io] += fabs(i->cc);
            }
        }
        else
        {
            // ityp vector norm
            for(monomial *i = ipoa; i < ipoa+illa; i++ )
            {
                const unsigned int io = DACECom.ieo[i->ii];
                onorm[io] += pown(fabs(i->cc), ityp);
            }
            for(unsigned int i = 0; i <= DACECom.nomax; i++)
                onorm[i] = pow(onorm[i], 1.0/ityp);
        }
    }
    else
    {
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
        unsigned int jj[DACE_STATIC_NVMAX];
#else
        unsigned int *jj = (unsigned int*)dacecalloc(DACECom.nvmax, sizeof(unsigned int));
#endif

        if(ityp == 0)
        {
            // max norm
            for(monomial *i = ipoa; i < ipoa+illa; i++ )
            {
                daceDecode(i->ii, jj);
                onorm[jj[ivar-1]] = fmax(onorm[jj[ivar-1]], fabs(i->cc));
            }
        }
        else if(ityp == 1)
        {
            // sum norm
            for(monomial *i = ipoa; i < ipoa+illa; i++ )
            {
                daceDecode(i->ii, jj);
                onorm[jj[ivar-1]] += fabs(i->cc);
            }
        }
        else
        {
            // ityp vector norm
            for(monomial *i = ipoa; i < ipoa+illa; i++ )
            {
                daceDecode(i->ii, jj);
                onorm[jj[ivar-1]] += pown(fabs(i->cc), ityp);
            }
            for(unsigned int i = 0; i <= DACECom.nomax; i++)
                onorm[i] = pow(onorm[i], 1.0/ityp);
        }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
        dacefree(jj);
#endif
    }
}

/*! Estimate order sorted norms of DA object ina up to given order.
   \param[in] ina Pointer to the DA object to take norm of
   \param[in] ivar Independent variable with respect to which to group.
     0 = group by monomial order
    >1 = group by given independent variable
   \param[in] ityp Type of norm to compute.
     0 = max norm
     1 = sum norm
    >1 = corresponding vector norm
   \param[in] nc Maximum order to estimate
   \param[out] c C array of length nc+1 containing the grouped estimates
   \param[out] err C array of length min(nc, nomax)+1 containing the residuals
    of the exponential fit at each order. If NULL is passed in, no residuals
    are computed and returned.
   \note If estimation is not possible, zero is returned for all
    requested orders. In most cases this is actually not too far off.
*/
void daceEstimate(const DACEDA *ina, const unsigned int ivar, const unsigned int ityp, double c[], double err[], const unsigned int nc)
{
    for(unsigned int i = 0; i <= nc; i++) c[i] = 0.0;

    // check input
    if(DACECom.nomax < 2)
    {
        daceSetError(__func__, DACE_ERROR, 51);
        return;
    }

    // get order sorted norms
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double onorm[DACE_STATIC_NOMAX+1];
#else
    double *onorm = (double*)dacecalloc(DACECom.nomax+1, sizeof(double));
#endif
    daceOrderedNorm(ina, ivar, ityp, onorm);

    // set up xtx^-1 and xty for linear least squares
    double ai[2] = {0.0};
    double xtx[2][2] = {{0.0}, {0.0}};
    for(unsigned int i = 1; i <= DACECom.nomax; i++)
    {
        if(!(onorm[i] <= DACECom_t.eps))
        {
            xtx[0][0] += i*i;
            xtx[0][1] -= i;
            xtx[1][1]++;
            ai[0] += log(onorm[i]);
            ai[1] += i*log(onorm[i]);
        }
    }

    if(xtx[1][1] < 2)
    {
        daceSetError(__func__, DACE_INFO, 63);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
        dacefree(onorm);
#endif
        return;
    }

    xtx[1][0] = xtx[0][1];
    const double det = xtx[0][0]*xtx[1][1]-xtx[0][1]*xtx[1][0];

    // compute ai=(xtx^-1)(xty)
    double a[2];
    a[0] = (ai[0]*xtx[0][0]+ai[1]*xtx[0][1])/det;
    a[1] = (ai[0]*xtx[1][0]+ai[1]*xtx[1][1])/det;

    // actually compute the order sorted norm estimates
    for(unsigned int i = 0; i <= nc; i++)
        c[i] = exp(a[0]+a[1]*i);

    // compute estimation error, if requested
    if(err != NULL)
        for(unsigned int i = 0; i <= umin(DACECom.nomax, nc); i++)
        {
            const double temp = onorm[i]-c[i];
            err[i] = temp > 0.0 ? temp : 0.0;
        }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(onorm);
#endif
}

/*! Compute an upper and lower bound of DA object ina over the domain [-1,1]^n.
   \param[in] ina Pointer to the DA object to bound
   \param[out] alo Pointer where to store the lower bound
   \param[out] aup Pointer where to store the upper bound
*/
void daceGetBounds(const DACEDA *ina, double *alo, double *aup)
{
    monomial *ipoa; unsigned int ilma, illa;
    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    // check special cases
    *alo = 0.0;
    *aup = 0.0;
    if(illa == 0) return;

    // constant part is special
    if(ipoa->ii == 0)
    {
        *alo = ipoa->cc;
        *aup = *alo;
        ipoa++;
        illa--;
    }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    unsigned int jj[DACE_STATIC_NVMAX];
#else
    unsigned int *jj = (unsigned int*)dacecalloc(DACECom.nvmax, sizeof(unsigned int));
#endif

    // bound each non-constant term
    for(monomial *i = ipoa; i < ipoa+illa; i++)
    {
        daceDecode(i->ii, jj);
        bool odd = false;
        for(unsigned int j = 0; j < DACECom.nvmax; j++)
        {
            if(jj[j] & 1)
            {
                odd = true;
                break;
            }
        }

        if(odd)
        {
            *aup += fabs(i->cc);
            *alo -= fabs(i->cc);
        }
        else
        {
            if(i->cc > 0.0)
               *aup += i->cc;
            else
               *alo += i->cc;
        }
    }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(jj);
#endif
}
/** @}*/
