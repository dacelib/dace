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
 *  dacebasic.c
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


/********************************************************************************
 *     DACE variable creation routines
 *********************************************************************************/

/*! Create a DA object to be ckon times the i-th independen variable.
   \param[in] ina Pointer to DA object to store the resulting DA in
   \param[in] i number of the independent variable to create
   \param[in] ckon coefficient of the independent DA variable created
   \note Independent DA variable indices are 1-based, i.e. the first independent
    variable is i=1. The case of i=0 corresponds to the constant part of the polynomial.
*/
void daceCreateVariable(DACEDA *ina, const unsigned int i, const double ckon)
{
    monomial *ipoa; unsigned int ilma, illa;
    unsigned int ic1, ic2, base;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    daceSetLength(ina, 0);
    if(i > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        return;
    }

    if(fabs(ckon) <= DACECom_t.eps)
    {
        return;
    }

    /* check the length, although it should never be zero in any case */
    if(ilma < 1)
    {
        daceSetError(__func__, DACE_ERROR, 21);
        return;
    }

    /* set up the exponents */
    ic1 = 0;
    ic2 = 0;

    if(i != 0)
    {
        base = DACECom.nomax+1;
        if(i > DACECom.nv1)
        {
            ic2 = npown(base, i-1-DACECom.nv1);
        }
        else
        {
            ic1 = npown(base, i-1);
        }
    }

    daceSetLength(ina, 1);
    ipoa->cc = ckon;
    ipoa->ii = DACECom.ia1[ic1]+DACECom.ia2[ic2];
}

/*! Create a DA object to be ckon times the monomial given by the exponents in jj[].
   \param[in] ina Pointer to DA object to store the resulting DA in
   \param[in] jj C array with nvmax exponents indicating the monomial to create
   \param[in] ckon coefficient of the monomial created
*/
void daceCreateMonomial(DACEDA *ina, const unsigned int jj[], const double ckon)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    if(ilma < 1)
    {
        daceSetError(__func__, DACE_ERROR, 21);
        daceSetLength(ina, 0);
        return;
    }

    if(fabs(ckon) <= DACECom_t.eps)
    {
        daceSetLength(ina, 0);
    }
    else
    {
        ipoa->ii = daceEncode(jj);
        ipoa->cc = ckon;
        daceSetLength(ina, 1);
    }
}

/*! Create a DA object with all coefficients set to the constant value ckon.
   \param[in] ina Pointer to DA object to store the resulting DA in
   \param[in] ckon coefficient of the monomials
*/
void daceCreateFilled(DACEDA *ina, const double ckon)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    unsigned int i = 0;
    for(monomial *ia = ipoa; ia < ipoa+ilma && i < DACECom.nmmax; ia++, i++)
    {
        ia->ii = i;
        ia->cc = ckon;
    }

    daceSetLength(ina, i);
}

/*! Create a DA object with constant part equal to ckon.
   \param[in] ina Pointer to DA object to store the resulting DA in
   \param[in] ckon coefficient of the constant part of the result
*/
void daceCreateConstant(DACEDA *ina, const double ckon)
{
    daceCreateVariable(ina, 0, ckon);
}

/*! Create a DA object with randomly filled coefficients.
   \param[in] ina Pointer to DA object to store the resulting DA in
   \param[in] cm The filling factor between -1.0 and 1.0.
    The absolute value of the filling factor determines the fraction of non-zero
    coefficients.
    If cm is positive, the values are weighted by order such that
    the coefficients decay exponentially with the order from 1.0 towards the
    machine epsilon in the highest order.
    If cm is negative, all coefficients are chosen to be between -1.0 and 1.0.
*/
void daceCreateRandom(DACEDA *ina, const double cm)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    monomial *ia = ipoa, *iamax = ipoa+ilma;
    if(cm < 0.0)
    {
        for(unsigned int i = 0; i < DACECom.nmmax && ia < iamax; i++)
        {
            if((DACECom.ieo[i] <= DACECom_t.nocut) && (daceRandom() < -cm))
            {
                ia->cc = 2.0*daceRandom()-1.0;
                ia->ii = i;
                ia++;
            }
        }
    }
    else
    {
        for(unsigned int i = 0; i < DACECom.nmmax && ia < iamax; i++)
        {
            if((DACECom.ieo[i] <= DACECom_t.nocut) && (daceRandom() < cm))
            {
                const double w = pow(DACECom.epsmac, (DACECom.ieo[i]/((double)DACECom_t.nocut)));
                ia->cc = w*(2.0*daceRandom()-1.0);
                ia->ii = i;
                ia++;
            }
        }
    }

    daceSetLength(ina, ia-ipoa);
}

/********************************************************************************
 *     DACE coefficient access routines
 *********************************************************************************/

/*! Extract the constant part from a DA object.
   \param[in] ina Pointer to DA object to extract constant part from
   \return Constant part of the given DA object
*/
double daceGetConstant(const DACEDA *ina)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    if(illa > 0 && ipoa->ii == 0)
    {
        return ipoa->cc;
    }
    else
    {
        return 0.0;
    }
}

/*! Extract the linear part of a DA object.
   \param[in] ina Pointer to DA object to extract linear part from
   \param[in] c C array of length nvmax containing the linear coefficients in order
*/
void daceGetLinear(const DACEDA *ina, double c[])
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    unsigned int jj[DACE_STATIC_NVMAX] = {0};
#else
    unsigned int *jj = (unsigned int*) dacecalloc(DACECom.nvmax, sizeof(unsigned int));
#endif

    for(unsigned int i = 0; i < DACECom.nvmax; i++)
    {
        jj[i] = 1;
        c[i] = daceGetCoefficient(ina, jj);
        jj[i] = 0;
    }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(jj);
#endif
}

/*! Extract coefficient of a monomial in a DA object.
   \param[in] ina Pointer to DA object to extract monomial coefficient from
   \param[in] jj C array of nvmax exponents identifying the monomial
   \return The coefficient of the given monomial in the DA object
*/
double daceGetCoefficient(const DACEDA *ina, const unsigned int jj[])
{
    return daceGetCoefficient0(ina, daceEncode(jj));
}

/*! Extract coefficient of a monomial in a DA object.
   \param[in] ina Pointer to DA object to extract monomial coefficient from
   \param[in] ic DA coding integer of the monomial to extract
   \return The coefficient of the given monomial in the DA object
*/
double daceGetCoefficient0(const DACEDA *ina, const unsigned int ic)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    /* Check for zero length vector */
    if(illa == 0) return 0.0;

    /* Determine if monomial is inside first and last monomials of A */
    monomial *iu = ipoa;
    monomial *iz = ipoa+(illa-1);
    if(ic == iu->ii)
    {
        return iu->cc;
    }
    else if(ic == iz->ii)
    {
        return iz->cc;
    }
    else if(ic < iu->ii || ic > iz->ii)
    {
        return 0.0;
    }

    /* binary search for proper monomial */
    while(iz-iu > 1)
    {
        monomial *i = iu+((iz-iu)/2);
        if(i->ii < ic)
        {
            iu = i;
        }
        else if(i->ii > ic)
        {
            iz = i;
        }
        else
        {
            return i->cc;
        }
    }

    return 0.0;
}

/*! Set coefficient of a monomial in a DA object.
   \param[in] ina Pointer to DA object to set monomial in
   \param[in] jj C array of nvmax exponents identifying the monomial
   \param[in] cjj Value of the corresponding coefficient
*/
void daceSetCoefficient(DACEDA *ina, const unsigned int jj[], const double cjj)
{
    daceSetCoefficient0(ina, daceEncode(jj), cjj);
}

/*! Set coefficient of a monomial in a DA object.
   \param[in] ina Pointer to DA object to set monomial in
   \param[in] ic DA coding integer of the monomial to set
   \param[in] cjj Value of the corresponding coefficient
*/
void daceSetCoefficient0(DACEDA *ina, const unsigned int ic, const double cjj)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *i;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    bool insert = true;

    if(illa == 0)
    {
        i = ipoa;
    }
    else
    {
        monomial *iu = ipoa;
        monomial *iz = ipoa+(illa-1);
        if(ic == iu->ii)
        {
            i = iu;
            insert = false;
        }
        else if(ic < iu->ii)
        {
            i = iu;
        }
        else if(ic == iz->ii)
        {
            i = iz;
            insert = false;
        }
        else if(ic > iz->ii)
        {
            i = iz+1;
        }
        else
        {
            while(iz-iu > 1)
            {
                i = iu+((iz-iu)/2);

                if(i->ii < ic)
                {
                    iu = i;
                }
                else if(i->ii > ic)
                {
                    iz = i;
                }
                else
                {
                    insert = false;
                    iz = i;
                    break;
                }
            }
            i = iz;
        }
    }

    if(insert)
    {
        // insert a new monomial before the i-th monomial
        if(fabs(cjj) <= DACECom_t.eps) return;

        if(illa+1 > ilma)
        {
            daceSetError(__func__, DACE_ERROR, 21);
            return;
        }

        //memmove(i+1, i, (ipoa+illa - i)*sizeof(monomial));
        for(monomial *ii = ipoa+illa; ii > i; ii--)
            *ii = *(ii-1);

        i->cc = cjj;
        i->ii = ic;

        daceSetLength(ina, illa+1);
    }
    else
    {
        // replace the i-th monomial
        if(!(fabs(cjj) <= DACECom_t.eps))
        {
            i->cc = cjj;
        }
        else
        {
            //memmove(i, i+1, (ipoa+illa - i-1)*sizeof(monomial));
            for(monomial *ii = i; ii < ipoa+(illa-1); ii++)
                *ii = *(ii+1);
            daceSetLength(ina, illa-1);
        }
    }
}

/*! Extract coefficient at position npos (starting with 1) in the list of
    non-zero coefficients in the DA object and return its exponents and
    coefficient. If the monomial does not exist, the value 0.0 is returned.
   \param[in] ina Pointer to DA object to extract monomial from
   \param[in] npos Index of the monomial to extract
   \param[out] jj C array of nvmax elements for returning the exponents of the monomial
   \param[out] cjj Pointer where to store the value of the coefficient of the monomial
*/
void daceGetCoefficientAt(const DACEDA *ina, const unsigned int npos, unsigned int jj[], double *cjj)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    if(npos > 0 && npos <= illa)
    {
        *cjj = ipoa[npos-1].cc;
        daceDecode(ipoa[npos-1].ii, jj);
    }
    else
    {
        *cjj = 0.0;
        for(unsigned int j = 0; j < DACECom.nvmax; j++) jj[j] = 0;
    }
}

/*! Return the number of non-zero monomials in a DA object.
   \param[in] ina Pointer to DA object to get length of
   \return Number of non-zero monomials
*/
unsigned int daceGetLength(const DACEDA *ina)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    return illa;
}

/*! Copy content of one DA object into another DA object.
   \param[in] ina Pointer to DA object to copy from
   \param[in] inb Pointer to DA object to copy to
*/
void daceCopy(const DACEDA *ina, DACEDA *inb)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipob; unsigned int ilmb, illb;

    if(daceIsSameObject(ina, inb)) return;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inb, &ipob, &ilmb, &illb);

    if(illa > ilmb)
    {
        daceSetError(__func__, DACE_ERROR, 21);
        illa = ilmb;
    }

    memmove(ipob, ipoa, illa*sizeof(monomial));

    daceSetLength(inb, illa);
}

/*! Copy content of one DA object into another DA object filtering out terms
    below a certain threshold.
   \param[in] ina Pointer to DA object to copy from
   \param[in] inb Pointer to DA object to copy to
   \note This routine is slightly worse than non-filtering version (about 10%)
*/
void daceCopyFiltering(const DACEDA *ina, DACEDA *inb)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipob; unsigned int ilmb, illb;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inb, &ipob, &ilmb, &illb);

    monomial *ib = ipob;
    if(illa <= ilmb)
    {
        for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
        {
            if(fabs(ia->cc) <= DACECom_t.eps || DACECom.ieo[ia->ii] > DACECom_t.nocut)
                continue;
            *ib = *ia;
            ib++;
        }
    }
    else
    {
        for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
        {
            if(fabs(ia->cc) <= DACECom_t.eps || DACECom.ieo[ia->ii] > DACECom_t.nocut)
                continue;
            if(ib >= ipob+ilmb)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            *ib = *ia;
            ib++;
        }
    }

    daceSetLength(inb, ib-ipob);
}

/*! Truncate a DA object to contain only terms of order larger or equal to imin
    and less than or equal imax.
   \param[in] ina Pointer to DA object to trim
   \param[in] imin Minimum order to keep
   \param[in] imax Maximum order to keep
   \param[in] inc Pointer to DA object to store the truncated result in
*/
void daceTrim(const DACEDA *ina, const unsigned int imin, const unsigned int imax, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    monomial *const icmax = ipoc+ilmc;
    monomial *ic = ipoc;
    for(monomial *i = ipoa; i < ipoa+illa; i++)
    {
        const unsigned int io = DACECom.ieo[i->ii];
        if((io > imax)||(io < imin))
        {
            continue;
        }
        if(ic >= icmax)
        {
            daceSetError(__func__, DACE_ERROR, 21);
            break;
        }
        *ic = *i;
        ic++;
    }

    daceSetLength(inc, ic-ipoc);
}

/*! Copy monomials from a DA object ina to DA object inb if the same monomial
    is non-zero in DA object inc, while filtering out terms below the current
    cutoff.
   \param[in] ina Pointer to DA object to filter
   \param[in] inb Pointer to DA object to store the filtered result in
   \param[in] inc Pointer to DA object providing the filter template
*/
void daceFilter(const DACEDA *ina, DACEDA *inb, const DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipob; unsigned int ilmb, illb;
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inb, &ipob, &ilmb, &illb);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    monomial *ib = ipob;
    monomial *ic = ipoc;
    if(illa <= ilmb)
    {
        for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
        {
            // skip forward in inc and exit if it's over
            while(ic < ipoc+illc && ia->ii > ic->ii)
                ic++;
            if(ic >= ipoc+illc) break;

            if(ia->ii < ic->ii || fabs(ia->cc) <= DACECom_t.eps || DACECom.ieo[ia->ii] > DACECom_t.nocut)
                continue;

            *ib = *ia;
            ib++;
        }
    }
    else
    {
        for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
        {
            // skip forward in inc and exit if it's over
            while(ic < ipoc+illc && ia->ii > ic->ii)
                ic++;
            if(ic >= ipoc+illc) break;

            if(ia->ii < ic->ii || fabs(ia->cc) <= DACECom_t.eps || DACECom.ieo[ia->ii] > DACECom_t.nocut)
                continue;

            if(ib >= ipob+ilmb)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }

            *ib = *ia;
            ib++;
        }
    }

    daceSetLength(inb, ib-ipob);
}

/*! Check each coefficient of DA object ina to see if any of them are NANs (not a number).
   \param[in] ina Pointer to DA object to check
   \return True (non-zero) if any of the coefficients of ina is NAN
*/
unsigned int daceIsNan(const DACEDA *ina)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
    {
        if(isnan(ia->cc))
            return true;
    }
    return false;
}

/*! Check each coefficient of DA object ina to see if any of them are INF (infinity).
   \param[in] ina Pointer to DA object to check
   \return True (non-zero) if any of the coefficients of ina is INF
*/
unsigned int daceIsInf(const DACEDA *ina)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
    {
        if(isinf(ia->cc))
            return true;
    }
    return false;
}
/** @}*/
