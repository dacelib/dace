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
 *  dacemath.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/** \addtogroup DACE Core 
 *  @{
 */

// MS C library needs this to trigger it to define math constants
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#include "DA/dacebase.h"
#include "DA/daceaux.h"

// define various math constants in case they have not been defined by math.h
// these are non-standard C, but most C libraries have them
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef M_PI_2
#define M_PI_2 (1.57079632679489661923)
#endif


/********************************************************************************
 *     Basic DACE arithmetic operations
 *********************************************************************************/

/*! Perform addition of two DA objects.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] inb Pointer to the first DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina or inb.
 */
void daceAdd(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    if(!daceIsSameObject(ina, inc) && !daceIsSameObject(inb, inc))
    {
        daceWeightedSum(ina, 1.0, inb, 1.0, inc);
    }
    else
    {
        DACEDA idaadd;
        daceAllocateDA(&idaadd, 0);
        daceWeightedSum(ina, 1.0, inb, 1.0, &idaadd);
        daceCopy(&idaadd, inc);
        daceFreeDA(&idaadd);
    }
}

/*! Perform subtraction of two DA objects.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] inb Pointer to the first DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina or inb.
 */
void daceSubtract(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    if(!daceIsSameObject(ina, inc) && !daceIsSameObject(inb, inc))
    {
        daceWeightedSum(ina, 1.0, inb, -1.0, inc);
    }
    else
    {
        DACEDA idasub;
        daceAllocateDA(&idasub, 0);
        daceWeightedSum(ina, 1.0, inb, -1.0, &idasub);
        daceCopy(&idasub, inc);
        daceFreeDA(&idasub);
    }
}

/*! Perform multiplication of two DA objects.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] inb Pointer to the first DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina or inb.
 */
void daceMultiply(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
// These should use thread local storage (TLS) for multithread safe implementations
// see https://en.wikipedia.org/wiki/Thread-local_storage
#ifdef DACE_STATIC_MEMORY
    static DACE_THREAD_LOCAL double cc[DACE_STATIC_NMMAX] = {0};
    static DACE_THREAD_LOCAL extended_monomial emb[DACE_STATIC_NMMAX];
    static DACE_THREAD_LOCAL extended_monomial *ipbeg[DACE_STATIC_NOMAX+1];
    static DACE_THREAD_LOCAL extended_monomial *ipend[DACE_STATIC_NOMAX+1];
	static DACE_THREAD_LOCAL unsigned int nomax = 0;
	static DACE_THREAD_LOCAL unsigned int nvmax = 0;

    // make sure static memory is correctly allocated
    if(nomax != DACECom.nomax || nvmax != DACECom.nvmax)
    {
		nomax = DACECom.nomax;
		nvmax = DACECom.nvmax;
		ipbeg[0] = &emb[0];
		for(unsigned int i = 1; i <= DACECom.nomax; i++)
			ipbeg[i] = emb + daceCountMonomials(i - 1, DACECom.nvmax);
    }
#else
    static DACE_THREAD_LOCAL double *cc = NULL;
    static DACE_THREAD_LOCAL extended_monomial *emb = NULL;
    static DACE_THREAD_LOCAL extended_monomial **ipbeg = NULL;
    static DACE_THREAD_LOCAL extended_monomial **ipend = NULL;
    static DACE_THREAD_LOCAL unsigned int nomax = 0;
	static DACE_THREAD_LOCAL unsigned int nvmax = 0;

    // make sure static memory is correctly allocated
	if(nomax != DACECom.nomax || nvmax != DACECom.nvmax)
	{
		nomax = DACECom.nomax;
		nvmax = DACECom.nvmax;
		dacefree(cc);
        dacefree(emb);
        dacefree(ipbeg);
        dacefree(ipend);
        cc = (double*) dacecalloc(DACECom.nmmax, sizeof(double));
        emb = (extended_monomial*) dacecalloc(DACECom.nmmax, sizeof(extended_monomial));
        ipbeg = (extended_monomial**) dacecalloc(DACECom.nomax+1, sizeof(extended_monomial*));
        ipend = (extended_monomial**) dacecalloc(DACECom.nomax+1, sizeof(extended_monomial*));
		ipbeg[0] = &emb[0];
		for(unsigned int i = 1; i <= DACECom.nomax; i++)
			ipbeg[i] = emb + daceCountMonomials(i - 1, DACECom.nvmax);
	}
#endif

    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipob; unsigned int ilmb, illb;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inb, &ipob, &ilmb, &illb);

    // sort so that ina is the short DA vector
    if(illa>illb)
    {
        unsigned int t1;
        t1 = illb; illb = illa; illa = t1;
        t1 = ilmb; ilmb = ilma; ilma = t1;
        monomial* t2;
        t2 = ipoa; ipoa = ipob; ipob = t2;
    }

    for(unsigned int i = 0; i <= DACECom_t.nocut; i++) ipend[i] = ipbeg[i];

    // sort vector b by order
    for(monomial *ib = ipob; ib < ipob+illb; ib++)
    {
        const unsigned int noib = DACECom.ieo[ib->ii];
        if(noib > DACECom_t.nocut) continue;
        ipend[noib]->i1 = DACECom.ie1[ib->ii];
        ipend[noib]->i2 = DACECom.ie2[ib->ii];
        ipend[noib]->cc = ib->cc;
        ipend[noib]++;
    }

    // perform actual multiplication
    for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
    {
        const unsigned int i1ia = DACECom.ie1[ia->ii];
        const unsigned int i2ia = DACECom.ie2[ia->ii];
        const double ccia = ia->cc;
        for(int noib = DACECom_t.nocut-DACECom.ieo[ia->ii]; noib >= 0; noib--)
        {
            for(extended_monomial *ib = ipbeg[noib]; ib < ipend[noib]; ib++)
            {
                const unsigned int ic = DACECom.ia1[i1ia+ib->i1] + DACECom.ia2[i2ia+ib->i2];
                cc[ic] += ccia*ib->cc;
            }
        }
    }

    dacePack(cc, inc);
}

/*! Multiply two DA vectors component-wise, i.e. each monomial of ina with the corresponding monomial of inb
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] inb Pointer to the first DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina or inb.
   \sa daceEvalMonomials
 */
void daceMultiplyMonomials(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
	monomial *ipoa; unsigned int ilma, illa;
	monomial *ipob; unsigned int ilmb, illb;
	monomial *ipoc; unsigned int ilmc, illc;

	daceVariableInformation(ina, &ipoa, &ilma, &illa);
	daceVariableInformation(inb, &ipob, &ilmb, &illb);
	daceVariableInformation(inc, &ipoc, &ilmc, &illc);

	monomial *ib = ipob, *ic = ipoc;
	monomial *const ibmax = ipob + ilmb, *const icmax = ipoc + ilmc;

	for (monomial *i = ipoa; i < ipoa + illa; i++)
	{
		while (ib->ii < i->ii && ib < ibmax)
			ib++;
		if (ib == ibmax) break;
		if (ib->ii == i->ii)
		{
			if (ic >= icmax)
			{
				daceSetError(__func__, DACE_ERROR, 21);
				break;
			}
			ic->cc = i->cc*ib->cc;
			ic->ii = i->ii;
			ic++;
		}
	}
}

/*! Perform division of two DA objects.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] inb Pointer to the first DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina or inb.
 */
void daceDivide(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    DACEDA idadiv;

    daceAllocateDA(&idadiv, 0);
    daceMultiplicativeInverse(inb, &idadiv);
    daceMultiply(ina, &idadiv, inc);
    daceFreeDA(&idadiv);
}

/*! Square a DA object.
   \param[in] ina Pointer to the DA object to square
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceSquare(const DACEDA *ina, DACEDA *inb)
{
    daceMultiply(ina, ina, inb);
}

/*! Add constant to a DA object.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] ckon Constant value to add
   \param[out] inb Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inb can be the same as ina.
 */
void daceAddDouble(const DACEDA *ina, const double ckon, DACEDA *inb)
{
    if(!daceIsSameObject(ina, inb))
        daceCopy(ina, inb);

    daceSetCoefficient0(inb, 0, daceGetConstant(inb)+ckon);
}

/*! Subtract DA object from constant.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] ckon Constant value to subtract from
   \param[out] inb Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inb can be the same as ina.
 */
void daceDoubleSubtract(const DACEDA *ina, const double ckon, DACEDA *inb)
{
    daceMultiplyDouble(ina, -1.0, inb);
    daceSetCoefficient0(inb, 0, daceGetConstant(inb)+ckon);
}

/*! Subtract constant from a DA object.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] ckon Constant value to subtract
   \param[out] inb Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inb can be the same as ina.
 */
void daceSubtractDouble(const DACEDA *ina, const double ckon, DACEDA *inb)
{
    daceAddDouble(ina, -ckon, inb);
}

/*! Multiply constant and DA object.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] ckon Constant value to multiply by
   \param[out] inb Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inb can be the same as ina.
 */
void daceMultiplyDouble(const DACEDA *ina, const double ckon, DACEDA *inb)
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
            if(DACECom.ieo[ia->ii] > DACECom_t.nocut)
                continue;

            const double c = ia->cc*ckon;
            if(fabs(c) < DACECom_t.eps)
                continue;

            ib->cc = c;
            ib->ii = ia->ii;
            ib++;
        }
    }
    else
    {
        monomial *const ibmax = ipob+ilmb;
        for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
        {
            if(DACECom.ieo[ia->ii] > DACECom_t.nocut)
                continue;

            const double c = ia->cc*ckon;
            if(fabs(c) < DACECom_t.eps)
                continue;

            if(ib >= ibmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ib->cc = c;
            ib->ii = ia->ii;
            ib++;
        }
    }

    daceSetLength(inb, ib-ipob);
}

/*! Divide DA object by a constant.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] ckon Constant value to divide by
   \param[out] inb Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inb can be the same as ina.
 */
void daceDivideDouble(const DACEDA *ina, const double ckon, DACEDA *inb)
{
    if(ckon == 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 41);
        daceCreateConstant(inb, 0.0);
        return;
    }

    daceMultiplyDouble(ina, 1.0/ckon, inb);
}

/*! Divide constant by DA object.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] ckon Constant value to divide
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceDoubleDivide(const DACEDA *ina, const double ckon, DACEDA *inc)
{
    daceMultiplicativeInverse(ina, inc);
    daceMultiplyDouble(inc, ckon, inc);
}

/*! Divide a DA vector by a single variable to some power, if possible.
   \param[in] ina Pointer to the DA object to operate on
   \param[in] var Number of the independent variable by which to divide
   \param[in] p Power of independent variable
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceDivideByVariable(const DACEDA *ina, const unsigned int var, const unsigned int p, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    if(var < 1 || var > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);      
        daceCreateConstant(inc, 0.0);
        return;
    }

    // treat a few special cases
    if(p == 0)
    {
        // dividing by 1
        daceCopy(ina, inc);
        return;
    }
    else if(illa == 0)
    {
        // dividing 0 by anything
        daceCreateConstant(inc, 0.0);
        return;
    }
    else if(p > DACECom.nomax)
    {
        // dividing non-zero DA by too high a power
        daceSetError(__func__, DACE_ERROR, 42);      
        daceCreateConstant(inc, 0.0);
        return;
    }

    const unsigned int ibase = DACECom.nomax+1;
    unsigned int j = var-1;
    if(var > DACECom.nv1)
        j = j-DACECom.nv1;
    const unsigned int idiv = npown(ibase, j);

    monomial *ic = ipoc;
    monomial *const icmax = ipoc+ilmc;

    if(var > DACECom.nv1)
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic2/idiv)%ibase;
            if(ipow < p)
            {
                daceSetError(__func__, DACE_ERROR, 42);      
                daceCreateConstant(inc, 0.0);
                return;
            }
            if(ic >= icmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ic->ii = DACECom.ia1[ic1] + DACECom.ia2[ic2-p*idiv];
            ic->cc = i->cc;
            ic++;
        }
    }
    else
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic1/idiv)%ibase;
            if(ipow < p)
            {
                daceSetError(__func__, DACE_ERROR, 42);
                daceCreateConstant(inc, 0.0);
                return;
            }
            if(ic >= icmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ic->ii = DACECom.ia1[ic1-p*idiv] + DACECom.ia2[ic2];
            ic->cc = i->cc;
            ic++;
        }
    }

    daceSetLength(inc, ic-ipoc);
}

/*! Derivative of DA object with respect to a given independent variable.
   \param[in] idif Number of the independent variable with respect to which the
    derivative is taken
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceDifferentiate(const unsigned int idif, const DACEDA *ina, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    if(idif < 1 || idif > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        daceCreateConstant(inc, 0.0);
        return;
    }

    const unsigned int ibase = DACECom.nomax+1;
    unsigned int j = idif-1;
    if(idif > DACECom.nv1)
        j = j-DACECom.nv1;
    const unsigned int idiv = npown(ibase, j);

    monomial *ic = ipoc;
    monomial *const icmax = ipoc+ilmc;

    if(idif > DACECom.nv1)
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic2/idiv)%ibase;
            if(ipow == 0 || DACECom.ieo[i->ii] > DACECom_t.nocut+1)
                continue;
            if(ic >= icmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ic->ii = DACECom.ia1[ic1] + DACECom.ia2[ic2-idiv];
            ic->cc = i->cc*ipow;
            ic++;
        }
    }
    else
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic1/idiv)%ibase;
            if(ipow == 0 || DACECom.ieo[i->ii] > DACECom_t.nocut+1)
                continue;
            if(ic >= icmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ic->ii = DACECom.ia1[ic1-idiv] + DACECom.ia2[ic2];
            ic->cc = i->cc*ipow;
            ic++;
        }
    }

    daceSetLength(inc, ic-ipoc);
}

/*! Integral of DA object with respect to a given independent variable.
   \param[in] idif Number of the independent variable with respect to which the
    integral is taken
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceIntegrate(const unsigned int iint, const DACEDA *ina, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    if(iint < 1 || iint > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        daceCreateConstant(inc, 0.0);
        return;
    }

    const unsigned int ibase = DACECom.nomax+1;
    unsigned int j = iint-1;
    if(iint > DACECom.nv1)
        j = j-DACECom.nv1;
    const unsigned int idiv = npown(ibase, j);

    monomial *ic = ipoc;
    monomial *const icmax = ipoc+ilmc;

    if(iint > DACECom.nv1)
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            if(DACECom.ieo[i->ii] >= DACECom_t.nocut)
                continue;
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic2/idiv)%ibase;
            const double ccc = i->cc/(ipow+1);
            if(fabs(ccc) < DACECom_t.eps)
                continue;
            if(ic >= icmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ic->ii = DACECom.ia1[ic1] + DACECom.ia2[ic2+idiv];
            ic->cc = ccc;
            ic = ic+1;
        }
    }
    else
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            if(DACECom.ieo[i->ii] >= DACECom_t.nocut)
                continue;
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic1/idiv)%ibase;
            const double ccc = i->cc/(ipow+1);
            if(fabs(ccc) < DACECom_t.eps)
                continue;
            if(ic >= icmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ic->ii = DACECom.ia1[ic1+idiv] + DACECom.ia2[ic2];
            ic->cc = ccc;
            ic = ic+1;
        }
    }

    daceSetLength(inc, ic-ipoc);
}


/********************************************************************************
 *     DACE intrinsic function routines
 *********************************************************************************/

/*! Truncate the constant part of a DA object to an integer.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceTruncate(const DACEDA *ina, DACEDA *inc)
{
    daceCopy(ina, inc);
    daceSetCoefficient0(inc, 0, rint(daceGetConstant(inc)));
}

/*! Round the constant part of a DA object to an integer.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceRound(const DACEDA *ina, DACEDA *inc)
{
    daceCopy(ina, inc);
    daceSetCoefficient0(inc, 0, round(daceGetConstant(inc)));
}

/*! Modulo the constant part of a DA object by p.
   \param[in] ina Pointer to the DA object to operate on
   \param[in] p Value with respect to which to compute the modulo
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceModulo(const DACEDA *ina, const double p, DACEDA *inc)
{
    daceCopy(ina, inc);
    daceSetCoefficient0(inc, 0, fmod(daceGetConstant(inc),p));
}

/*! Raise a DA object to the p-th power.
   \param[in] ina Pointer to the DA object to operate on
   \param[in] p Power to which to raise the DA object
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void dacePowerDouble(const DACEDA *ina, const double p, DACEDA *inc)
{
    // check simple cases
    if(p == 0.0)
    {
        daceCreateConstant(inc, 1.0);
        return;
    }
    else if(p == (int)p)
    {
        dacePower(ina, (int)p, inc);
        return;
    }

    const double a0 = daceGetConstant(ina);
    if(a0 <= 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 43); 
        daceCreateConstant(inc, 0.0);
        return;
    }

#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double *xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = pow(a0, p);
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
        xf[i] = xf[i-1]/i*(p-(i-1));

    daceDivideDouble(ina, a0, inc);     // more accurate than including a0 in series (uses non-linear part in EvaluateSeries)
    daceEvaluateSeries(inc, xf, inc);
#ifndef DACE_STATIC_MEMORY
   dacefree(xf);
#endif
}

/*! Raise a DA object to the p-th integer power.
   \param[in] ina Pointer to the DA object to operate on
   \param[in] p Power to which to raise the DA object
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void dacePower(const DACEDA *ina, const int np, DACEDA *inc)
{
    DACEDA itemp;

    // handle some common simple cases directly
    switch(np)
    {
        case 0:
            daceCreateConstant(inc, 1.0);
            return;

        case 1:
            daceCopy(ina, inc);
            return;

        case -1:
            daceMultiplicativeInverse(ina, inc);
            return;
    }

    // handle all other cases, again with common special cases hard coded
    switch(abs(np))
    {
        case 2:
            daceSquare(ina, inc);
            break;

        case 3:
            daceAllocateDA(&itemp, 0);
            daceSquare(ina, &itemp);
            daceMultiply(ina, &itemp, inc);
            daceFreeDA(&itemp);
            break;

        case 4:
            daceAllocateDA(&itemp, 0);
            daceSquare(ina, &itemp);
            daceSquare(&itemp, inc);
            daceFreeDA(&itemp);
            break;

        default:
            daceAllocateDA(&itemp, 0);
            daceCopy(ina, &itemp);
            daceCreateConstant(inc, 1.0);
            unsigned int inp = abs(np);
            while(inp)
            {
                if(inp & 1u)
                    daceMultiply(inc, &itemp, inc);
                inp >>= 1;
                if(inp)
                    daceSquare(&itemp, &itemp);
            }
            daceFreeDA(&itemp);
    }

    if(np < 0)
        daceMultiplicativeInverse(inc, inc);
}

/*! Take the np-th root of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[in] np Root to take of the DA object
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceRoot(const DACEDA *ina, const int np, DACEDA *inc)
{
    if(np == 0)
    {
        daceSetError(__func__, DACE_ERROR, 44);
        daceCreateConstant(inc, 0.0);
        return;
    }

    const double a0 = daceGetConstant(ina);
    const unsigned int iodd = abs(np) & 1u;

    if((iodd == 0) && (a0 <= 0.0))
    {
        daceSetError(__func__, DACE_ERROR, 45);
        daceCreateConstant(inc, 0.0);
        return;
    }
    else if((iodd == 1) && (a0 == 0.0))
    {
        daceSetError(__func__, DACE_ERROR, 46);
        daceCreateConstant(inc, 0.0);
        return;
    }

    double cr = 1.0/np;
#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double *xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = copysign(pow(fabs(a0), cr), a0);
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = xf[i-1]/i*cr;
        cr--;
    }

    daceDivideDouble(ina, a0, inc);     // more accurate than including a0 in series (uses non-linear part in EvaluateSeries)
    daceEvaluateSeries(inc, xf, inc);
#ifndef DACE_STATIC_MEMORY
   dacefree(xf);
#endif
}

/*! Compute the multiplicative inverse of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceMultiplicativeInverse(const DACEDA *ina, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    
    if(a0 == 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 41);
        daceCreateConstant(inc, 0.0);
        return;
    }
    
    if(DACECom_t.nocut < 5)
    {
        // lower orders: compute series directly
        daceMultiplicativeInverse0(ina, inc, a0);
    }
    else
    {
        // higher orders: use iteration
        const unsigned int nocut = DACECom_t.nocut;
        DACECom_t.nocut = 2;
        daceMultiplicativeInverse0(ina, inc, a0);
        DACEDA temp;
        daceAllocateDA(&temp, 0);
        for(unsigned int ord = 3; ord <= nocut; ord *= 2)
        {
            DACECom_t.nocut = umin(nocut, 2*ord-1);
            daceMultiply(ina, inc, &temp);
            daceDoubleSubtract(&temp, 2.0, &temp);
            daceMultiply(inc, &temp, inc);
        }
        daceFreeDA(&temp);
    }
}

/*! Compute the multiplicative inverse of a DA object using series expansion.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \param[in] a0 Constant part of ina
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceMultiplicativeInverse0(const DACEDA *ina, DACEDA *inc, const double a0)
{
    daceDivideDouble(ina, a0, inc);

#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double *xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = 1.0/a0;
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
        xf[i] = -xf[i-1];

    daceEvaluateSeries(inc, xf, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(xf);
#endif
}

/*! Compute the square root of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceSquareRoot(const DACEDA *ina, DACEDA *inc)
{
    daceRoot(ina, 2, inc);
}

/*! Compute the inverse square root of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceInverseSquareRoot(const DACEDA *ina, DACEDA *inc)
{
    daceRoot(ina, -2, inc);
}

/*! Compute the cubic root of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceCubicRoot(const DACEDA *ina, DACEDA *inc)
{
    daceRoot(ina, 3, inc);
}

/*! Compute the inverse cubic root of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceInverseCubicRoot(const DACEDA *ina, DACEDA *inc)
{
    daceRoot(ina, -3, inc);
}

/*! Compute the hypothenuse of two DA objects.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] inb Pointer to the second DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina or inb.
 */
void daceHypotenuse(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    DACEDA itemp1, itemp2;

    daceAllocateDA(&itemp1, 0);
    daceAllocateDA(&itemp2, 0);
    daceSquare(ina, &itemp1);
    daceSquare(inb, &itemp2);
    daceAdd(&itemp1, &itemp2, inc);
    daceRoot(inc, 2, inc);
    daceFreeDA(&itemp2);
    daceFreeDA(&itemp1);
}

/*! Compute the exponential of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceExponential(const DACEDA *ina, DACEDA *inc)
{
#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = exp(daceGetConstant(ina));
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
        xf[i] = xf[i-1]/i;

    daceEvaluateSeries(ina, xf, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(xf);
#endif
}

/*! Compute the natural logarithm root of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceLogarithm(const DACEDA *ina, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    if(a0 <= 0)
    {
        daceSetError(__func__, DACE_ERROR, 47);
        daceCreateConstant(inc, 0.0);
        return;
    }

#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    daceDivideDouble(ina, a0, inc);
    xf[0] = log(a0);
    xf[1] = 1.0;
    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = -xf[i-1]/i*(i-1);
    }

    daceEvaluateSeries(inc, xf, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(xf);
#endif
}

/*! Compute the logarithm with respect to base b of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[in] b Base of the logarithm to use
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceLogarithmBase(const DACEDA *ina, const double b, DACEDA *inc)
{
    if(b <= 0)
    {
        daceSetError(__func__, DACE_ERROR, 48);
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceLogarithm(ina, inc);
    daceMultiplyDouble(inc, 1.0/log(b), inc);
}

/*! Compute the decadic logarithm of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceLogarithm10(const DACEDA *ina, DACEDA *inc)
{
    daceLogarithmBase(ina, 10.0, inc);
}

/*! Compute the binary logarithm of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceLogarithm2(const DACEDA *ina, DACEDA *inc)
{
    daceLogarithmBase(ina, 2.0, inc);
}

/*! Compute the sine of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceSine(const DACEDA *ina, DACEDA *inc)
{
#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif
    const double a0 = daceGetConstant(ina);

    xf[0] = sin(a0);
    xf[1] = cos(a0);
    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = -xf[i-2]/(i*(i-1));
    }

    daceEvaluateSeries(ina, xf, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(xf);
#endif
}

/*! Compute the cosine of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceCosine(const DACEDA *ina, DACEDA *inc)
{
#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif
    const double a0 = daceGetConstant(ina);

    xf[0] = cos(a0);
    xf[1] = -sin(a0);
    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = -xf[i-2]/(i*(i-1));
    }

    daceEvaluateSeries(ina, xf, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(xf);
#endif
}

/*! Compute the tangent of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceTangent(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;

    if(cos(daceGetConstant(ina)) == 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 49);
        daceCreateConstant(inc, 0.0);
        return;
    }
    
    daceAllocateDA(&itemp, 0);
    daceSine(ina, &itemp);
    daceCosine(ina, inc);
    daceDivide(&itemp, inc, inc);
    daceFreeDA(&itemp);
}

/*! Compute the arcsine of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceArcSine(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;

    if(fabs(daceGetConstant(ina)) >= 1.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }
    
    daceAllocateDA(&itemp, 0);
    daceSquare(ina, &itemp);
    daceDoubleSubtract(&itemp, 1.0, &itemp);
    daceSquareRoot(&itemp, &itemp);
    daceDivide(ina, &itemp, inc);
    daceArcTangent(inc, inc);
    daceFreeDA(&itemp);
}

/*! Compute the arccosine of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceArcCosine(const DACEDA *ina, DACEDA *inc)
{
    if(fabs(daceGetConstant(ina)) >= 1.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }
    daceArcSine(ina, inc);
    daceDoubleSubtract(inc, M_PI_2, inc);
}

/*! Compute the arctangent of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceArcTangent(const DACEDA *ina, DACEDA *inc)
{
    DACEDA iarg;
#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1] = {0};
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif
    const double a0 = daceGetConstant(ina);

    daceAllocateDA(&iarg, 0);
    daceMultiplyDouble(ina, a0, &iarg);
    daceAddDouble(&iarg, 1.0, &iarg);
    daceSubtractDouble(ina, a0, inc);
    daceDivide(inc, &iarg, &iarg);
    
    double s = 1.0;
    xf[0] = atan(a0);
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i+=2)
    {
        xf[i] = s/i;
        s = -s;
    }

    daceEvaluateSeries(&iarg, xf, inc);
    daceFreeDA(&iarg);
#ifndef DACE_STATIC_MEMORY
    dacefree(xf);
#endif
}

/*! Arctangent of ina/inb with proper sign in [-pi, pi]. This function follows
    the C standard atan2(y,x) function syntax.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] ina Pointer to the second DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceArcTangent2(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    const double cx = daceGetConstant(inb);
    const double cy = daceGetConstant(ina);

    if(cx == 0.0 && cy == 0.0)
    {
        daceCreateConstant(inc, 0.0);
    }
    else
    {
        if(fabs(cy) > fabs(cx))
        {
            daceDivide(inb, ina, inc);
            daceArcTangent(inc, inc);
            if(cy < 0.0)
            {
                daceDoubleSubtract(inc, -M_PI_2, inc);
            }
            else
            {
                daceDoubleSubtract(inc, M_PI_2, inc);
            }
        }
        else
        {
            daceDivide(ina, inb, inc);
            daceArcTangent(inc, inc);
            if(cx < 0.0)
            {
                if(cy > 0.0)
                {
                    daceAddDouble(inc, M_PI, inc);
                }
                else
                {
                    daceAddDouble(inc, -M_PI, inc);
                }
            }
        }
    }
}

/*! Compute the hyperbolic sine of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceHyperbolicSine(const DACEDA *ina, DACEDA *inc)
{
#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    const double a0 = daceGetConstant(ina);
    xf[0] = sinh(a0);
    xf[1] = cosh(a0);
    
    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = xf[i-2]/(i*(i-1));
    }

    daceEvaluateSeries(ina, xf, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(xf);
#endif
}

/*! Compute the hyperbolic cosine of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceHyperbolicCosine(const DACEDA *ina, DACEDA *inc)
{
#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    const double a0 = daceGetConstant(ina);
    xf[0] = cosh(a0);
    xf[1] = sinh(a0);
    
    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = xf[i-2]/(i*(i-1));
    }

    daceEvaluateSeries(ina, xf, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(xf);
#endif
}

/*! Compute the hyperbolic tangent of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceHyperbolicTangent(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;
    
    daceAllocateDA(&itemp, 0);
    daceHyperbolicSine(ina, &itemp);
    daceHyperbolicCosine(ina, inc);
    daceDivide(&itemp, inc, inc);
    daceFreeDA(&itemp);
}

/*! Compute the hyperbolic arcsince of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceHyperbolicArcSine(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;
    
    daceAllocateDA(&itemp, 0);
    daceSquare(ina, inc);
    daceAddDouble(inc, 1.0, &itemp);
    daceSquareRoot(&itemp, inc);
    daceAdd(ina, inc, &itemp);
    daceLogarithm(&itemp, inc);
    daceFreeDA(&itemp);
}

/*! Compute the hyperbolic arccosine of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceHyperbolicArcCosine(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;
    
    if(daceGetConstant(ina) <= 1.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }
    
    daceAllocateDA(&itemp, 0);
    daceSquare(ina, inc);
    daceSubtractDouble(inc, 1.0, &itemp);
    daceSquareRoot(&itemp, inc);
    daceAdd(ina, inc, &itemp);
    daceLogarithm(&itemp, inc);
    daceFreeDA(&itemp);
}

/*! Compute the hyperbolic arctangent of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceHyperbolicArcTangent(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;
    
    if(fabs(daceGetConstant(ina)) >= 1.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }
    
    daceAllocateDA(&itemp, 0);
    daceAddDouble(ina, 1.0, &itemp);
    daceDoubleSubtract(ina, 1.0, inc);
    daceDivide(&itemp, inc, inc);
    daceLogarithm(inc, &itemp);
    daceMultiplyDouble(&itemp, 0.5, inc);
    daceFreeDA(&itemp);
}

/*! Compute the error function of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceErrorFunction(const DACEDA *ina, DACEDA *inc)
{
#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    const double a0 = daceGetConstant(ina);
    double factor = 2.0*exp(-a0*a0)/sqrt(M_PI);
    xf[0] = erf(a0);
    xf[1] = factor;
    double Hi2 = 1.0;       // Hermite polynomial H_{i-2} = H_0
    double Hi1 = 2.0*a0;    // Hermite polynomial H_{i-1} = H_1

    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        factor /= -((double)i);
        xf[i] = factor*Hi1;
        const double temp = 2.0*a0*Hi1 - 2.0*(i-1)*Hi2;     // recursion relation: H_i = 2*x*H_{i-1} - 2*(i-1)*H_{i-2}
        Hi2 = Hi1;
        Hi1 = temp;
    }

    daceEvaluateSeries(ina, xf, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(xf);
#endif
}

/*! Compute the complementary error function of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceComplementaryErrorFunction(const DACEDA *ina, DACEDA *inc)
{
#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    const double a0 = daceGetConstant(ina);
    double factor = -2.0*exp(-a0*a0)/sqrt(M_PI);
    xf[0] = erfc(a0);
    xf[1] = factor;
    double Hi2 = 1.0;       // Hermite polynomial H_{i-2} = H_0
    double Hi1 = 2.0*a0;    // Hermite polynomial H_{i-1} = H_1

    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
		factor /= -((double)i);
		xf[i] = factor*Hi1;
        const double temp = 2.0*a0*Hi1 - 2.0*(i-1)*Hi2;     // recursion relation: H_i = 2*x*H_{i-1} - 2*(i-1)*H_{i-2}
        Hi2 = Hi1;
        Hi1 = temp;
    }

    daceEvaluateSeries(ina, xf, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(xf);
#endif
}

/*! Compute the Bessel function J_n of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceBesselJFunction(const DACEDA *ina, const int n, DACEDA *inc)
{
#ifdef DACE_STATIC_MEMORY
    double bz[2*DACE_STATIC_NOMAX+1];
#else
    double* bz = (double*) dacecalloc(2*DACECom_t.nocut+1, sizeof(double));
#endif

    const double a0 = daceGetConstant(ina);
    for(unsigned int i = 0; i < 2*DACECom_t.nocut+1; i++)
        bz[i] = jn(n-DACECom_t.nocut+i, a0);		// XXX: Mauro: MSVC complains that "warning C4996: 'jn': The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _jn. See online help for details." But at least on my Mac these don't seem to exist. I recommend we leave it as is (remove comment after fix).

    daceEvaluateBesselFunction(ina, bz, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(bz);
#endif
}

/*! Compute the Bessel function Y_n of a DA object.
   \param[in] ina Pointer to the DA object to operate on
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceBesselYFunction(const DACEDA *ina, const int n, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    if(a0 <= 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

#ifdef DACE_STATIC_MEMORY
    double bz[2*DACE_STATIC_NOMAX+1];
#else
    double* bz = (double*) dacecalloc(2*DACECom_t.nocut+1, sizeof(double));
#endif

    for(unsigned int i = 0; i < 2*DACECom_t.nocut+1; i++)
        bz[i] = yn(n-DACECom_t.nocut+i, a0);		// XXX: Mauro, same here as above. MSVC recommends using _yn

    daceEvaluateBesselFunction(ina, bz, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(bz);
#endif
}

/*! Evaluate a Bessel function with coefficients bz with the non-constant part of ina.
   \param[in] ina Pointer to the DA object to operate on
   \param[in] bz C array of nocut+1 elements containing the Bessel coefficients of the polynomial
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceEvaluateBesselFunction(const DACEDA *ina, const double bz[], DACEDA *inc)
{
#ifdef DACE_STATIC_MEMORY
    double xf[DACE_STATIC_NOMAX+1];
    double binomial[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
    double* binomial = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = bz[DACECom_t.nocut];
    binomial[0] = 1.0;
    double factor = 1.0;
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
    {
        factor *= 0.5/i;
        // calculate binomial coefficients i choose j
        binomial[i] = 1.0;
        for(unsigned int j = i-1; j > 0; j--)
            binomial[j] += binomial[j-1];
        // Calculate n-th derivative of Bessel function C, see http://dlmf.nist.gov/10.6
        // bz contains values of C_{n-o} to C_{n+o} of constant part of ina
        double sign = 1.0;
        xf[i] = 0.0;
        for(unsigned int j = 0; j <= i; j++)
        {
            xf[i] += binomial[j]*sign*bz[DACECom_t.nocut-i+2*j];
            sign *= -1.0;
        }
        xf[i] *= factor;
    }

    daceEvaluateSeries(ina, xf, inc);
#ifndef DACE_STATIC_MEMORY
    dacefree(binomial);
    dacefree(xf);
#endif
}

/*! Evaluate a polynomial with coefficients xf with the non-constant part of ina.
   \param[in] ina Pointer to the DA object to operate on
   \param[in] xf C array of nocut+1 elements containing the coefficients of the polynomial
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is aliasing safe, i.e. inc can be the same as ina.
 */
void daceEvaluateSeries(const DACEDA *ina, const double xf[], DACEDA *inc)
{
    DACEDA inon;
    const unsigned int nocut = DACECom_t.nocut;

    daceAllocateDA(&inon, 0);
    daceCopy(ina, &inon);
    daceSetCoefficient0(&inon, 0, 0.0);

    DACECom_t.nocut = 1;
    daceMultiplyDouble(&inon, xf[nocut], inc);
    daceAddDouble(inc, xf[nocut-1], inc);

    // evaluate series
    for(int i = nocut-2; i >= 0; i--)
    {
        DACECom_t.nocut = nocut-i;
        daceMultiply(&inon, inc, inc);
        daceAddDouble(inc, xf[i], inc);
    }

    DACECom_t.nocut = nocut;
    daceFreeDA(&inon);
}

/*! Compute the weighted sum of two DA objects.
   \param[in] ina Pointer to the first DA object to operate on
   \param[in] afac Weighting factor to multiply ina by
   \param[in] inb Pointer to the second DA object to operate on
   \param[in] bfac Weighting factor to multiply inb by
   \param[out] inc Pointer to the DA object to store the result in
   \note This routine is NOT aliasing safe! So inc MUST BE DIFFERENT from ina
    and inb.
 */
void daceWeightedSum(const DACEDA *ina, const double afac, const DACEDA *inb, const double bfac, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipob; unsigned int ilmb, illb;
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inb, &ipob, &ilmb, &illb);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    monomial *ia = ipoa, *ib = ipob, *ic = ipoc;
    monomial *const iamax = ipoa+illa, *const ibmax = ipob+illb, *const icmax = ipoc+ilmc;

    if(illa > 0 && illb > 0)
    {
        // both polynomials have coefficients, merge until one runs out
        unsigned int ja = ia->ii;
        unsigned int jb = ib->ii;
        while(true)
        {
            if(ja == jb)
            {
                // add the two terms
                if(DACECom.ieo[ja] <= DACECom_t.nocut)
                {
                    const double ccc = ia->cc*afac + ib->cc*bfac;
                    if(fabs(ccc) >= DACECom_t.eps)
                    {
                        if(ic >= icmax)
                        {
                            daceSetError(__func__, DACE_ERROR, 21);
                            daceSetLength(inc, ilmc);
                            return;
                        }
                        ic->cc = ccc;
                        ic->ii = ia->ii;
                        ic++;
                    }
                }
                ia++; ib++;
                if(ia >= iamax || ib >= ibmax) break;
                ja = ia->ii;
                jb = ib->ii;
            }
            else if(ja < jb)
            {
                // store term a
                if(DACECom.ieo[ja] <= DACECom_t.nocut)
                {
                    const double ccc = ia->cc*afac;
                    if(fabs(ccc) >= DACECom_t.eps)
                    {
                        if(ic >= icmax)
                        {
                            daceSetError(__func__, DACE_ERROR, 21);
                            daceSetLength(inc, ilmc);
                            return;
                        }
                        ic->cc = ccc;
                        ic->ii = ia->ii;
                        ic++;
                    }
                }
                ia++;
                if(ia >= iamax) break;
                ja = ia->ii;
            }
            else
            {
                // store term b
                if(DACECom.ieo[jb] <= DACECom_t.nocut)
                {
                    const double ccc = ib->cc*bfac;
                    if(fabs(ccc) >= DACECom_t.eps)
                    {
                        if(ic >= icmax)
                        {
                            daceSetError(__func__, DACE_ERROR, 21);
                            daceSetLength(inc, ilmc);
                            return;
                        }
                        ic->cc = ccc;
                        ic->ii = ib->ii;
                        ic++;
                    }
                }
                ib++;
                if(ib >= ibmax) break;
                jb = ib->ii;
            }
        }
    }

    // copy any remaining terms from either ina or inb
    monomial *ismin, *ismax;
    double fac;
    if(ia < iamax)
    {
        ismin = ia;
        ismax = iamax;
        fac = afac;
    }
    else
    {
        ismin = ib;
        ismax = ibmax;
        fac = bfac;
    }

    for(monomial *is = ismin; is < ismax; is++)
    {
        if(DACECom.ieo[is->ii] <= DACECom_t.nocut)
        {
            const double ccc = is->cc*fac;
            if(fabs(ccc) >= DACECom_t.eps)
            {
                if(ic >= icmax)
                {
                    daceSetError(__func__, DACE_ERROR, 21);
                    daceSetLength(inc, ilmc);
                    return;
                }
                ic->cc = ccc;
                ic->ii = is->ii;
                ic++;
            }
        }
    }

    daceSetLength(inc, ic-ipoc);
}
/** @}*/
