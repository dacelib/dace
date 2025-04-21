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
 *  daceeval.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/** \addtogroup DACE Core
 *  @{
 */

#include <math.h>

#include "dace/config.h"
#include "dace/dacebase.h"
#include "dace/daceaux.h"


/********************************************************************************
 *     DACE polynomial evaluation routines
 *********************************************************************************/

/*! Evaluate DA object ina by providing the value to use for each monomial in
    DA object inb.
    This is equivalent to a monomial-wise DA dot product.
   \param[in] ina Pointer to first DA object to evaluate
   \param[in] inb Pointer to second DA object to provide monomial values
   \return Monomial-wise dot product of both DA objects
   \sa daceMultiplyMonomials()
*/
double daceEvalMonomials(const DACEDA *ina, const DACEDA *inb)
{
	monomial *ipoa; unsigned int ilma, illa;
	monomial *ipob; unsigned int ilmb, illb;

	daceVariableInformation(ina, &ipoa, &ilma, &illa);
	daceVariableInformation(inb, &ipob, &ilmb, &illb);

	monomial *ib = ipob;
	monomial *const ibmax = ipob + ilmb;

	double res = 0.0;
	for(monomial *i = ipoa; i < ipoa + illa; i++)
	{
		while(ib->ii < i->ii && ib < ibmax)
			ib++;
		if(ib == ibmax) break;
		if(ib->ii == i->ii) res += ib->cc*i->cc;
	}

	return res;
}

/*! Perform partial evaluation of DA object ina by replacing independent variable
    number nvar by the value val.
   \param[in] ina Pointer to DA object to evaluate
   \param[in] nvar Number of the independent variable to replace (one-based)
   \param[in] val Value to replace independent variable with
   \param[out] inc Pointer to DA object to store the result of the partial evaluation
*/
void daceEvalVariable(const DACEDA *ina, const unsigned int nvar, const double val, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;

    if(nvar < 1 || nvar > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    const unsigned int ibase = DACECom.nomax+1;
    unsigned int j = nvar-1;
    if(nvar > DACECom.nv1)
        j -= DACECom.nv1;
    const unsigned int idiv = npown(ibase, j);

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double p[DACE_STATIC_NOMAX+1];
#else
    double *p = dacecalloc(DACECom.nomax+1, sizeof(double));
#endif
    p[0] = 1.0;
    for(unsigned int i = 1; i <= DACECom.nomax; i++)
        p[i] = p[i-1]*val;

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double cc[DACE_STATIC_NMMAX] = {0};
#else
    double *cc = dacecalloc(DACECom.nmmax, sizeof(double));
#endif

    if(nvar > DACECom.nv1)
    {
        for(monomial* i = ipoa; i < ipoa+illa; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            unsigned int ipow = (ic2/idiv)%ibase;
            j = DACECom.ia1[ic1]+DACECom.ia2[ic2-ipow*idiv];
            if(DACECom.ieo[j] <= DACECom_t.nocut)
                cc[j] += i->cc*p[ipow];
        }
    }
    else
    {
        for(monomial* i = ipoa; i < ipoa+illa; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            unsigned int ipow = (ic1/idiv)%ibase;
            j = DACECom.ia1[ic1-ipow*idiv]+DACECom.ia2[ic2];
            if(DACECom.ieo[j] <= DACECom_t.nocut)
                cc[j] += i->cc*p[ipow];
        }
    }

    dacePack(cc, inc);

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(cc);
    dacefree(p);
#endif
}

/*! Replace independent variable with index from by val times the independent
    variable with index to.
   \param[in] ina Pointer to DA object to evaluate
   \param[in] from Number of the independent variable to replace
   \param[in] to Number of the independent variable to change to
   \param[in] val Value to scale new independent variable with
   \param[out] inc Pointer to DA object to store the result of the replacement
*/
void daceReplaceVariable(const DACEDA *ina, const unsigned int from, const unsigned int to, const double val, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;

    if(from < 1 || from > DACECom.nvmax || to < 1 || to > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        daceCreateConstant(inc, 0.0);
        return;
    }
    else if(from == to)
    {
        daceScaleVariable(ina, from, val, inc);
        return;
    }

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    unsigned int p[DACE_STATIC_NOMAX+1];
    double pows[DACE_STATIC_NOMAX+1];
    double cc[DACE_STATIC_NMMAX] = {0};
#else
    unsigned int *p = dacecalloc(DACECom.nomax+1, sizeof(unsigned int));
    double *pows = dacecalloc(DACECom.nomax+1, sizeof(double));
    double *cc = dacecalloc(DACECom.nmmax, sizeof(double));
#endif

    pows[0] = 1.0;
    for(unsigned int i = 0; i < DACECom.nomax; i++)
        pows[i+1] = pows[i]*val;

    for(monomial* i = ipoa; i < ipoa+illa; i++)
    {
        daceDecode(i->ii, p);
        p[to] += p[from];
        const double c = pows[p[from]]*i->cc;
        p[from] = 0;
        cc[daceEncode(p)] += c;
    }

    dacePack(cc, inc);

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(cc);
    dacefree(p);
    dacefree(pows);
#endif
}

/*! Scale independent variable nvar by val.
   \param[in] ina Pointer to DA object to scale
   \param[in] nvar Number of the independent variable to scale
   \param[in] val Value to scale independent variable with
   \param[out] inc Pointer to DA object to store the result of the scaling
*/
void daceScaleVariable(const DACEDA *ina, const unsigned int nvar, const double val, DACEDA *inc)
{
    monomial *ipoc; unsigned int ilmc, illc;

    if(nvar < 1 || nvar > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceCopy(ina, inc);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);


#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double pows[DACE_STATIC_NOMAX+1];
#else
    double *pows = (double*) dacecalloc(DACECom.nomax+1, sizeof(double));
#endif
    pows[0] = 1.0;
    for(unsigned int i = 0; i < DACECom.nomax; i++) pows[i+1] = pows[i]*val;

    const unsigned int ibase = DACECom.nomax+1;
    unsigned int j = nvar-1;
    if(nvar > DACECom.nv1)
        j -= DACECom.nv1;
    const unsigned int idiv = npown(ibase, j);

    if(nvar > DACECom.nv1)
    {
        for(monomial* i = ipoc; i < ipoc+illc; i++)
        {
            const unsigned int ic2 = DACECom.ie2[i->ii];
            unsigned int ipow = (ic2/idiv)%ibase;
            i->cc *= pows[ipow];
        }
    }
    else
    {
        for(monomial* i = ipoc; i < ipoc+illc; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            unsigned int ipow = (ic1/idiv)%ibase;
            i->cc *= pows[ipow];
        }
    }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(pows);
#endif
}

/*! Translate independent variable nvar to (a*x + c).
   \param[in] ina Pointer to DA object to translate
   \param[in] nvar Number of the independent variable to translate
   \param[in] a Linear value to scale independent variable by
   \param[in] c Constant value to translate independent variable by
   \param[out] inc Pointer to DA object to store the result of the translation
*/
void daceTranslateVariable(const DACEDA *ina, const unsigned int nvar, const double a, const double c, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    if(nvar < 1 || nvar > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        daceCreateConstant(inc, 0.0);
        return;
    }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    unsigned int p[DACE_STATIC_NOMAX+1];
    double cc[DACE_STATIC_NMMAX] = {0};
    double powa[DACE_STATIC_NOMAX+1];
    double powc[DACE_STATIC_NOMAX+1];
    double binomial[(DACE_STATIC_NOMAX+1)*(DACE_STATIC_NOMAX+1)];
#else
    unsigned int *p = dacecalloc(DACECom.nomax+1, sizeof(unsigned int));
    double *cc = dacecalloc(DACECom.nmmax, sizeof(double));
    double *powa = dacecalloc(DACECom.nomax+1, sizeof(double));
    double *powc = dacecalloc(DACECom.nomax+1, sizeof(double));
    double *binomial = dacecalloc(((size_t)DACECom.nomax+1)*((size_t)DACECom.nomax+1), sizeof(double));
#endif

    // precompute powers of a and c
    powa[0] = 1.0;
    powc[0] = 1.0;
    for(unsigned int i = 0; i < DACECom.nomax; i++)
    {
        powa[i+1] = powa[i]*a;
        powc[i+1] = powc[i]*c;
    }

    // pre-compute binomial coefficients n choose k
    for(unsigned int n = 0; n <= DACECom.nomax; n++)
    {
        binomial[n*(DACECom.nomax+1)+0] = 1.0;
        binomial[n*(DACECom.nomax+1)+n] = 1.0;
        for(unsigned int k = 1; k < n; k++)
            binomial[n*(DACECom.nomax+1)+k] = binomial[(n-1)*(DACECom.nomax+1)+k-1] + binomial[(n-1)*(DACECom.nomax+1)+k];
    }

    for(monomial* i = ipoa; i < ipoa+illa; i++)
    {
        daceDecode(i->ii, p);
        unsigned int n = p[nvar-1];

        // shortcut the case when the monomial doesn't depend on nvar
        if(n == 0)
        {
            cc[i->ii] += i->cc;
            continue;
        }

        double *nk = binomial+n*(DACECom.nomax+1);
        double *pa = powa+n;
        double *pc = powc;
        for(unsigned int k = 0; k <= n; k++, nk++, pa--, pc++, p[nvar-1]--)
            cc[daceEncode(p)] += i->cc*(*nk)*(*pa)*(*pc);
    }

    dacePack(cc, inc);

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(cc);
    dacefree(p);
    dacefree(powa);
    dacefree(powc);
    dacefree(binomial);
#endif
}

/*! Compute an evaluation tree to efficiently evaluate several DA objects.
   \param[in] das C array of pointers to DA objects to evaluate
   \param[in] count Number of DA objects in das[]
   \param[out] ac C array of doubles containing compiled coefficients
   \param[out] nterm Pointer where to store the total number of terms in evaluation tree
   \param[out] nvar Pointer where to store the total number of variables in evaluation tree
   \param[out] nord Pointer where to store the maximum order in evaluation tree
   \sa DACE::compiledDA
   \note See the C++ interface for a detailed description of the resulting data
    arrays from this routine.
*/
void daceEvalTree(const DACEDA *das[], const unsigned int count, double ac[], unsigned int *nterm, unsigned int *nvar, unsigned int *nord)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    unsigned int nc[DACE_STATIC_NMMAX] = {0};
    unsigned int p[DACE_STATIC_NVMAX];
    unsigned int stack[DACE_STATIC_NOMAX];
#else
    unsigned int *nc = dacecalloc(DACECom.nmmax, sizeof(unsigned int));
    unsigned int *p = dacecalloc(DACECom.nvmax, sizeof(unsigned int));
    unsigned int *stack = dacecalloc(DACECom.nomax, sizeof(unsigned int));
#endif

    // mark all used monomials as new
    for(unsigned int i = 0; i < count; i++)
    {
        monomial* ipoa; unsigned int ilma, illa;
        daceVariableInformation(das[i], &ipoa, &ilma, &illa);
        for(monomial* j = ipoa; j < ipoa+illa; j++)
            nc[j->ii] = 2;
    }

    // make sure each term has a parent
    nc[0] = 1;  // constant part is the root, doesn't need a parent
    for(unsigned int i = 1; i < DACECom.nmmax; i++)
    {
        if(nc[i] != 2 ) continue;
        nc[i] = 1;
        daceDecode(i, p);
        // generate a ancestor tree for this entry
        int parent = 0;
        while(parent >= 0)
        {
            parent = -1;
            // find a parent
            for(unsigned int j = 0; j < DACECom.nvmax; j++)
            {
                if(p[j] == 0) continue;
                p[j]--;
                if(nc[daceEncode(p)] != 0)
                {
                    // parent already exists => done
                    parent = -1;
                    break;
                }
                p[j]++;
                parent = j;
            }
            // no parent found => create foster parent
            if(parent >= 0)
            {
                p[parent]--;
                nc[daceEncode(p)] = 1;
            }
        }
    }

    // constant terms are always stored
    nc[0] = 3;
    *nord = 0;
    *nvar = 0;
    *nterm = 1;
    ac[0] = ac[1] = 0.0;
    double *a = ac+2;
    for(unsigned int i = 0; i < count; i++)
        *(a++) = daceGetConstant(das[i]);     // evil C pointer magic: a++ evaluates to value of a, then increases it

    // higher order terms
    p[0] = 1;
    for(unsigned int i = 1; i < DACECom.nvmax; i++) p[i] = 0;
    int sp = 0;
    stack[sp] = 0;
    while(sp >= 0)
    {
        const unsigned int ic = daceEncode(p);
        if(nc[ic] == 1)
        {
            // store entry
            nc[ic] = 3;
            *nord = umax(*nord, sp+1);
            *nvar = umax(*nvar, stack[sp]+1);
            (*nterm)++;
            *(a++) = sp+1;          // +1 because old Fortran code returned 1 based indices
            *(a++) = stack[sp]+1;   // same here
            for(unsigned int i = 0; i < count; i++)
                *(a++) = daceGetCoefficient0(das[i], ic);

            // step forward if we can
            if(sp < (int)(DACECom.nomax-1))
            {
                sp++;
                stack[sp] = 0;
                p[0]++;
                continue;
            }
        }
        if(stack[sp] < DACECom.nvmax-1)
        {
            // step sideways
            p[stack[sp]]--;
            stack[sp]++;
            p[stack[sp]]++;
        }
        else
        {
            // step back
            p[stack[sp]]--;
            sp--;
        }
    }

#ifdef WITH_DEBUG
    // check if all monomials were properly encoded (should never fail except if above algorithm is implemented incorrectly)
    for(unsigned int i = 0; i < DACECom.nmmax; i++)
        if(nc[i] != 3 && nc[i] != 0)
        {
            daceSetError(__func__, DACE_PANIC, 8);
            exit(1);
        }
#endif

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(nc);
    dacefree(p);
    dacefree(stack);
#endif
}
/** @}*/
