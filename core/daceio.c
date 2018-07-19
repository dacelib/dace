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
 *  daceio.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/** \addtogroup DACE Core 
 *  @{
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dace/config.h"
#include "dace/dacebase.h"
#include "dace/daceaux.h"

/********************************************************************************
 *     DACE input/output routines
 *********************************************************************************/

/*! Print the DA object ina to string strs (of line length DACE_STRLEN).
   \param[in] ina Pointer to the DA object to be printed
   \param[out] strs C array of size (nmmax+2)*DACE_STRLEN containing the
    zero-terminated lines of length DACE_STRLEN
   \param[out] nstrs Pointer where to store the final number of strings printed
   \sa daceRead
   \sa dacePrint
   \note The format of the output is written in DACE format which is loosely
    based on the format used by COSY INFINITY but is not fully compatible to it.
 */
void daceWrite(const DACEDA *ina, char *strs, unsigned int *nstrs)
{
    const char* BEGSTR = "     I  COEFFICIENT              ORDER EXPONENTS";
    const char* ENDSTR = "------------------------------------------------";
    const char* ZEROSTR = "        ALL COEFFICIENTS ZERO";
    monomial *ipoa; unsigned int ilma, illa;
    char* s = strs;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    if(illa > DACECom.nmmax)
    {
       daceSetError(__func__, DACE_ERROR, 33);
       illa = DACECom.nmmax;
    }

    *nstrs = 0;

    if(illa == 0)
    {
        snprintf(s, DACE_STRLEN, "%s", ZEROSTR);
        (*nstrs)++;
        s += DACE_STRLEN;
    }
    else
    {
        snprintf(s, DACE_STRLEN, "%s", BEGSTR);
        (*nstrs)++;
        s += DACE_STRLEN;
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
        unsigned int j[DACE_STATIC_NVMAX];
#else
        unsigned int *j = dacecalloc(DACECom.nvmax, sizeof(unsigned int));
#endif
        unsigned int iout = 1;
        for(unsigned int ioa = 0; ioa <= DACECom.nomax; ioa++ )
        {
            for(monomial* m = ipoa; m < ipoa+illa; m++)
            {
                if(DACECom.ieo[m->ii] != ioa) continue;
                daceDecode(m->ii, j);
                snprintf(s, DACE_STRLEN, "%6u  %24.16e%4u ", iout, m->cc, ioa);
                for(unsigned int i = 0; i < DACECom.nvmax; i++)
                    snprintf(s+37+3*i, DACE_STRLEN-37-3*i, " %2u", j[i]);
                s += DACE_STRLEN;
                (*nstrs)++;
                iout++;
            }
        }
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
        dacefree(j);
#endif
    }

    snprintf(s, DACE_STRLEN, "%s", ENDSTR);
    (*nstrs)++;
}

/*! Read a DA object ina from a human readable string representation in strs,
    containing nstrs contiguous zero-terminated lines of line length DACE_STRLEN.
   \param[out] ina Pointer to the DA object to read into
   \param[in] strs C array of size nstrs*DACE_STRLEN containing the
    zero-terminated lines of length DACE_STRLEN
   \param[in] nstrs Number of lines in strs
   \sa daceWrite
   \sa dacePrint
   \note This routine can read both DACE output strings as well as some COSY INFINITY
    formated strings. COSY INFINITY input is limited to less than 16 variables
    (i.e. a single line per coefficient).
 */
void daceRead(DACEDA *ina, char* strs, unsigned int nstrs)
{
    const char* BEGSTR = "     I  COEFFICIENT              ORDER EXPONENTS";
    const char* ENDSTR = "------------------------------------------------";
    const char* ZEROSTR = "        ALL COEFFICIENTS ZERO";
    const char* COSYBEGSTR = "     I  COEFFICIENT            ORDER EXPONENTS";
    const char* COSYZEROSTR = "        ALL COMPONENTS ZERO";
    const unsigned int DACE_COEFF_LEN = 24;
    const unsigned int COSY_COEFF_LEN = 22;
    char *str = strs;

    daceCreateConstant(ina, 0.0);

    if(nstrs < 1)
    {
        daceSetError(__func__, DACE_ERROR, 34);
        return;
    }
    else if(strncmp(str, ZEROSTR, strlen(ZEROSTR)) == 0 ||
            strncmp(str, COSYZEROSTR, strlen(COSYZEROSTR)) == 0)
    {
        // looks like an empty da, nothing to do
        return;
    }
    else
    {
        bool cosy;
        unsigned int coefflen;

        if(strncmp(str, BEGSTR, strlen(BEGSTR)) == 0)
        {
            // looks like a DACE DA vector
            coefflen = DACE_COEFF_LEN;
            cosy = false;
        }
        else if(strncmp(str, COSYBEGSTR, strlen(COSYBEGSTR)) == 0)
        {
            // looks like a COSY DA vector
            coefflen = COSY_COEFF_LEN;
            cosy = true;
        }
        else
        {
            // no idea what this input is supposed to be
            daceSetError(__func__, DACE_ERROR, 32);
            return;
        }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
        double cc[DACE_STATIC_NMMAX] = {0};
        unsigned int j[DACE_STATIC_NVMAX];
#else
        double *cc = dacecalloc(DACECom.nmmax, sizeof(double));
        unsigned int *j = dacecalloc(DACECom.nvmax, sizeof(unsigned int));
#endif
        str += DACE_STRLEN;
        for(unsigned int iin = 1; iin < nstrs; iin++, str += DACE_STRLEN)
        {
            char *s = str;
            unsigned int len = strlen(s);
            char temp;

            // check if it's the end of data line by comparing characters 5-35 (works for both COSY and DACE format)
            if(len < 4)
            {
                daceSetError(__func__, DACE_ERROR, 32);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
                dacefree(cc);
                dacefree(j);
#endif
                return;
            }
            if(strncmp(s+4, ENDSTR, 31) == 0) break;

            // read fields one by one

            // line number
            if(len < 6)
            {
                daceSetError(__func__, DACE_ERROR, 32);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
                dacefree(cc);
                dacefree(j);
#endif
                return;
            }
            temp = s[6];
            s[6] = '\0';
            const unsigned int ii = (unsigned int)strtoul(s, NULL, 10);
            s[6] = temp;
            s += 6;
            len -= 6;

            // skip 2 spaces
            if(len < 2)
            {
                daceSetError(__func__, DACE_ERROR, 32);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
                dacefree(cc);
                dacefree(j);
#endif
                return;
            }
            s += 2;
            len -= 2;

            // coefficient
            if(len < coefflen)
            {
                daceSetError(__func__, DACE_ERROR, 32);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
                dacefree(cc);
                dacefree(j);
#endif
                return;
            }
            temp = s[coefflen];
            s[coefflen] = '\0';
            const double c = strtod(s, NULL);
            s[coefflen] = temp;
            s += coefflen;
            len -= coefflen;

            // order
            if(len < 4)
            {
                daceSetError(__func__, DACE_ERROR, 32);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
                dacefree(cc);
                dacefree(j);
#endif
                return;
            }
            temp = s[4];
            s[4] = '\0';
            const unsigned int io1 = (unsigned int)strtoul(s, NULL, 10);
            s[4] = temp;
            s += 4;
            len -= 4;

            // skip 1 space
            if(len < 1)
            {
                daceSetError(__func__, DACE_ERROR, 32);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
                dacefree(cc);
                dacefree(j);
#endif
                return;
            }
            s += 1;
            len -= 1;

            // read exponents one by one
            for(unsigned int i = 0; i < DACECom.nvmax; i++)
            {
                j[i] = 0;   // in case the read below fails

                // single space (in COSY case only in front of every other variable)
                if(len >= 1 && (!cosy || (i&1) == 0))
                {
                    s += 1;
                    len -= 1;
                }

                // single exponent
                if(len >= 2)
                {
                    temp = s[2];
                    s[2] = '\0';
                    j[i] = (unsigned int) strtoul(s, NULL, 10);
                    s[2] = temp;
                    s += 2;
                    len -= 2;
                }
            }

            // check line numbers
            if(ii != iin)
                daceSetError(__func__, DACE_INFO, 64);

            // check order and hence number of variables
            unsigned int io = 0;
            for(unsigned int i = 0; i < DACECom.nvmax; i++) io += j[i];
            if(io != io1)
            {
                daceSetError(__func__, DACE_INFO, 65);
                continue;
            }

            // check cutoff order
            if(io > DACECom_t.nocut) continue;

            // actually add the coefficient to the scratch DA
            const unsigned int icc = daceEncode(j);
            if(cc[icc] != 0.0)
                daceSetError(__func__, DACE_INFO, 66);
            cc[icc] += c;
        }

        dacePack(cc, ina);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
        dacefree(cc);
        dacefree(j);
#endif
    }
}

/*! Print a DA object ina to the standard output.
   \param[out] ina Pointer to the DA object to printed
   \sa daceWrite
   \sa daceRead
 */
void dacePrint(const DACEDA *ina)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    if(illa == 0)
    {
        printf("        ALL COEFFICIENTS ZERO\n");
    }
    else
    {
        printf("     I  COEFFICIENT              ORDER EXPONENTS\n");
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
        unsigned int j[DACE_STATIC_NVMAX];
#else
        unsigned int *j = dacecalloc(DACECom.nvmax, sizeof(unsigned int));
#endif
        unsigned int iout = 1;
        for(unsigned int ioa = 0; ioa <= DACECom.nomax; ioa++)
        {
            for(monomial* m = ipoa; m < ipoa+illa; m++)
            {
                if(DACECom.ieo[m->ii] != ioa) continue;
                daceDecode(m->ii, j);
                printf("%6u  %24.16e%6u ", iout, m->cc, ioa);
                for(unsigned int i = 0; i < DACECom.nvmax; i++)
                    printf(" %2u", j[i]);
                printf("\n");
                iout++;
            }
        }
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
        dacefree(j);
#endif
    }

    printf("------------------------------------------------\n");
}


/*! Export a DA object in a binary format.
    The binary data is not supposed to be modified and its format is considered
    internal to the DACE. It is guaranteed that a binary representation can
    be read back into a DA object even with differently initialized settings.
   \param[in] ina The DA object to export
   \param[inout] blob Pointer to memory where to store the data
   \param[inout] size On input contains the size (in bytes) of the memory
    pointed to by blob. On output contains the actual amount of memory used.
   \return Returns 0 is the entire DA object was exported successfully or the
    number of truncated monomials if the object was truncated.
   \note If blob is NULL, the value returned in size is the size (in bytes) of
    the memory needed to store the entire DA object.
*/
unsigned int daceExportBlob(const DACEDA *ina, void *blob, unsigned int *size)
{
    monomial *ipoa; unsigned int ilma, illa;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);

    if(blob == NULL)
    {
        // even if illa == 0 the size we request is still at least that of a daceblob
        *size = sizeof(struct daceblob)+illa*sizeof(extended_monomial);
        if(illa > 0)
            *size -= sizeof(extended_monomial);
        return 0;
    }

    const unsigned int maxsize = *size;
    if(maxsize < sizeof(struct daceblob))
    {
        *size = 0;
        return 1;
    }
    const unsigned int maxnum = (maxsize-sizeof(struct daceblob))/sizeof(extended_monomial)+1;  // maximum number of monomials we can store

    struct daceblob *data = (struct daceblob*) blob;
    data->magic = DACE_BINARY_MAGIC;
    data->no = DACECom.nomax;
    data->nv1 = DACECom.nv1;
    data->nv2 = DACECom.nv2;
    data->len = umin(maxnum, illa);

    for(unsigned int i = 0; i < data->len; i++, ipoa++)
    {
        data->monomials[i].cc = ipoa->cc;
        data->monomials[i].i1 = DACECom.ie1[ipoa->ii];
        data->monomials[i].i2 = DACECom.ie2[ipoa->ii];
    }

    *size = sizeof(struct daceblob)+data->len*sizeof(extended_monomial);
    if(illa > 0)
        *size -= sizeof(extended_monomial);
    return illa - data->len;
}

/*! Determine the total size (in byte) of a DACE blob.
   \param[in] blob Pointer to memory where the data is stored
   \note If called with blob==NULL, the routine will return the minimum size
    of data that must be read in order to determine the total size. A user can
    therefore call this routine twice: first with NULL to determine the size
    of the blob header to read, and then a second time with the header to determine
    the size of the remaining data.
   \return The total size (in bytes) of the DACE blob, or 0 on error (e.g. when
    the data pointed to by blob is not a DACE blob at all). If blob is NULL, the
    minimum size (in bytes) of a blob (i.e. the blob header size) is returned.
*/
unsigned int daceBlobSize(const void *blob)
{
    if(blob == NULL)
    {
        return sizeof(struct daceblob);
    }
    else
    {
        struct daceblob *data = (struct daceblob*)blob;
        if(data->magic != DACE_BINARY_MAGIC)
            return 0;
        unsigned int len = data->len;
        if(len > 0) len--;  // subtract one which is included in the header
        return sizeof(struct daceblob) + len*sizeof(extended_monomial);
    }
}

/*! Import a DA object in a binary format.
   \param[in] blob Pointer to memory where the data is stored
   \param[in] ina The DA object to import into
   \note This routine will silently truncate orders above the currently set
    maximum computation order as well as any extra variables present.
*/
void daceImportBlob(const void *blob, DACEDA *inc)
{
    struct daceblob *data = (struct daceblob*) blob;

    if(data->magic != DACE_BINARY_MAGIC)
    {
        daceSetError(__func__, DACE_ERROR, 31); 
        daceCreateConstant(inc, 0.0);
        return;
    }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double cc[DACE_STATIC_NMMAX] = {0};
#else
    double *cc = dacecalloc(DACECom.nmmax, sizeof(double));
#endif
    const unsigned int nv = data->nv1 + data->nv2;
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    if(nv > DACE_STATIC_NVMAX)
    {
        daceSetError(__func__, DACE_ERROR, 23); 
        daceCreateConstant(inc, 0.0);
        return;
    }
    unsigned int p[DACE_STATIC_NVMAX];
#else
    unsigned int *p = dacecalloc(umax(nv, DACECom.nvmax), sizeof(unsigned int));
#endif
    //bool truncated = false;

    for(unsigned int i = 0; i < data->len; i++)
    {
        unsigned int order = daceDecodeExponents(data->monomials[i].i1, data->no, data->nv1, p);
        order += daceDecodeExponents(data->monomials[i].i2, data->no, data->nv2, p+data->nv1);

        // order of variables outside of current setup
        unsigned int extravar = 0;
        for(unsigned int j = DACECom.nvmax; j < nv; j++)
            extravar += p[j];

        if(order <= DACECom.nomax && extravar == 0)
            cc[daceEncode(p)] = data->monomials[i].cc;
        //else
        //    truncated = true;
    }

    dacePack(cc, inc);

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(cc);
    dacefree(p);
#endif
}
/** @}*/
