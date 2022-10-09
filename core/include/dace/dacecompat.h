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
 *  dacecompat.h
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/*
    Compatibility file mapping old unreadable names to new DACE function names.
    Where possible, function names are #define'ed to map to new names.
    Where necessary, small shims with the old interface are provided in dacecompat.c.
*/

/** \addtogroup DACE Core
 *  @{
 */
#ifndef DINAMICA_DACECOMPAT_H_
#define DINAMICA_DACECOMPAT_H_

#ifdef __cplusplus
    extern "C" {
#endif  /* __cplusplus */

#define daceini(no, nv) daceInitialize(no, nv)
#define dacever(imaj, imin, icos) daceGetVersion(imaj, imin, icos)
#define daceseteps(deps) daceSetEpsilon(deps)
#define dacegeteps() daceGetEpsilon()
#define dacegetepsmac() daceGetMachineEpsilon()
#define dacegetnomax() daceGetMaxOrder()
#define dacegetnvmax() daceGetMaxVariables()
#define dacegetnmmax() daceGetMaxMonomials()
#define dacegetnot() daceGetTruncationOrder()
#define dacesetnot(fnot) daceSetTruncationOrder(fnot)

//#define dacegeterr() daceGetError()        // do not redefine this routine, it causes problems with the inlining of a function of the same name we do in C++ interface!
#define dacegetxerr() daceGetErrorX()
#define dacegetyyerr() daceGetErrorYY()
#define daceclrerr() daceClearError()

#define daceall(inc, len) daceAllocateDA(inc, len)
#define dacedal(inc) daceFreeDA(ARG(inc))
#define dacememdump(iunit) daceMemoryDump()

#define dacevar(ina, i, ckon) daceCreateVariable(ARG(ina), i, ckon)
#define dacecoef(ina, jj, ckon) daceCreateMonomial(ARG(ina), jj, ckon)
#define dacecon(ina, ckon) daceCreateConstant(ARG(ina), ckon)
#define daceran(ina, cm) daceCreateRandom(ARG(ina), cm)

#define daceconst(ina) daceGetConstant(ARG(ina))
#define dacelinear(ina, c) daceGetLinear(ARG(ina), c)
#define dacepok(ina, jj, cjj) daceSetCoefficient(ARG(ina), jj, cjj)
DACE_API void dacepek(const DACEDA REF(ina), const unsigned int jj[], double REF(cjj));
#define dacelist(ina, npos, jj, cjj) daceGetCoefficientAt(ARG(ina), npos, jj, cjj)
DACE_API void dacesize(const DACEDA REF(ina), unsigned int REF(size));

#define dacecop(ina, inb) daceCopy(ARG(ina), ARG(inb))

#define dacetrim(ina, imin, imax, inc) daceTrim(ARG(ina), imin, imax, ARG(inc))

DACE_API void daceabs(const DACEDA REF(ina), double REF(anorm));
DACE_API void dacenorm(const DACEDA REF(ina), const unsigned int ityp, double REF(anorm));
#define daceonorm(ina, ivar, ityp, onorm) daceOrderedNorm(ARG(ina), ivar, ityp, onorm)
#define daceest(ina, ivar, ityp, c, nc) daceEstimate(ARG(ina), ivar, ityp, c, NULL, nc)
#define dacebound(ina, alo, aup) daceGetBounds(ARG(ina), alo, aup)

#define daceplug(ina, nvar, val, inc) daceEvalVariable(ARG(ina), nvar, val, ARG(inc))
DACE_API void dacetree(const DACEDA das[], const unsigned int count, double ac[], unsigned int REF(nterm), unsigned int REF(nvar), unsigned int REF(nord));
#define dacewrite(ina, strs, nstrs) daceWrite(ARG(ina), strs, nstrs)
#define daceread(ina, strs, nstrs) daceRead(ARG(ina), strs, nstrs)

#define daceadd(ina, inb, inc) daceAdd(ARG(ina), ARG(inb), ARG(inc))
#define dacesub(ina, inb, inc) daceSubtract(ARG(ina), ARG(inb), ARG(inc))
#define dacemul(ina, inb, inc) daceMultiply(ARG(ina), ARG(inb), ARG(inc))
#define dacediv(ina, inb, inc) daceDivide(ARG(ina), ARG(inb), ARG(inc))
#define dacesqr(ina, inb) daceSquare(ARG(ina), ARG(inb))
#define dacecadd(ina, ckon, inb) daceAddDouble(ARG(ina), ckon, ARG(inb))
#define dacecsub(ina, ckon, inb) daceDoubleSubtract(ARG(ina), ckon, ARG(inb))
#define dacesubc(ina, ckon, inb) daceSubtractDouble(ARG(ina), ckon, ARG(inb))
#define dacecmul(ina, ckon, inb) daceMultiplyDouble(ARG(ina), ckon, ARG(inb))
#define dacecdiv(ina, ckon, inb) daceDoubleDivide(ARG(ina), ckon, ARG(inb))
#define dacedivc(ina, ckon, inb) daceDivideDouble(ARG(ina), ckon, ARG(inb))
#define daceder(idif, ina, inc) daceDifferentiate(idif, ARG(ina), ARG(inc))
#define daceint(iint, ina, inc) daceIntegrate(iint, ARG(ina), ARG(inc))

#define dacetrunc(ina, inc) daceTruncate(ARG(ina), ARG(inc))
#define daceround(ina, inc) daceRound(ARG(ina), ARG(inc))
#define dacemod(ina, p, inc) daceModulo(ARG(ina), p, ARG(inc))
#define dacepow(ina, np, inc) dacePower(ARG(ina), np, ARG(inc))
#define daceroot(ina, np, inc) daceRoot(ARG(ina), np, ARG(inc))
#define daceminv(ina, inc) daceMultiplicativeInverse(ARG(ina), ARG(inc))
#define dacesqrt(ina, inc) daceSquareRoot(ARG(ina), ARG(inc))
#define daceisrt(ina, inc) daceInverseSquareRoot(ARG(ina), ARG(inc))
#define daceexp(ina, inc) daceExponential(ARG(ina), ARG(inc))
#define dacelog(ina, inc) daceLogarithm(ARG(ina), ARG(inc))
#define dacelogb(ina, b, inc) daceLogarithmBase(ARG(ina), b, ARG(inc))
#define dacesin(ina, inc) daceSine(ARG(ina), ARG(inc))
#define dacecos(ina, inc) daceCosine(ARG(ina), ARG(inc))
#define dacetan(ina, inc) daceTangent(ARG(ina), ARG(inc))
#define daceasin(ina, inc) daceArcSine(ARG(ina), ARG(inc))
#define daceacos(ina, inc) daceArcCosine(ARG(ina), ARG(inc))
#define daceatan(ina, inc) daceArcTangent(ARG(ina), ARG(inc))
#define daceatan2(ina, inb, inc) daceArcTangent2(ARG(ina), ARG(inb), ARG(inc))
#define dacesinh(ina, inc) daceHyperbolicSine(ARG(ina), ARG(inc))
#define dacecosh(ina, inc) daceHyperbolicCosine(ARG(ina), ARG(inc))
#define dacetanh(ina, inc) daceHyperbolicTangent(ARG(ina), ARG(inc))
#define daceasinh(ina, inc) daceHyperbolicArcSine(ARG(ina), ARG(inc))
#define daceacosh(ina, inc) daceHyperbolicArcCosine(ARG(ina), ARG(inc))
#define daceatanh(ina, inc) daceHyperbolicArcTangent(ARG(ina), ARG(inc))

#ifdef __cplusplus
    }   // extern "C"
#endif  /* __cplusplus */
/** @}*/
#endif /* DINAMICA_DACECOMPAT_H_ */
