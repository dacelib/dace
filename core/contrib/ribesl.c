/* ribesl.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/** \addtogroup DACE Contrib
 *  @{
 */

/// @cond

#include "f2c.h"

/* Subroutine */ int ribesl_(const doublereal *x, const doublereal *alpha, const integer *nb,
	const integer *ize, doublereal *b, integer *ncalc)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal ensig = 1e16;
    static doublereal rtnsig = 1e-4;
    static doublereal enmten = 8.9e-308;
    static doublereal two = 2.;
    static doublereal zero = 0.;
    static doublereal half = .5;
    static doublereal const__ = 1.585;
    static integer nsig = 16;
    static doublereal xlarge = 1e4;
    static doublereal exparg = 709.;
    static doublereal enten = 1e308;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), pow_di(const doublereal *, const integer *), pow_dd(
	    const doublereal *, const doublereal *), exp(doublereal);

    /* Local variables */
    integer k, l, n;
    doublereal p, em, en, sum;
    integer nend, magx;
    doublereal pold;
    integer nbmx;
    doublereal test, empal, halfx, tempa, tempb, tempc, psave, plast, tover,
	    emp2al;
    extern doublereal dgamma_(doublereal *);
    doublereal psavel;
    integer nstart;

/* ------------------------------------------------------------------- */

/*  http://www.netlib.org/specfun/ribesl */

/*  This routine calculates Bessel functions I SUB(N+ALPHA) (X) */
/*  for non-negative argument X, and non-negative order N+ALPHA, */
/*  with or without exponential scaling. */


/* Explanation of variables in the calling sequence */

/* X     - Working precision non-negative real argument for which */
/*         I's or exponentially scaled I's (I*EXP(-X)) */
/*         are to be calculated.  If I's are to be calculated, */
/*         X must be less than EXPARG (see below). */
/* ALPHA - Working precision fractional part of order for which */
/*         I's or exponentially scaled I's (I*EXP(-X)) are */
/*         to be calculated.  0 .LE. ALPHA .LT. 1.0. */
/* NB    - Integer number of functions to be calculated, NB .GT. 0. */
/*         The first function calculated is of order ALPHA, and the */
/*         last is of order (NB - 1 + ALPHA). */
/* IZE   - Integer type.  IZE = 1 if unscaled I's are to calculated, */
/*         and 2 if exponentially scaled I's are to be calculated. */
/* B     - Working precision output vector of length NB.  If the routine */
/*         terminates normally (NCALC=NB), the vector B contains the */
/*         functions I(ALPHA,X) through I(NB-1+ALPHA,X), or the */
/*         corresponding exponentially scaled functions. */
/* NCALC - Integer output variable indicating possible errors. */
/*         Before using the vector B, the user should check that */
/*         NCALC=NB, i.e., all orders have been calculated to */
/*         the desired accuracy.  See error returns below. */


/* ******************************************************************* */
/* ******************************************************************* */

/* Explanation of machine-dependent constants */

/*   beta   = Radix for the floating-point system */
/*   minexp = Smallest representable power of beta */
/*   maxexp = Smallest power of beta that overflows */
/*   it     = Number of bits in the mantissa of a working precision */
/*            variable */
/*   NSIG   = Decimal significance desired.  Should be set to */
/*            INT(LOG10(2)*it+1).  Setting NSIG lower will result */
/*            in decreased accuracy while setting NSIG higher will */
/*            increase CPU time without increasing accuracy.  The */
/*            truncation error is limited to a relative error of */
/*            T=.5*10**(-NSIG). */
/*   ENTEN  = 10.0 ** K, where K is the largest integer such that */
/*            ENTEN is machine-representable in working precision */
/*   ENSIG  = 10.0 ** NSIG */
/*   RTNSIG = 10.0 ** (-K) for the smallest integer K such that */
/*            K .GE. NSIG/4 */
/*   ENMTEN = Smallest ABS(X) such that X/4 does not underflow */
/*   XLARGE = Upper limit on the magnitude of X when IZE=2.  Bear */
/*            in mind that if ABS(X)=N, then at least N iterations */
/*            of the backward recursion will be executed.  The value */
/*            of 10.0 ** 4 is used on every machine. */
/*   EXPARG = Largest working precision argument that the library */
/*            EXP routine can handle and upper limit on the */
/*            magnitude of X when IZE=1; approximately */
/*            LOG(beta**maxexp) */


/*     Approximate values for some important machines are: */

/*                        beta       minexp      maxexp       it */

/*  CRAY-1        (S.P.)    2        -8193        8191        48 */
/*  Cyber 180/855 */
/*    under NOS   (S.P.)    2         -975        1070        48 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (S.P.)    2         -126         128        24 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (D.P.)    2        -1022        1024        53 */
/*  IBM 3033      (D.P.)   16          -65          63        14 */
/*  VAX           (S.P.)    2         -128         127        24 */
/*  VAX D-Format  (D.P.)    2         -128         127        56 */
/*  VAX G-Format  (D.P.)    2        -1024        1023        53 */


/*                        NSIG       ENTEN       ENSIG      RTNSIG */

/* CRAY-1        (S.P.)    15       1.0E+2465   1.0E+15     1.0E-4 */
/* Cyber 180/855 */
/*   under NOS   (S.P.)    15       1.0E+322    1.0E+15     1.0E-4 */
/* IEEE (IBM/XT, */
/*   SUN, etc.)  (S.P.)     8       1.0E+38     1.0E+8      1.0E-2 */
/* IEEE (IBM/XT, */
/*   SUN, etc.)  (D.P.)    16       1.0D+308    1.0D+16     1.0D-4 */
/* IBM 3033      (D.P.)     5       1.0D+75     1.0D+5      1.0D-2 */
/* VAX           (S.P.)     8       1.0E+38     1.0E+8      1.0E-2 */
/* VAX D-Format  (D.P.)    17       1.0D+38     1.0D+17     1.0D-5 */
/* VAX G-Format  (D.P.)    16       1.0D+307    1.0D+16     1.0D-4 */


/*                         ENMTEN      XLARGE   EXPARG */

/* CRAY-1        (S.P.)   1.84E-2466   1.0E+4    5677 */
/* Cyber 180/855 */
/*   under NOS   (S.P.)   1.25E-293    1.0E+4     741 */
/* IEEE (IBM/XT, */
/*   SUN, etc.)  (S.P.)   4.70E-38     1.0E+4      88 */
/* IEEE (IBM/XT, */
/*   SUN, etc.)  (D.P.)   8.90D-308    1.0D+4     709 */
/* IBM 3033      (D.P.)   2.16D-78     1.0D+4     174 */
/* VAX           (S.P.)   1.17E-38     1.0E+4      88 */
/* VAX D-Format  (D.P.)   1.17D-38     1.0D+4      88 */
/* VAX G-Format  (D.P.)   2.22D-308    1.0D+4     709 */

/* ******************************************************************* */
/* ******************************************************************* */

/* Error returns */

/*  In case of an error,  NCALC .NE. NB, and not all I's are */
/*  calculated to the desired accuracy. */

/*  NCALC .LT. 0:  An argument is out of range. For example, */
/*     NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE. EXPARG. */
/*     In this case, the B-vector is not calculated, and NCALC is */
/*     set to MIN0(NB,0)-1 so that NCALC .NE. NB. */

/*  NB .GT. NCALC .GT. 0: Not all requested function values could */
/*     be calculated accurately.  This usually occurs because NB is */
/*     much larger than ABS(X).  In this case, B(N) is calculated */
/*     to the desired accuracy for N .LE. NCALC, but precision */
/*     is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish */
/*     for N .GT. NCALC (because it is too small to be represented), */
/*     and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K */
/*     significant figures of B(N) can be trusted. */


/* Intrinsic functions required are: */

/*     DBLE, EXP, DGAMMA, GAMMA, INT, MAX, MIN, REAL, SQRT */


/* Acknowledgement */

/*  This program is based on a program written by David J. */
/*  Sookne (2) that computes values of the Bessel functions J or */
/*  I of real argument and integer order.  Modifications include */
/*  the restriction of the computation to the I Bessel function */
/*  of non-negative real argument, the extension of the computation */
/*  to arbitrary positive order, the inclusion of optional */
/*  exponential scaling, and the elimination of most underflow. */
/*  An earlier version was published in (3). */

/* References: "A Note on Backward Recurrence Algorithms," Olver, */
/*              F. W. J., and Sookne, D. J., Math. Comp. 26, 1972, */
/*              pp 941-947. */

/*             "Bessel Functions of Real Argument and Integer Order," */
/*              Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp */
/*              125-132. */

/*             "ALGORITHM 597, Sequence of Modified Bessel Functions */
/*              of the First Kind," Cody, W. J., Trans. Math. Soft., */
/*              1983, pp. 242-245. */

/*  Latest modification: May 30, 1989 */

/*  Modified by: W. J. Cody and L. Stoltz */
/*               Applied Mathematics Division */
/*               Argonne National Laboratory */
/*               Argonne, IL  60439 */

/* ------------------------------------------------------------------- */
/* S    REAL              GAMMA, */
/* ------------------------------------------------------------------- */
/*  Mathematical constants */
/* ------------------------------------------------------------------- */
/* S    DATA ONE,TWO,ZERO,HALF,CONST/1.0E0,2.0E0,0.0E0,0.5E0,1.585E0/ */
    /* Parameter adjustments */
    --b;

    /* Function Body */
/* ------------------------------------------------------------------- */
/*  Machine-dependent parameters */
/* ------------------------------------------------------------------- */
/* S    DATA NSIG,XLARGE,EXPARG /8,1.0E4,88.0E0/ */
/* S    DATA ENTEN,ENSIG,RTNSIG/1.0E38,1.0E8,1.0E-2/ */
/* S    DATA ENMTEN/4.7E-38/ */
/* ------------------------------------------------------------------- */
/*  Statement functions for conversion */
/* ------------------------------------------------------------------- */
/* S    CONV(N) = REAL(N) */
/* S    FUNC(X) = GAMMA(X) */
/* ------------------------------------------------------------------- */
/* Check for X, NB, OR IZE out of range. */
/* ------------------------------------------------------------------- */
    if (*nb > 0 && *x >= zero && *alpha >= zero && *alpha < one && ((*ize == 1
	    && *x <= exparg) || (*ize == 2 && *x <= xlarge))) {
/* ------------------------------------------------------------------- */
/* Use 2-term ascending series for small X */
/* ------------------------------------------------------------------- */
	*ncalc = *nb;
	magx = (integer) (*x);
	if (*x >= rtnsig) {
/* ------------------------------------------------------------------- */
/* Initialize the forward sweep, the P-sequence of Olver */
/* ------------------------------------------------------------------- */
	    nbmx = *nb - magx;
	    n = magx + 1;
	    i__1 = n + n;
	    en = (doublereal) i__1 + (*alpha + *alpha);
	    plast = one;
	    p = en / *x;
/* ------------------------------------------------------------------- */
/* Calculate general significance test */
/* ------------------------------------------------------------------- */
	    test = ensig + ensig;
	    if (magx << 1 > nsig * 5) {
		test = sqrt(test * p);
	    } else {
		test /= pow_di(&const__, &magx);
	    }
	    if (nbmx >= 3) {
/* ------------------------------------------------------------------- */
/* Calculate P-sequence until N = NB-1.  Check for possible overflow. */
/* ------------------------------------------------------------------- */
		tover = enten / ensig;
		nstart = magx + 2;
		nend = *nb - 1;
		i__1 = nend;
		for (k = nstart; k <= i__1; ++k) {
		    n = k;
		    en += two;
		    pold = plast;
		    plast = p;
		    p = en * plast / *x + pold;
		    if (p > tover) {
/* ------------------------------------------------------------------- */
/* To avoid overflow, divide P-sequence by TOVER.  Calculate */
/* P-sequence until ABS(P) .GT. 1. */
/* ------------------------------------------------------------------- */
			tover = enten;
			p /= tover;
			plast /= tover;
			psave = p;
			psavel = plast;
			nstart = n + 1;
L60:
			++n;
			en += two;
			pold = plast;
			plast = p;
			p = en * plast / *x + pold;
			if (p <= one) {
			    goto L60;
			}
			tempb = en / *x;
/* ------------------------------------------------------------------- */
/* Calculate backward test, and find NCALC, the highest N */
/* such that the test is passed. */
/* ------------------------------------------------------------------- */
			test = pold * plast / ensig;
			test *= half - half / (tempb * tempb);
			p = plast * tover;
			--n;
			en -= two;
			nend = min(*nb,n);
			i__2 = nend;
			for (l = nstart; l <= i__2; ++l) {
			    *ncalc = l;
			    pold = psavel;
			    psavel = psave;
			    psave = en * psavel / *x + pold;
			    if (psave * psavel > test) {
				goto L90;
			    }
/* L80: */
			}
			*ncalc = nend + 1;
L90:
			--(*ncalc);
			goto L120;
		    }
/* L100: */
		}
		n = nend;
		i__1 = n + n;
		en = (doublereal) i__1 + (*alpha + *alpha);
/* ------------------------------------------------------------------- */
/* Calculate special significance test for NBMX .GT. 2. */
/* ------------------------------------------------------------------- */
/* Computing MAX */
		d__1 = test, d__2 = sqrt(plast * ensig) * sqrt(p + p);
		test = max(d__1,d__2);
	    }
/* ------------------------------------------------------------------- */
/* Calculate P-sequence until significance test passed. */
/* ------------------------------------------------------------------- */
L110:
	    ++n;
	    en += two;
	    pold = plast;
	    plast = p;
	    p = en * plast / *x + pold;
	    if (p < test) {
		goto L110;
	    }
/* ------------------------------------------------------------------- */
/* Initialize the backward recursion and the normalization sum. */
/* ------------------------------------------------------------------- */
L120:
	    ++n;
	    en += two;
	    tempb = zero;
	    tempa = one / p;
	    em = (doublereal) n - one;
	    empal = em + *alpha;
	    emp2al = em - one + (*alpha + *alpha);
	    sum = tempa * empal * emp2al / em;
	    nend = n - *nb;
	    if (nend < 0) {
/* ------------------------------------------------------------------- */
/* N .LT. NB, so store B(N) and set higher orders to zero. */
/* ------------------------------------------------------------------- */
		b[n] = tempa;
		nend = -nend;
		i__1 = nend;
		for (l = 1; l <= i__1; ++l) {
/* L130: */
		    b[n + l] = zero;
		}
	    } else {
		if (nend > 0) {
/* ------------------------------------------------------------------- */
/* Recur backward via difference equation, calculating (but */
/* not storing) B(N), until N = NB. */
/* ------------------------------------------------------------------- */
		    i__1 = nend;
		    for (l = 1; l <= i__1; ++l) {
			--n;
			en -= two;
			tempc = tempb;
			tempb = tempa;
			tempa = en * tempb / *x + tempc;
			em -= one;
			emp2al -= one;
			if (n == 1) {
			    goto L150;
			}
			if (n == 2) {
			    emp2al = one;
			}
			empal -= one;
			sum = (sum + tempa * empal) * emp2al / em;
/* L140: */
		    }
		}
/* ------------------------------------------------------------------- */
/* Store B(NB) */
/* ------------------------------------------------------------------- */
L150:
		b[n] = tempa;
		if (*nb <= 1) {
		    sum = sum + sum + tempa;
		    goto L230;
		}
/* ------------------------------------------------------------------- */
/* Calculate and Store B(NB-1) */
/* ------------------------------------------------------------------- */
		--n;
		en -= two;
		b[n] = en * tempa / *x + tempb;
		if (n == 1) {
		    goto L220;
		}
		em -= one;
		emp2al -= one;
		if (n == 2) {
		    emp2al = one;
		}
		empal -= one;
		sum = (sum + b[n] * empal) * emp2al / em;
	    }
	    nend = n - 2;
	    if (nend > 0) {
/* ------------------------------------------------------------------- */
/* Calculate via difference equation and store B(N), until N = 2. */
/* ------------------------------------------------------------------- */
		i__1 = nend;
		for (l = 1; l <= i__1; ++l) {
		    --n;
		    en -= two;
		    b[n] = en * b[n + 1] / *x + b[n + 2];
		    em -= one;
		    emp2al -= one;
		    if (n == 2) {
			emp2al = one;
		    }
		    empal -= one;
		    sum = (sum + b[n] * empal) * emp2al / em;
/* L200: */
		}
	    }
/* ------------------------------------------------------------------- */
/* Calculate B(1) */
/* ------------------------------------------------------------------- */
	    b[1] = two * empal * b[2] / *x + b[3];
L220:
	    sum = sum + sum + b[1];
/* ------------------------------------------------------------------- */
/* Normalize.  Divide all B(N) by sum. */
/* ------------------------------------------------------------------- */
L230:
	    if (*alpha != zero) {
		d__1 = one + *alpha;
		d__2 = *x * half;
		d__3 = -(*alpha);
		sum = sum * dgamma_(&d__1) * pow_dd(&d__2, &d__3);
	    }
	    if (*ize == 1) {
		sum *= exp(-(*x));
	    }
	    tempa = enmten;
	    if (sum > one) {
		tempa *= sum;
	    }
	    i__1 = *nb;
	    for (n = 1; n <= i__1; ++n) {
		if (b[n] < tempa) {
		    b[n] = zero;
		}
		b[n] /= sum;
/* L260: */
	    }
	    return 0;
/* ------------------------------------------------------------------- */
/* Two-term ascending series for small X. */
/* ------------------------------------------------------------------- */
	} else {
	    tempa = one;
	    empal = one + *alpha;
	    halfx = zero;
	    if (*x > enmten) {
		halfx = half * *x;
	    }
	    if (*alpha != zero) {
		tempa = pow_dd(&halfx, alpha) / dgamma_(&empal);
	    }
	    if (*ize == 2) {
		tempa *= exp(-(*x));
	    }
	    tempb = zero;
	    if (*x + one > one) {
		tempb = halfx * halfx;
	    }
	    b[1] = tempa + tempa * tempb / empal;
	    if (*x != zero && b[1] == zero) {
		*ncalc = 0;
	    }
	    if (*nb > 1) {
		if (*x == zero) {
		    i__1 = *nb;
		    for (n = 2; n <= i__1; ++n) {
			b[n] = zero;
/* L310: */
		    }
		} else {
/* ------------------------------------------------------------------- */
/* Calculate higher-order functions. */
/* ------------------------------------------------------------------- */
		    tempc = halfx;
		    tover = (enmten + enmten) / *x;
		    if (tempb != zero) {
			tover = enmten / tempb;
		    }
		    i__1 = *nb;
		    for (n = 2; n <= i__1; ++n) {
			tempa /= empal;
			empal += one;
			tempa *= tempc;
			if (tempa <= tover * empal) {
			    tempa = zero;
			}
			b[n] = tempa + tempa * tempb / empal;
			if (b[n] == zero && *ncalc > n) {
			    *ncalc = n - 1;
			}
/* L340: */
		    }
		}
	    }
	}
    } else {
	*ncalc = min(*nb,0) - 1;
    }
    return 0;
/* ---------- Last line of RIBESL ---------- */
} /* ribesl_ */

/// @endcond
/** @}*/