/* rjbesl.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int rjbesl_(const doublereal *x, const doublereal *alpha, const integer *nb,
	doublereal *b, integer *ncalc)
{
    /* Initialized data */

    static doublereal pi2 = .636619772367581343075535;
    static doublereal four = 4.;
    static doublereal twofiv = 25.;
    static doublereal one30 = 130.;
    static doublereal three5 = 35.;
    static doublereal enten = 1e308;
    static doublereal ensig = 1e16;
    static doublereal rtnsig = 1e-4;
    static doublereal enmten = 8.9e-308;
    static doublereal xlarge = 1e4;
    static doublereal fact[25] = { 1.,1.,2.,6.,24.,120.,720.,5040.,40320.,
	    362880.,3628800.,39916800.,479001600.,6227020800.,87178291200.,
	    1.307674368e12,2.0922789888e13,3.55687428096e14,6.402373705728e15,
	    1.21645100408832e17,2.43290200817664e18,5.109094217170944e19,
	    1.12400072777760768e21,2.585201673888497664e22,
	    6.2044840173323943936e23 };
    static doublereal twopi1 = 6.28125;
    static doublereal twopi2 = .001935307179586476925286767;
    static doublereal zero = 0.;
    static doublereal eighth = .125;
    static doublereal half = .5;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(const doublereal *, const doublereal *), sqrt(doublereal), d_int(
	    const doublereal *), sin(doublereal), cos(doublereal);

    /* Local variables */
    integer i__, j, k, l, m, n;
    doublereal p, s, t, z__, t1, em, en, xc, xk, xm, gnu, xin, sum, capp;
    integer nend;
    doublereal capq;
    integer magx;
    doublereal pold;
    integer nbmx;
    doublereal vcos, test, vsin, alpem, halfx, tempa, tempb, tempc, psave,
	    plast, tover, alp2em;
    extern doublereal dgamma_(const doublereal *);
    doublereal psavel;
    integer nstart;

/* --------------------------------------------------------------------- */
/* This routine calculates Bessel functions J sub(N+ALPHA) (X) */
/*   for non-negative argument X, and non-negative order N+ALPHA. */


/*  Explanation of variables in the calling sequence. */

/*   X     - working precision non-negative real argument for which */
/*           J's are to be calculated. */
/*   ALPHA - working precision fractional part of order for which */
/*           J's or exponentially scaled J'r (J*exp(X)) are */
/*           to be calculated.  0 <= ALPHA < 1.0. */
/*   NB  - integer number of functions to be calculated, NB > 0. */
/*           The first function calculated is of order ALPHA, and the */
/*           last is of order (NB - 1 + ALPHA). */
/*   B  - working precision output vector of length NB.  If RJBESL */
/*           terminates normally (NCALC=NB), the vector B contains the */
/*           functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the */
/*           corresponding exponentially scaled functions. */
/*   NCALC - integer output variable indicating possible errors. */
/*           Before using the vector B, the user should check that */
/*           NCALC=NB, i.e., all orders have been calculated to */
/*           the desired accuracy.  See Error Returns below. */


/* ******************************************************************* */
/* ******************************************************************* */

/*  Explanation of machine-dependent constants */

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
/*   XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N, */
/*            then at least N iterations of the backward recursion */
/*            will be executed.  The value of 10.0 ** 4 is used on */
/*            every machine. */


/*     Approximate values for some important machines are: */


/*                            it    NSIG    ENTEN       ENSIG */

/*   CRAY-1        (S.P.)     48     15    1.0E+2465   1.0E+15 */
/*   Cyber 180/855 */
/*     under NOS   (S.P.)     48     15    1.0E+322    1.0E+15 */
/*   IEEE (IBM/XT, */
/*     SUN, etc.)  (S.P.)     24      8    1.0E+38     1.0E+8 */
/*   IEEE (IBM/XT, */
/*     SUN, etc.)  (D.P.)     53     16    1.0D+308    1.0D+16 */
/*   IBM 3033      (D.P.)     14      5    1.0D+75     1.0D+5 */
/*   VAX           (S.P.)     24      8    1.0E+38     1.0E+8 */
/*   VAX D-Format  (D.P.)     56     17    1.0D+38     1.0D+17 */
/*   VAX G-Format  (D.P.)     53     16    1.0D+307    1.0D+16 */


/*                           RTNSIG      ENMTEN      XLARGE */

/*   CRAY-1        (S.P.)    1.0E-4    1.84E-2466   1.0E+4 */
/*   Cyber 180/855 */
/*     under NOS   (S.P.)    1.0E-4    1.25E-293    1.0E+4 */
/*   IEEE (IBM/XT, */
/*     SUN, etc.)  (S.P.)    1.0E-2    4.70E-38     1.0E+4 */
/*   IEEE (IBM/XT, */
/*     SUN, etc.)  (D.P.)    1.0E-4    8.90D-308    1.0D+4 */
/*   IBM 3033      (D.P.)    1.0E-2    2.16D-78     1.0D+4 */
/*   VAX           (S.P.)    1.0E-2    1.17E-38     1.0E+4 */
/*   VAX D-Format  (D.P.)    1.0E-5    1.17D-38     1.0D+4 */
/*   VAX G-Format  (D.P.)    1.0E-4    2.22D-308    1.0D+4 */

/* ******************************************************************* */
/* ******************************************************************* */

/*  Error returns */

/*    In case of an error,  NCALC .NE. NB, and not all J's are */
/*    calculated to the desired accuracy. */

/*    NCALC .LT. 0:  An argument is out of range. For example, */
/*       NBES .LE. 0, ALPHA .LT. 0 or .GT. 1, or X is too large. */
/*       In this case, B(1) is set to zero, the remainder of the */
/*       B-vector is not calculated, and NCALC is set to */
/*       MIN(NB,0)-1 so that NCALC .NE. NB. */

/*    NB .GT. NCALC .GT. 0: Not all requested function values could */
/*       be calculated accurately.  This usually occurs because NB is */
/*       much larger than ABS(X).  In this case, B(N) is calculated */
/*       to the desired accuracy for N .LE. NCALC, but precision */
/*       is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish */
/*       for N .GT. NCALC (because it is too small to be represented), */
/*       and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K */
/*       significant figures of B(N) can be trusted. */


/*  Intrinsic and other functions required are: */

/*     ABS, AINT, COS, DBLE, GAMMA (or DGAMMA), INT, MAX, MIN, */

/*     REAL, SIN, SQRT */


/*  Acknowledgement */

/*   This program is based on a program written by David J. Sookne */
/*   (2) that computes values of the Bessel functions J or I of real */
/*   argument and integer order.  Modifications include the restriction */
/*   of the computation to the J Bessel function of non-negative real */
/*   argument, the extension of the computation to arbitrary positive */
/*   order, and the elimination of most underflow. */

/*  References: "A Note on Backward Recurrence Algorithms," Olver, */
/*               F. W. J., and Sookne, D. J., Math. Comp. 26, 1972, */
/*               pp 941-947. */

/*              "Bessel Functions of Real Argument and Integer Order," */
/*               Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp */
/*               125-132. */

/*  Latest modification: March 19, 1990 */

/*  Author: W. J. Cody */
/*          Applied Mathematics Division */
/*          Argonne National Laboratory */
/*          Argonne, IL  60439 */

/* --------------------------------------------------------------------- */
/* S    REAL               GAMMA, */
/* --------------------------------------------------------------------- */
/*  Mathematical constants */

/*   PI2    - 2 / PI */
/*   TWOPI1 - first few significant digits of 2 * PI */
/*   TWOPI2 - (2*PI - TWOPI) to working precision, i.e., */
/*            TWOPI1 + TWOPI2 = 2 * PI to extra precision. */
/* --------------------------------------------------------------------- */
/* S    DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535E0,6.28125E0, */
/* S   1 1.935307179586476925286767E-3/ */
/* S    DATA ZERO, EIGHTH, HALF, ONE /0.0E0,0.125E0,0.5E0,1.0E0/ */
/* S    DATA TWO, THREE, FOUR, TWOFIV /2.0E0,3.0E0,4.0E0,25.0E0/ */
/* S    DATA ONE30, THREE5 /130.0E0,35.0E0/ */
    /* Parameter adjustments */
    --b;

    /* Function Body */
/* --------------------------------------------------------------------- */
/*  Machine-dependent parameters */
/* --------------------------------------------------------------------- */
/* S    DATA ENTEN, ENSIG, RTNSIG /1.0E38,1.0E8,1.0E-2/ */
/* S    DATA ENMTEN, XLARGE /1.2E-37,1.0E4/ */
/* --------------------------------------------------------------------- */
/*     Factorial(N) */
/* --------------------------------------------------------------------- */
/* S    DATA FACT /1.0E0,1.0E0,2.0E0,6.0E0,24.0E0,1.2E2,7.2E2,5.04E3, */
/* S   1 4.032E4,3.6288E5,3.6288E6,3.99168E7,4.790016E8,6.2270208E9, */
/* S   2 8.71782912E10,1.307674368E12,2.0922789888E13,3.55687428096E14, */
/* S   3 6.402373705728E15,1.21645100408832E17,2.43290200817664E18, */
/* S   4 5.109094217170944E19,1.12400072777760768E21, */
/* S   5 2.585201673888497664E22,6.2044840173323943936E23/ */
/* --------------------------------------------------------------------- */
/* Statement functions for conversion and the gamma function. */
/* --------------------------------------------------------------------- */
/* S    CONV(I) = REAL(I) */
/* S    FUNC(X) = GAMMA(X) */
/* --------------------------------------------------------------------- */
/* Check for out of range arguments. */
/* --------------------------------------------------------------------- */
    magx = (integer) (*x);
    if (*nb > 0 && *x >= zero && *x <= xlarge && *alpha >= zero && *alpha <
	    one) {
/* --------------------------------------------------------------------- */
/* Initialize result array to zero. */
/* --------------------------------------------------------------------- */
	*ncalc = *nb;
	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[i__] = zero;
/* L20: */
	}
/* --------------------------------------------------------------------- */
/* Branch to use 2-term ascending series for small X and asymptotic */
/* form for large X when NB is not too large. */
/* --------------------------------------------------------------------- */
	if (*x < rtnsig) {
/* --------------------------------------------------------------------- */
/* Two-term ascending series for small X. */
/* --------------------------------------------------------------------- */
	    tempa = one;
	    alpem = one + *alpha;
	    halfx = zero;
	    if (*x > enmten) {
		halfx = half * *x;
	    }
	    if (*alpha != zero) {
		tempa = pow_dd(&halfx, alpha) / (*alpha * dgamma_(alpha));
	    }
	    tempb = zero;
	    if (*x + one > one) {
		tempb = -halfx * halfx;
	    }
	    b[1] = tempa + tempa * tempb / alpem;
	    if (*x != zero && b[1] == zero) {
		*ncalc = 0;
	    }
	    if (*nb != 1) {
		if (*x <= zero) {
		    i__1 = *nb;
		    for (n = 2; n <= i__1; ++n) {
			b[n] = zero;
/* L30: */
		    }
		} else {
/* --------------------------------------------------------------------- */
/* Calculate higher order functions. */
/* --------------------------------------------------------------------- */
		    tempc = halfx;
		    tover = (enmten + enmten) / *x;
		    if (tempb != zero) {
			tover = enmten / tempb;
		    }
		    i__1 = *nb;
		    for (n = 2; n <= i__1; ++n) {
			tempa /= alpem;
			alpem += one;
			tempa *= tempc;
			if (tempa <= tover * alpem) {
			    tempa = zero;
			}
			b[n] = tempa + tempa * tempb / alpem;
			if (b[n] == zero && *ncalc > n) {
			    *ncalc = n - 1;
			}
/* L50: */
		    }
		}
	    }
	} else if (*x > twofiv && *nb <= magx + 1) {
/* --------------------------------------------------------------------- */
/* Asymptotic series for X .GT. 21.0. */
/* --------------------------------------------------------------------- */
	    xc = sqrt(pi2 / *x);
/* Computing 2nd power */
	    d__1 = eighth / *x;
	    xin = d__1 * d__1;
	    m = 11;
	    if (*x >= three5) {
		m = 8;
	    }
	    if (*x >= one30) {
		m = 4;
	    }
	    xm = four * (doublereal) m;
/* --------------------------------------------------------------------- */
/* Argument reduction for SIN and COS routines. */
/* --------------------------------------------------------------------- */
	    d__1 = *x / (twopi1 + twopi2) + half;
	    t = d_int(&d__1);
	    z__ = *x - t * twopi1 - t * twopi2 - (*alpha + half) / pi2;
	    vsin = sin(z__);
	    vcos = cos(z__);
	    gnu = *alpha + *alpha;
	    for (i__ = 1; i__ <= 2; ++i__) {
		s = (xm - one - gnu) * (xm - one + gnu) * xin * half;
		t = (gnu - (xm - three)) * (gnu + (xm - three));
		capp = s * t / fact[m * 2];
		t1 = (gnu - (xm + one)) * (gnu + (xm + one));
		capq = s * t1 / fact[(m << 1) + 1];
		xk = xm;
		k = m + m;
		t1 = t;
		i__1 = m;
		for (j = 2; j <= i__1; ++j) {
		    xk -= four;
		    s = (xk - one - gnu) * (xk - one + gnu);
		    t = (gnu - (xk - three)) * (gnu + (xk - three));
		    capp = (capp + one / fact[k - 2]) * s * t * xin;
		    capq = (capq + one / fact[k - 1]) * s * t1 * xin;
		    k += -2;
		    t1 = t;
/* L70: */
		}
		capp += one;
		capq = (capq + one) * (gnu * gnu - one) * (eighth / *x);
		b[i__] = xc * (capp * vcos - capq * vsin);
		if (*nb == 1) {
		    goto L300;
		}
		t = vsin;
		vsin = -vcos;
		vcos = t;
		gnu += two;
/* L80: */
	    }
/* --------------------------------------------------------------------- */
/* If  NB .GT. 2, compute J(X,ORDER+I)  I = 2, NB-1 */
/* --------------------------------------------------------------------- */
	    if (*nb > 2) {
		gnu = *alpha + *alpha + two;
		i__1 = *nb;
		for (j = 3; j <= i__1; ++j) {
		    b[j] = gnu * b[j - 1] / *x - b[j - 2];
		    gnu += two;
/* L90: */
		}
	    }
/* --------------------------------------------------------------------- */
/* Use recurrence to generate results.  First initialize the */
/* calculation of P*S. */
/* --------------------------------------------------------------------- */
	} else {
	    nbmx = *nb - magx;
	    n = magx + 1;
	    i__1 = n + n;
	    en = (doublereal) i__1 + (*alpha + *alpha);
	    plast = one;
	    p = en / *x;
/* --------------------------------------------------------------------- */
/* Calculate general significance test. */
/* --------------------------------------------------------------------- */
	    test = ensig + ensig;
	    if (nbmx >= 3) {
/* --------------------------------------------------------------------- */
/* Calculate P*S until N = NB-1.  Check for possible overflow. */
/* --------------------------------------------------------------------- */
		tover = enten / ensig;
		nstart = magx + 2;
		nend = *nb - 1;
		i__1 = nstart + nstart;
		en = (doublereal) i__1 - two + (*alpha + *alpha);
		i__1 = nend;
		for (k = nstart; k <= i__1; ++k) {
		    n = k;
		    en += two;
		    pold = plast;
		    plast = p;
		    p = en * plast / *x - pold;
		    if (p > tover) {
/* --------------------------------------------------------------------- */
/* To avoid overflow, divide P*S by TOVER.  Calculate P*S until */
/* ABS(P) .GT. 1. */
/* --------------------------------------------------------------------- */
			tover = enten;
			p /= tover;
			plast /= tover;
			psave = p;
			psavel = plast;
			nstart = n + 1;
L100:
			++n;
			en += two;
			pold = plast;
			plast = p;
			p = en * plast / *x - pold;
			if (p <= one) {
			    goto L100;
			}
			tempb = en / *x;
/* --------------------------------------------------------------------- */
/* Calculate backward test and find NCALC, the highest N such that */
/* the test is passed. */
/* --------------------------------------------------------------------- */
			test = pold * plast * (half - half / (tempb * tempb));
			test /= ensig;
			p = plast * tover;
			--n;
			en -= two;
			nend = min(*nb,n);
			i__2 = nend;
			for (l = nstart; l <= i__2; ++l) {
			    pold = psavel;
			    psavel = psave;
			    psave = en * psavel / *x - pold;
			    if (psave * psavel > test) {
				*ncalc = l - 1;
				goto L190;
			    }
/* L110: */
			}
			*ncalc = nend;
			goto L190;
		    }
/* L130: */
		}
		n = nend;
		i__1 = n + n;
		en = (doublereal) i__1 + (*alpha + *alpha);
/* --------------------------------------------------------------------- */
/* Calculate special significance test for NBMX .GT. 2. */
/* --------------------------------------------------------------------- */
/* Computing MAX */
		d__1 = test, d__2 = sqrt(plast * ensig) * sqrt(p + p);
		test = max(d__1,d__2);
	    }
/* --------------------------------------------------------------------- */
/* Calculate P*S until significance test passes. */
/* --------------------------------------------------------------------- */
L140:
	    ++n;
	    en += two;
	    pold = plast;
	    plast = p;
	    p = en * plast / *x - pold;
	    if (p < test) {
		goto L140;
	    }
/* --------------------------------------------------------------------- */
/* Initialize the backward recursion and the normalization sum. */
/* --------------------------------------------------------------------- */
L190:
	    ++n;
	    en += two;
	    tempb = zero;
	    tempa = one / p;
	    m = (n << 1) - (n / 2 << 2);
	    sum = zero;
	    i__1 = n / 2;
	    em = (doublereal) i__1;
	    alpem = em - one + *alpha;
	    alp2em = em + em + *alpha;
	    if (m != 0) {
		sum = tempa * alpem * alp2em / em;
	    }
	    nend = n - *nb;
	    if (nend > 0) {
/* --------------------------------------------------------------------- */
/* Recur backward via difference equation, calculating (but not */
/* storing) B(N), until N = NB. */
/* --------------------------------------------------------------------- */
		i__1 = nend;
		for (l = 1; l <= i__1; ++l) {
		    --n;
		    en -= two;
		    tempc = tempb;
		    tempb = tempa;
		    tempa = en * tempb / *x - tempc;
		    m = 2 - m;
		    if (m != 0) {
			em -= one;
			alp2em = em + em + *alpha;
			if (n == 1) {
			    goto L210;
			}
			alpem = em - one + *alpha;
			if (alpem == zero) {
			    alpem = one;
			}
			sum = (sum + tempa * alp2em) * alpem / em;
		    }
/* L200: */
		}
	    }
/* --------------------------------------------------------------------- */
/* Store B(NB). */
/* --------------------------------------------------------------------- */
L210:
	    b[n] = tempa;
	    if (nend >= 0) {
		if (*nb <= 1) {
		    alp2em = *alpha;
		    if (*alpha + one == one) {
			alp2em = one;
		    }
		    sum += b[1] * alp2em;
		    goto L250;
		} else {
/* --------------------------------------------------------------------- */
/* Calculate and store B(NB-1). */
/* --------------------------------------------------------------------- */
		    --n;
		    en -= two;
		    b[n] = en * tempa / *x - tempb;
		    if (n == 1) {
			goto L240;
		    }
		    m = 2 - m;
		    if (m != 0) {
			em -= one;
			alp2em = em + em + *alpha;
			alpem = em - one + *alpha;
			if (alpem == zero) {
			    alpem = one;
			}
			sum = (sum + b[n] * alp2em) * alpem / em;
		    }
		}
	    }
	    nend = n - 2;
	    if (nend != 0) {
/* --------------------------------------------------------------------- */
/* Calculate via difference equation and store B(N), until N = 2. */
/* --------------------------------------------------------------------- */
		i__1 = nend;
		for (l = 1; l <= i__1; ++l) {
		    --n;
		    en -= two;
		    b[n] = en * b[n + 1] / *x - b[n + 2];
		    m = 2 - m;
		    if (m != 0) {
			em -= one;
			alp2em = em + em + *alpha;
			alpem = em - one + *alpha;
			if (alpem == zero) {
			    alpem = one;
			}
			sum = (sum + b[n] * alp2em) * alpem / em;
		    }
/* L230: */
		}
	    }
/* --------------------------------------------------------------------- */
/* Calculate B(1). */
/* --------------------------------------------------------------------- */
	    b[1] = two * (*alpha + one) * b[2] / *x - b[3];
L240:
	    em -= one;
	    alp2em = em + em + *alpha;
	    if (alp2em == zero) {
		alp2em = one;
	    }
	    sum += b[1] * alp2em;
/* --------------------------------------------------------------------- */
/* Normalize.  Divide all B(N) by sum. */
/* --------------------------------------------------------------------- */
L250:
	    if (*alpha + one != one) {
		d__1 = *x * half;
		d__2 = -(*alpha);
		sum = sum * dgamma_(alpha) * pow_dd(&d__1, &d__2);
	    }
	    tempa = enmten;
	    if (sum > one) {
		tempa *= sum;
	    }
	    i__1 = *nb;
	    for (n = 1; n <= i__1; ++n) {
		if ((d__1 = b[n], abs(d__1)) < tempa) {
		    b[n] = zero;
		}
		b[n] /= sum;
/* L260: */
	    }
	}
/* --------------------------------------------------------------------- */
/* Error return -- X, NB, or ALPHA is out of range. */
/* --------------------------------------------------------------------- */
    } else {
	b[1] = zero;
	*ncalc = min(*nb,0) - 1;
    }
/* --------------------------------------------------------------------- */
/* Exit */
/* --------------------------------------------------------------------- */
L300:
    return 0;
/* ---------- Last line of RJBESL ---------- */
} /* rjbesl_ */

/// @endcond
/** @}*/