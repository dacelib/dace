/* rkbesl.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int rkbesl_(const doublereal *x, const doublereal *alpha, const integer *nb,
	const integer *ize, doublereal *bk, integer *ncalc)
{
    /* Initialized data */

    static doublereal half = .5;
    static doublereal sqxmin = 1.49e-154;
    static doublereal xinf = 1.79e308;
    static doublereal xmin = 2.23e-308;
    static doublereal xmax = 705.342;
    static doublereal p[8] = { .805629875690432845,20.4045500205365151,
	    157.705605106676174,536.671116469207504,900.382759291288778,
	    730.923886650660393,229.299301509425145,.822467033424113231 };
    static doublereal q[7] = { 29.4601986247850434,277.577868510221208,
	    1206.70325591027438,2762.91444159791519,3443.74050506564618,
	    2210.63190113378647,572.267338359892221 };
    static doublereal r__[5] = { -.48672575865218401848,13.079485869097804016,
	    -101.96490580880537526,347.65409106507813131,
	    3.495898124521934782e-4 };
    static doublereal s[4] = { -25.579105509976461286,212.57260432226544008,
	    -610.69018684944109624,422.69668805777760407 };
    static doublereal t[6] = { 1.6125990452916363814e-10,
	    2.5051878502858255354e-8,2.7557319615147964774e-6,
	    1.9841269840928373686e-4,.0083333333333334751799,
	    .16666666666666666446 };
    static doublereal estm[6] = { 52.0583,5.7607,2.7782,14.4303,185.3004,
	    9.3715 };
    static doublereal one = 1.;
    static doublereal estf[7] = { 41.8341,7.1075,6.4306,42.511,1.35633,
	    84.5096,20. };
    static doublereal two = 2.;
    static doublereal zero = 0.;
    static doublereal four = 4.;
    static doublereal tinyx = 1e-10;
    static doublereal a = .11593151565841244881;
    static doublereal d__ = .797884560802865364;
    static doublereal eps = 2.22e-16;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal), exp(doublereal), sinh(doublereal), sqrt(
	    doublereal), d_int(const doublereal *);

    /* Local variables */
    doublereal c__;
    integer i__, j, k, m;
    doublereal d1, d2, d3, f0, f1, f2, p0, q0, t1, t2, dm, ex, bk1, bk2, enu;
    integer iend;
    doublereal x2by4, twox, blpha;
    integer itemp;
    doublereal ratio, wminf, twonu;
    integer mplus1;

/* ------------------------------------------------------------------- */

/*  This FORTRAN 77 routine calculates modified Bessel functions */
/*  of the second kind, K SUB(N+ALPHA) (X), for non-negative */
/*  argument X, and non-negative order N+ALPHA, with or without */
/*  exponential scaling. */

/*  Explanation of variables in the calling sequence */

/*  Description of output values .. */

/* X     - Working precision non-negative real argument for which */
/*         K's or exponentially scaled K's (K*EXP(X)) */
/*         are to be calculated.  If K's are to be calculated, */
/*         X must not be greater than XMAX (see below). */
/* ALPHA - Working precision fractional part of order for which */
/*         K's or exponentially scaled K's (K*EXP(X)) are */
/*         to be calculated.  0 .LE. ALPHA .LT. 1.0. */
/* NB    - Integer number of functions to be calculated, NB .GT. 0. */
/*         The first function calculated is of order ALPHA, and the */
/*         last is of order (NB - 1 + ALPHA). */
/* IZE   - Integer type.  IZE = 1 if unscaled K's are to be calculated, */
/*         and 2 if exponentially scaled K's are to be calculated. */
/* BK    - Working precision output vector of length NB.  If the */
/*         routine terminates normally (NCALC=NB), the vector BK */
/*         contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X), */
/*         or the corresponding exponentially scaled functions. */
/*         If (0 .LT. NCALC .LT. NB), BK(I) contains correct function */
/*         values for I .LE. NCALC, and contains the ratios */
/*         K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array. */
/* NCALC - Integer output variable indicating possible errors. */
/*         Before using the vector BK, the user should check that */
/*         NCALC=NB, i.e., all orders have been calculated to */
/*         the desired accuracy.  See error returns below. */


/* ******************************************************************* */
/* ******************************************************************* */

/* Explanation of machine-dependent constants */

/*   beta   = Radix for the floating-point system */
/*   minexp = Smallest representable power of beta */
/*   maxexp = Smallest power of beta that overflows */
/*   EPS    = The smallest positive floating-point number such that */
/*            1.0+EPS .GT. 1.0 */
/*   XMAX   = Upper limit on the magnitude of X when IZE=1;  Solution */
/*            to equation: */
/*               W(X) * (1-1/8X+9/128X**2) = beta**minexp */
/*            where  W(X) = EXP(-X)*SQRT(PI/2X) */
/*   SQXMIN = Square root of beta**minexp */
/*   XINF   = Largest positive machine number; approximately */
/*            beta**maxexp */
/*   XMIN   = Smallest positive machine number; approximately */
/*            beta**minexp */


/*     Approximate values for some important machines are: */

/*                          beta       minexp      maxexp      EPS */

/*  CRAY-1        (S.P.)      2        -8193        8191    7.11E-15 */
/*  Cyber 180/185 */
/*    under NOS   (S.P.)      2         -975        1070    3.55E-15 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (S.P.)      2         -126         128    1.19E-7 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (D.P.)      2        -1022        1024    2.22D-16 */
/*  IBM 3033      (D.P.)     16          -65          63    2.22D-16 */
/*  VAX           (S.P.)      2         -128         127    5.96E-8 */
/*  VAX D-Format  (D.P.)      2         -128         127    1.39D-17 */
/*  VAX G-Format  (D.P.)      2        -1024        1023    1.11D-16 */


/*                         SQXMIN       XINF        XMIN      XMAX */

/* CRAY-1        (S.P.)  6.77E-1234  5.45E+2465  4.59E-2467 5674.858 */
/* Cyber 180/855 */
/*   under NOS   (S.P.)  1.77E-147   1.26E+322   3.14E-294   672.788 */
/* IEEE (IBM/XT, */
/*   SUN, etc.)  (S.P.)  1.08E-19    3.40E+38    1.18E-38     85.337 */
/* IEEE (IBM/XT, */
/*   SUN, etc.)  (D.P.)  1.49D-154   1.79D+308   2.23D-308   705.342 */
/* IBM 3033      (D.P.)  7.35D-40    7.23D+75    5.40D-79    177.852 */
/* VAX           (S.P.)  5.42E-20    1.70E+38    2.94E-39     86.715 */
/* VAX D-Format  (D.P.)  5.42D-20    1.70D+38    2.94D-39     86.715 */
/* VAX G-Format  (D.P.)  7.46D-155   8.98D+307   5.57D-309   706.728 */

/* ******************************************************************* */
/* ******************************************************************* */

/* Error returns */

/*  In case of an error, NCALC .NE. NB, and not all K's are */
/*  calculated to the desired accuracy. */

/*  NCALC .LT. -1:  An argument is out of range. For example, */
/*       NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE. */
/*       XMAX.  In this case, the B-vector is not calculated, */
/*       and NCALC is set to MIN0(NB,0)-2  so that NCALC .NE. NB. */
/*  NCALC = -1:  Either  K(ALPHA,X) .GE. XINF  or */
/*       K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) .GE. XINF.  In this case, */
/*       the B-vector is not calculated.  Note that again */
/*       NCALC .NE. NB. */

/*  0 .LT. NCALC .LT. NB: Not all requested function values could */
/*       be calculated accurately.  BK(I) contains correct function */
/*       values for I .LE. NCALC, and contains the ratios */
/*       K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array. */


/* Intrinsic functions required are: */

/*     ABS, AINT, EXP, INT, LOG, MAX, MIN, SINH, SQRT */


/* Acknowledgement */

/*  This program is based on a program written by J. B. Campbell */
/*  (2) that computes values of the Bessel functions K of real */
/*  argument and real order.  Modifications include the addition */
/*  of non-scaled functions, parameterization of machine */
/*  dependencies, and the use of more accurate approximations */
/*  for SINH and SIN. */

/* References: "On Temme's Algorithm for the Modified Bessel */
/*              Functions of the Third Kind," Campbell, J. B., */
/*              TOMS 6(4), Dec. 1980, pp. 581-586. */

/*             "A FORTRAN IV Subroutine for the Modified Bessel */
/*              Functions of the Third Kind of Real Order and Real */
/*              Argument," Campbell, J. B., Report NRC/ERB-925, */
/*              National Research Council, Canada. */

/*  Latest modification: May 30, 1989 */

/*  Modified by: W. J. Cody and L. Stoltz */
/*               Applied Mathematics Division */
/*               Argonne National Laboratory */
/*               Argonne, IL  60439 */

/* ------------------------------------------------------------------- */
/* S    REAL */
/* --------------------------------------------------------------------- */
/*  Mathematical constants */
/*    A = LOG(2.D0) - Euler's constant */
/*    D = SQRT(2.D0/PI) */
/* --------------------------------------------------------------------- */
/* S    DATA HALF,ONE,TWO,ZERO/0.5E0,1.0E0,2.0E0,0.0E0/ */
/* S    DATA FOUR,TINYX/4.0E0,1.0E-10/ */
/* S    DATA A/ 0.11593151565841244881E0/,D/0.797884560802865364E0/ */
    /* Parameter adjustments */
    --bk;

    /* Function Body */
/* --------------------------------------------------------------------- */
/*  Machine dependent parameters */
/* --------------------------------------------------------------------- */
/* S    DATA EPS/1.19E-7/,SQXMIN/1.08E-19/,XINF/3.40E+38/ */
/* S    DATA XMIN/1.18E-38/,XMAX/85.337E0/ */
/* --------------------------------------------------------------------- */
/*  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA */
/*                                         + Euler's constant */
/*         Coefficients converted from hex to decimal and modified */
/*         by W. J. Cody, 2/26/82 */
/*  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA) */
/*  T    - Approximation for SINH(Y)/Y */
/* --------------------------------------------------------------------- */
/* S    DATA P/ 0.805629875690432845E00,    0.204045500205365151E02, */
/* S   1        0.157705605106676174E03,    0.536671116469207504E03, */
/* S   2        0.900382759291288778E03,    0.730923886650660393E03, */
/* S   3        0.229299301509425145E03,    0.822467033424113231E00/ */
/* S    DATA Q/ 0.294601986247850434E02,    0.277577868510221208E03, */
/* S   1        0.120670325591027438E04,    0.276291444159791519E04, */
/* S   2        0.344374050506564618E04,    0.221063190113378647E04, */
/* S   3        0.572267338359892221E03/ */
/* S    DATA R/-0.48672575865218401848E+0,  0.13079485869097804016E+2, */
/* S   1       -0.10196490580880537526E+3,  0.34765409106507813131E+3, */
/* S   2        0.34958981245219347820E-3/ */
/* S    DATA S/-0.25579105509976461286E+2,  0.21257260432226544008E+3, */
/* S   1       -0.61069018684944109624E+3,  0.42269668805777760407E+3/ */
/* S    DATA T/ 0.16125990452916363814E-9, 0.25051878502858255354E-7, */
/* S   1        0.27557319615147964774E-5, 0.19841269840928373686E-3, */
/* S   2        0.83333333333334751799E-2, 0.16666666666666666446E+0/ */
/* S    DATA ESTM/5.20583E1, 5.7607E0, 2.7782E0, 1.44303E1, 1.853004E2, */
/* S   1          9.3715E0/ */
/* S    DATA ESTF/4.18341E1, 7.1075E0, 6.4306E0, 4.25110E1, 1.35633E0, */
/* S   1          8.45096E1, 2.0E1/ */
/* --------------------------------------------------------------------- */
    ex = *x;
    enu = *alpha;
    *ncalc = min(*nb,0) - 2;
    if (*nb > 0 && (enu >= zero && enu < one) && (*ize >= 1 && *ize <= 2) && (
	    *ize != 1 || ex <= xmax) && ex > zero) {
	k = 0;
	if (enu < sqxmin) {
	    enu = zero;
	}
	if (enu > half) {
	    k = 1;
	    enu -= one;
	}
	twonu = enu + enu;
	iend = *nb + k - 1;
	c__ = enu * enu;
	d3 = -c__;
	if (ex <= one) {
/* --------------------------------------------------------------------- */
/*  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA */
/*                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA */
/* --------------------------------------------------------------------- */
	    d1 = zero;
	    d2 = p[0];
	    t1 = one;
	    t2 = q[0];
	    for (i__ = 2; i__ <= 7; i__ += 2) {
		d1 = c__ * d1 + p[i__ - 1];
		d2 = c__ * d2 + p[i__];
		t1 = c__ * t1 + q[i__ - 1];
		t2 = c__ * t2 + q[i__];
/* L10: */
	    }
	    d1 = enu * d1;
	    t1 = enu * t1;
	    f1 = log(ex);
	    f0 = a + enu * (p[7] - enu * (d1 + d2) / (t1 + t2)) - f1;
	    q0 = exp(-enu * (a - enu * (p[7] + enu * (d1 - d2) / (t1 - t2)) -
		    f1));
	    f1 = enu * f0;
	    p0 = exp(f1);
/* --------------------------------------------------------------------- */
/*  Calculation of F0 = */
/* --------------------------------------------------------------------- */
	    d1 = r__[4];
	    t1 = one;
	    for (i__ = 1; i__ <= 4; ++i__) {
		d1 = c__ * d1 + r__[i__ - 1];
		t1 = c__ * t1 + s[i__ - 1];
/* L20: */
	    }
	    if (abs(f1) <= half) {
		f1 *= f1;
		d2 = zero;
		for (i__ = 1; i__ <= 6; ++i__) {
		    d2 = f1 * d2 + t[i__ - 1];
/* L30: */
		}
		d2 = f0 + f0 * f1 * d2;
	    } else {
		d2 = sinh(f1) / enu;
	    }
	    f0 = d2 - enu * d1 / (t1 * p0);
	    if (ex <= tinyx) {
/* -------------------------------------------------------------------- */
/*  X.LE.1.0E-10 */
/*  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X) */
/* -------------------------------------------------------------------- */
		bk[1] = f0 + ex * f0;
		if (*ize == 1) {
		    bk[1] -= ex * bk[1];
		}
		ratio = p0 / f0;
		c__ = ex * xinf;
		if (k != 0) {
/* -------------------------------------------------------------------- */
/*  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X), */
/*  ALPHA .GE. 1/2 */
/* -------------------------------------------------------------------- */
		    *ncalc = -1;
		    if (bk[1] >= c__ / ratio) {
			goto L500;
		    }
		    bk[1] = ratio * bk[1] / ex;
		    twonu += two;
		    ratio = twonu;
		}
		*ncalc = 1;
		if (*nb == 1) {
		    goto L500;
		}
/* -------------------------------------------------------------------- */
/*  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1 */
/* -------------------------------------------------------------------- */
		*ncalc = -1;
		i__1 = *nb;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    if (ratio >= c__) {
			goto L500;
		    }
		    bk[i__] = ratio / ex;
		    twonu += two;
		    ratio = twonu;
/* L80: */
		}
		*ncalc = 1;
		goto L420;
	    } else {
/* -------------------------------------------------------------------- */
/*  1.0E-10 .LT. X .LE. 1.0 */
/* -------------------------------------------------------------------- */
		c__ = one;
		x2by4 = ex * ex / four;
		p0 = half * p0;
		q0 = half * q0;
		d1 = -one;
		d2 = zero;
		bk1 = zero;
		bk2 = zero;
		f1 = f0;
		f2 = p0;
L100:
		d1 += two;
		d2 += one;
		d3 = d1 + d3;
		c__ = x2by4 * c__ / d2;
		f0 = (d2 * f0 + p0 + q0) / d3;
		p0 /= d2 - enu;
		q0 /= d2 + enu;
		t1 = c__ * f0;
		t2 = c__ * (p0 - d2 * f0);
		bk1 += t1;
		bk2 += t2;
		if ((d__1 = t1 / (f1 + bk1), abs(d__1)) > eps || (d__2 = t2 /
			(f2 + bk2), abs(d__2)) > eps) {
		    goto L100;
		}
		bk1 = f1 + bk1;
		bk2 = two * (f2 + bk2) / ex;
		if (*ize == 2) {
		    d1 = exp(ex);
		    bk1 *= d1;
		    bk2 *= d1;
		}
		wminf = estf[0] * ex + estf[1];
	    }
	} else if (eps * ex > one) {
/* -------------------------------------------------------------------- */
/*  X .GT. ONE/EPS */
/* -------------------------------------------------------------------- */
	    *ncalc = *nb;
	    bk1 = one / (d__ * sqrt(ex));
	    i__1 = *nb;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		bk[i__] = bk1;
/* L110: */
	    }
	    goto L500;
	} else {
/* -------------------------------------------------------------------- */
/*  X .GT. 1.0 */
/* -------------------------------------------------------------------- */
	    twox = ex + ex;
	    blpha = zero;
	    ratio = zero;
	    if (ex <= four) {
/* -------------------------------------------------------------------- */
/*  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 .LE. X .LE. 4.0 */
/* -------------------------------------------------------------------- */
		d__1 = estm[0] / ex + estm[1];
		d2 = d_int(&d__1);
		m = (integer) d2;
		d1 = d2 + d2;
		d2 -= half;
		d2 *= d2;
		i__1 = m;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    d1 -= two;
		    d2 -= d1;
		    ratio = (d3 + d2) / (twox + d1 - ratio);
/* L120: */
		}
/* -------------------------------------------------------------------- */
/*  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward */
/*    recurrence and K(ALPHA,X) from the wronskian */
/* -------------------------------------------------------------------- */
		d__1 = estm[2] * ex + estm[3];
		d2 = d_int(&d__1);
		m = (integer) d2;
		c__ = abs(enu);
		d3 = c__ + c__;
		d1 = d3 - one;
		f1 = xmin;
		f0 = (two * (c__ + d2) / ex + half * ex / (c__ + d2 + one)) *
			xmin;
		i__1 = m;
		for (i__ = 3; i__ <= i__1; ++i__) {
		    d2 -= one;
		    f2 = (d3 + d2 + d2) * f0;
		    blpha = (one + d1 / d2) * (f2 + blpha);
		    f2 = f2 / ex + f1;
		    f1 = f0;
		    f0 = f2;
/* L130: */
		}
		f1 = (d3 + two) * f0 / ex + f1;
		d1 = zero;
		t1 = one;
		for (i__ = 1; i__ <= 7; ++i__) {
		    d1 = c__ * d1 + p[i__ - 1];
		    t1 = c__ * t1 + q[i__ - 1];
/* L140: */
		}
		p0 = exp(c__ * (a + c__ * (p[7] - c__ * d1 / t1) - log(ex))) /
			 ex;
		f2 = (c__ + half - ratio) * f1 / ex;
		bk1 = p0 + (d3 * f0 - f2 + f0 + blpha) / (f2 + f1 + f0) * p0;
		if (*ize == 1) {
		    bk1 *= exp(-ex);
		}
		wminf = estf[2] * ex + estf[3];
	    } else {
/* -------------------------------------------------------------------- */
/*  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward */
/*  recurrence, for  X .GT. 4.0 */
/* -------------------------------------------------------------------- */
		d__1 = estm[4] / ex + estm[5];
		dm = d_int(&d__1);
		m = (integer) dm;
		d2 = dm - half;
		d2 *= d2;
		d1 = dm + dm;
		i__1 = m;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    dm -= one;
		    d1 -= two;
		    d2 -= d1;
		    ratio = (d3 + d2) / (twox + d1 - ratio);
		    blpha = (ratio + ratio * blpha) / dm;
/* L160: */
		}
		bk1 = one / ((d__ + d__ * blpha) * sqrt(ex));
		if (*ize == 1) {
		    bk1 *= exp(-ex);
		}
		wminf = estf[4] * (ex - (d__1 = ex - estf[6], abs(d__1))) +
			estf[5];
	    }
/* -------------------------------------------------------------------- */
/*  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and */
/*    K(ALPHA+1,X)/K(ALPHA,X) */
/* -------------------------------------------------------------------- */
	    bk2 = bk1 + bk1 * (enu + half - ratio) / ex;
	}
/* -------------------------------------------------------------------- */
/*  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1, */
/*  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1 */
/* -------------------------------------------------------------------- */
	*ncalc = *nb;
	bk[1] = bk1;
	if (iend == 0) {
	    goto L500;
	}
	j = 2 - k;
	if (j > 0) {
	    bk[j] = bk2;
	}
	if (iend == 1) {
	    goto L500;
	}
/* Computing MIN */
	i__1 = (integer) (wminf - enu);
	m = min(i__1,iend);
	i__1 = m;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    t1 = bk1;
	    bk1 = bk2;
	    twonu += two;
	    if (ex < one) {
		if (bk1 >= xinf / twonu * ex) {
		    goto L195;
		}
		goto L187;
	    } else {
		if (bk1 / ex >= xinf / twonu) {
		    goto L195;
		}
	    }
L187:
	    bk2 = twonu / ex * bk1 + t1;
	    itemp = i__;
	    ++j;
	    if (j > 0) {
		bk[j] = bk2;
	    }
/* L190: */
	}
L195:
	m = itemp;
	if (m == iend) {
	    goto L500;
	}
	ratio = bk2 / bk1;
	mplus1 = m + 1;
	*ncalc = -1;
	i__1 = iend;
	for (i__ = mplus1; i__ <= i__1; ++i__) {
	    twonu += two;
	    ratio = twonu / ex + one / ratio;
	    ++j;
	    if (j > 1) {
		bk[j] = ratio;
	    } else {
		if (bk2 >= xinf / ratio) {
		    goto L500;
		}
		bk2 = ratio * bk2;
	    }
/* L410: */
	}
/* Computing MAX */
	i__1 = mplus1 - k;
	*ncalc = max(i__1,1);
	if (*ncalc == 1) {
	    bk[1] = bk2;
	}
	if (*nb == 1) {
	    goto L500;
	}
L420:
	j = *ncalc + 1;
	i__1 = *nb;
	for (i__ = j; i__ <= i__1; ++i__) {
	    if (bk[*ncalc] >= xinf / bk[i__]) {
		goto L500;
	    }
	    bk[i__] = bk[*ncalc] * bk[i__];
	    *ncalc = i__;
/* L430: */
	}
    }
L500:
    return 0;
/* ---------- Last line of RKBESL ---------- */
} /* rkbesl_ */

/// @endcond
/** @}*/