/* rybesl.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int rybesl_(const doublereal *x, const doublereal *alpha, const integer *nb,
	doublereal *by, integer *ncalc)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal piby2 = 1.5707963267948966192;
    static doublereal pi = 3.1415926535897932385;
    static doublereal sq2bpi = .79788456080286535588;
    static doublereal pim5 = .70796326794896619231;
    static doublereal onbpi = .31830988618379067154;
    static doublereal del = 1e-8;
    static doublereal xmin = 4.46e-308;
    static doublereal xinf = 1.79e308;
    static doublereal eps = 1.11e-16;
    static doublereal thresh = 16.;
    static doublereal half = .5;
    static doublereal xlarge = 1e8;
    static doublereal ch[21] = { -6.7735241822398840964e-24,
	    -6.1455180116049879894e-23,2.9017595056104745456e-21,
	    1.3639417919073099464e-19,2.3826220476859635824e-18,
	    -9.0642907957550702534e-18,-1.4943667065169001769e-15,
	    -3.3919078305362211264e-14,-1.7023776642512729175e-13,
	    9.1609750938768647911e-12,2.4230957900482704055e-10,
	    1.7451364971382984243e-9,-3.3126119768180852711e-8,
	    -8.6592079961391259661e-7,-4.9717367041957398581e-6,
	    7.6309597585908126618e-5,.0012719271366545622927,
	    .0017063050710955562222,-.07685284084478667369,
	    -.28387654227602353814,.92187029365045265648 };
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;
    static doublereal eight = 8.;
    static doublereal one5 = 15.;
    static doublereal ten9 = 19.;
    static doublereal fivpi = 15.707963267948966192;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_int(const doublereal *), sqrt(doublereal), sin(doublereal), cos(
	    doublereal), log(doublereal), pow_dd(const doublereal *, const doublereal *);

    /* Local variables */
    doublereal b, c__, d__, e, f, g, h__;
    integer i__, k;
    doublereal p, q, r__, s, d1, d2, q0, x2;
    integer na;
    doublereal pa, qa, en, ya, ex, pa1, qa1, en1, ya1, den, odd, aye, div,
	    dmu, xna, enu, alfa, ddiv, even, term, gamma, cosmu, sinmu,
	    twobyx;

/* ---------------------------------------------------------------------- */

/*  This routine calculates Bessel functions Y SUB(N+ALPHA) (X) */
/*  for non-negative argument X, and non-negative order N+ALPHA. */


/* Explanation of variables in the calling sequence */

/* X     - Working precision non-negative real argument for which */
/*         Y's are to be calculated. */
/* ALPHA - Working precision fractional part of order for which */
/*         Y's are to be calculated.  0 .LE. ALPHA .LT. 1.0. */
/* NB    - Integer number of functions to be calculated, NB .GT. 0. */
/*         The first function calculated is of order ALPHA, and the */
/*         last is of order (NB - 1 + ALPHA). */
/* BY    - Working precision output vector of length NB.  If the */
/*         routine terminates normally (NCALC=NB), the vector BY */
/*         contains the functions Y(ALPHA,X), ... , Y(NB-1+ALPHA,X), */
/*         If (0 .LT. NCALC .LT. NB), BY(I) contains correct function */
/*         values for I .LE. NCALC, and contains the ratios */
/*         Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array. */
/* NCALC - Integer output variable indicating possible errors. */
/*         Before using the vector BY, the user should check that */
/*         NCALC=NB, i.e., all orders have been calculated to */
/*         the desired accuracy.  See error returns below. */


/* ******************************************************************* */
/* ******************************************************************* */

/* Explanation of machine-dependent constants */

/*   beta   = Radix for the floating-point system */
/*   p      = Number of significant base-beta digits in the */
/*            significand of a floating-point number */
/*   minexp = Smallest representable power of beta */
/*   maxexp = Smallest power of beta that overflows */
/*   EPS    = beta ** (-p) */
/*   DEL    = Machine number below which sin(x)/x = 1; approximately */
/*            SQRT(EPS). */
/*   XMIN   = Smallest acceptable argument for RBESY; approximately */
/*            max(2*beta**minexp,2/XINF), rounded up */
/*   XINF   = Largest positive machine number; approximately */
/*            beta**maxexp */
/*   THRESH = Lower bound for use of the asymptotic form; approximately */
/*            AINT(-LOG10(EPS/2.0))+1.0 */
/*   XLARGE = Upper bound on X; approximately 1/DEL, because the sine */
/*            and cosine functions have lost about half of their */
/*            precision at that point. */


/*     Approximate values for some important machines are: */

/*                        beta    p     minexp      maxexp      EPS */

/*  CRAY-1        (S.P.)    2    48     -8193        8191    3.55E-15 */
/*  Cyber 180/185 */
/*    under NOS   (S.P.)    2    48      -975        1070    3.55E-15 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (S.P.)    2    24      -126         128    5.96E-8 */
/*  IEEE (IBM/XT, */
/*    SUN, etc.)  (D.P.)    2    53     -1022        1024    1.11D-16 */
/*  IBM 3033      (D.P.)   16    14       -65          63    1.39D-17 */
/*  VAX           (S.P.)    2    24      -128         127    5.96E-8 */
/*  VAX D-Format  (D.P.)    2    56      -128         127    1.39D-17 */
/*  VAX G-Format  (D.P.)    2    53     -1024        1023    1.11D-16 */


/*                         DEL      XMIN      XINF     THRESH  XLARGE */

/* CRAY-1        (S.P.)  5.0E-8  3.67E-2466 5.45E+2465  15.0E0  2.0E7 */
/* Cyber 180/855 */
/*   under NOS   (S.P.)  5.0E-8  6.28E-294  1.26E+322   15.0E0  2.0E7 */
/* IEEE (IBM/XT, */
/*   SUN, etc.)  (S.P.)  1.0E-4  2.36E-38   3.40E+38     8.0E0  1.0E4 */
/* IEEE (IBM/XT, */
/*   SUN, etc.)  (D.P.)  1.0D-8  4.46D-308  1.79D+308   16.0D0  1.0D8 */
/* IBM 3033      (D.P.)  1.0D-8  2.77D-76   7.23D+75    17.0D0  1.0D8 */
/* VAX           (S.P.)  1.0E-4  1.18E-38   1.70E+38     8.0E0  1.0E4 */
/* VAX D-Format  (D.P.)  1.0D-9  1.18D-38   1.70D+38    17.0D0  1.0D9 */
/* VAX G-Format  (D.P.)  1.0D-8  2.23D-308  8.98D+307   16.0D0  1.0D8 */

/* ******************************************************************* */
/* ******************************************************************* */

/* Error returns */

/*  In case of an error, NCALC .NE. NB, and not all Y's are */
/*  calculated to the desired accuracy. */

/*  NCALC .LT. -1:  An argument is out of range. For example, */
/*       NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE. */
/*       XMAX.  In this case, BY(1) = 0.0, the remainder of the */
/*       BY-vector is not calculated, and NCALC is set to */
/*       MIN0(NB,0)-2  so that NCALC .NE. NB. */
/*  NCALC = -1:  Y(ALPHA,X) .GE. XINF.  The requested function */
/*       values are set to 0.0. */
/*  1 .LT. NCALC .LT. NB: Not all requested function values could */
/*       be calculated accurately.  BY(I) contains correct function */
/*       values for I .LE. NCALC, and and the remaining NB-NCALC */
/*       array elements contain 0.0. */


/* Intrinsic functions required are: */

/*     DBLE, EXP, INT, MAX, MIN, REAL, SQRT */


/* Acknowledgement */

/*  This program draws heavily on Temme's Algol program for Y(a,x) */
/*  and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's */
/*  scheme is used for  x < THRESH, and Campbell's scheme is used */
/*  in the asymptotic region.  Segments of code from both sources */
/*  have been translated into Fortran 77, merged, and heavily modified. */
/*  Modifications include parameterization of machine dependencies, */
/*  use of a new approximation for ln(gamma(x)), and built-in */
/*  protection against over/underflow. */

/* References: "Bessel functions J_nu(x) and Y_nu(x) of real */
/*              order and real argument," Campbell, J. B., */
/*              Comp. Phy. Comm. 18, 1979, pp. 133-142. */

/*             "On the numerical evaluation of the ordinary */
/*              Bessel function of the second kind," Temme, */
/*              N. M., J. Comput. Phys. 21, 1976, pp. 343-350. */

/*  Latest modification: March 19, 1990 */

/*  Modified by: W. J. Cody */
/*               Applied Mathematics Division */
/*               Argonne National Laboratory */
/*               Argonne, IL  60439 */

/* ---------------------------------------------------------------------- */
/* S    REAL */
/* ---------------------------------------------------------------------- */
/*  Mathematical constants */
/*    FIVPI = 5*PI */
/*    PIM5 = 5*PI - 15 */
/*    ONBPI = 1/PI */
/*    PIBY2 = PI/2 */
/*    SQ2BPI = SQUARE ROOT OF 2/PI */
/* ---------------------------------------------------------------------- */
/* S    DATA ZERO,HALF,ONE,TWO,THREE/0.0E0,0.5E0,1.0E0,2.0E0,3.0E0/ */
/* S    DATA EIGHT,ONE5,TEN9/8.0E0,15.0E0,1.9E1/ */
/* S    DATA FIVPI,PIBY2/1.5707963267948966192E1,1.5707963267948966192E0/ */
/* S    DATA PI,SQ2BPI/3.1415926535897932385E0,7.9788456080286535588E-1/ */
/* S    DATA PIM5,ONBPI/7.0796326794896619231E-1,3.1830988618379067154E-1/ */
    /* Parameter adjustments */
    --by;

    /* Function Body */
/* ---------------------------------------------------------------------- */
/*  Machine-dependent constants */
/* ---------------------------------------------------------------------- */
/* S    DATA DEL,XMIN,XINF,EPS/1.0E-4,2.36E-38,3.40E38,5.96E-8/ */
/* S    DATA THRESH,XLARGE/8.0E0,1.0E4/ */
/* ---------------------------------------------------------------------- */
/*  Coefficients for Chebyshev polynomial expansion of */
/*         1/gamma(1-x), abs(x) .le. .5 */
/* ---------------------------------------------------------------------- */
/* S    DATA CH/-0.67735241822398840964E-23,-0.61455180116049879894E-22, */
/* S   1         0.29017595056104745456E-20, 0.13639417919073099464E-18, */
/* S   2         0.23826220476859635824E-17,-0.90642907957550702534E-17, */
/* S   3        -0.14943667065169001769E-14,-0.33919078305362211264E-13, */
/* S   4        -0.17023776642512729175E-12, 0.91609750938768647911E-11, */
/* S   5         0.24230957900482704055E-09, 0.17451364971382984243E-08, */
/* S   6        -0.33126119768180852711E-07,-0.86592079961391259661E-06, */
/* S   7        -0.49717367041957398581E-05, 0.76309597585908126618E-04, */
/* S   8         0.12719271366545622927E-02, 0.17063050710955562222E-02, */
/* S   9        -0.76852840844786673690E-01,-0.28387654227602353814E+00, */
/* S   A         0.92187029365045265648E+00/ */
/* ---------------------------------------------------------------------- */
    ex = *x;
    enu = *alpha;
    if (*nb > 0 && *x >= xmin && ex < xlarge && enu >= zero && enu < one) {
	d__1 = enu + half;
	xna = d_int(&d__1);
	na = (integer) xna;
	if (na == 1) {
	    enu -= xna;
	}
	if (enu == -half) {
	    p = sq2bpi / sqrt(ex);
	    ya = p * sin(ex);
	    ya1 = -p * cos(ex);
	} else if (ex < three) {
/* ---------------------------------------------------------------------- */
/*  Use Temme's scheme for small X */
/* ---------------------------------------------------------------------- */
	    b = ex * half;
	    d__ = -log(b);
	    f = enu * d__;
	    d__1 = -enu;
	    e = pow_dd(&b, &d__1);
	    if (abs(enu) < del) {
		c__ = onbpi;
	    } else {
		c__ = enu / sin(enu * pi);
	    }
/* ---------------------------------------------------------------------- */
/*  Computation of sinh(f)/f */
/* ---------------------------------------------------------------------- */
	    if (abs(f) < one) {
		x2 = f * f;
		en = ten9;
		s = one;
		for (i__ = 1; i__ <= 9; ++i__) {
		    s = s * x2 / en / (en - one) + one;
		    en -= two;
/* L80: */
		}
	    } else {
		s = (e - one / e) * half / f;
	    }
/* ---------------------------------------------------------------------- */
/*  Computation of 1/gamma(1-a) using Chebyshev polynomials */
/* ---------------------------------------------------------------------- */
	    x2 = enu * enu * eight;
	    aye = ch[0];
	    even = zero;
	    alfa = ch[1];
	    odd = zero;
	    for (i__ = 3; i__ <= 19; i__ += 2) {
		even = -(aye + aye + even);
		aye = -even * x2 - aye + ch[i__ - 1];
		odd = -(alfa + alfa + odd);
		alfa = -odd * x2 - alfa + ch[i__];
/* L40: */
	    }
	    even = (even * half + aye) * x2 - aye + ch[20];
	    odd = (odd + alfa) * two;
	    gamma = odd * enu + even;
/* ---------------------------------------------------------------------- */
/*  End of computation of 1/gamma(1-a) */
/* ---------------------------------------------------------------------- */
	    g = e * gamma;
	    e = (e + one / e) * half;
	    f = two * c__ * (odd * e + even * s * d__);
	    e = enu * enu;
	    p = g * c__;
	    q = onbpi / g;
	    c__ = enu * piby2;
	    if (abs(c__) < del) {
		r__ = one;
	    } else {
		r__ = sin(c__) / c__;
	    }
	    r__ = pi * c__ * r__ * r__;
	    c__ = one;
	    d__ = -b * b;
	    h__ = zero;
	    ya = f + r__ * q;
	    ya1 = p;
	    en = zero;
L100:
	    en += one;
	    if ((d__1 = g / (one + abs(ya)), abs(d__1)) + (d__2 = h__ / (one
		    + abs(ya1)), abs(d__2)) > eps) {
		f = (f * en + p + q) / (en * en - e);
		c__ = c__ * d__ / en;
		p /= en - enu;
		q /= en + enu;
		g = c__ * (f + r__ * q);
		h__ = c__ * p - en * g;
		ya += g;
		ya1 += h__;
		goto L100;
	    }
	    ya = -ya;
	    ya1 = -ya1 / b;
	} else if (ex < thresh) {
/* ---------------------------------------------------------------------- */
/*  Use Temme's scheme for moderate X */
/* ---------------------------------------------------------------------- */
	    c__ = (half - enu) * (half + enu);
	    b = ex + ex;
	    e = ex * onbpi * cos(enu * pi) / eps;
	    e *= e;
	    p = one;
	    q = -ex;
	    r__ = one + ex * ex;
	    s = r__;
	    en = two;
L200:
	    if (r__ * en * en < e) {
		en1 = en + one;
		d__ = (en - one + c__ / en) / s;
		p = (en + en - p * d__) / en1;
		q = (-b + q * d__) / en1;
		s = p * p + q * q;
		r__ *= s;
		en = en1;
		goto L200;
	    }
	    f = p / s;
	    p = f;
	    g = -q / s;
	    q = g;
L220:
	    en -= one;
	    if (en > zero) {
		r__ = en1 * (two - p) - two;
		s = b + en1 * q;
		d__ = (en - one + c__ / en) / (r__ * r__ + s * s);
		p = d__ * r__;
		q = d__ * s;
		e = f + one;
		f = p * e - g * q;
		g = q * e + p * g;
		en1 = en;
		goto L220;
	    }
	    f = one + f;
	    d__ = f * f + g * g;
	    pa = f / d__;
	    qa = -g / d__;
	    d__ = enu + half - p;
	    q += ex;
	    pa1 = (pa * q - qa * d__) / ex;
	    qa1 = (qa * q + pa * d__) / ex;
	    b = ex - piby2 * (enu + half);
	    c__ = cos(b);
	    s = sin(b);
	    d__ = sq2bpi / sqrt(ex);
	    ya = d__ * (pa * s + qa * c__);
	    ya1 = d__ * (qa1 * s - pa1 * c__);
	} else {
/* ---------------------------------------------------------------------- */
/*  Use Campbell's asymptotic scheme. */
/* ---------------------------------------------------------------------- */
	    na = 0;
	    d__1 = ex / fivpi;
	    d1 = d_int(&d__1);
	    i__ = (integer) d1;
	    dmu = ex - one5 * d1 - d1 * pim5 - (*alpha + half) * piby2;
	    if (i__ - (i__ / 2 << 1) == 0) {
		cosmu = cos(dmu);
		sinmu = sin(dmu);
	    } else {
		cosmu = -cos(dmu);
		sinmu = -sin(dmu);
	    }
	    ddiv = eight * ex;
	    dmu = *alpha;
	    den = sqrt(ex);
	    for (k = 1; k <= 2; ++k) {
		p = cosmu;
		cosmu = sinmu;
		sinmu = -p;
		d1 = (two * dmu - one) * (two * dmu + one);
		d2 = zero;
		div = ddiv;
		p = zero;
		q = zero;
		q0 = d1 / div;
		term = q0;
		for (i__ = 2; i__ <= 20; ++i__) {
		    d2 += eight;
		    d1 -= d2;
		    div += ddiv;
		    term = -term * d1 / div;
		    p += term;
		    d2 += eight;
		    d1 -= d2;
		    div += ddiv;
		    term = term * d1 / div;
		    q += term;
		    if (abs(term) <= eps) {
			goto L320;
		    }
/* L310: */
		}
L320:
		p += one;
		q += q0;
		if (k == 1) {
		    ya = sq2bpi * (p * cosmu - q * sinmu) / den;
		} else {
		    ya1 = sq2bpi * (p * cosmu - q * sinmu) / den;
		}
		dmu += one;
/* L350: */
	    }
	}
	if (na == 1) {
	    h__ = two * (enu + one) / ex;
	    if (h__ > one) {
		if (abs(ya1) > xinf / h__) {
		    h__ = zero;
		    ya = zero;
		}
	    }
	    h__ = h__ * ya1 - ya;
	    ya = ya1;
	    ya1 = h__;
	}
/* ---------------------------------------------------------------------- */
/*  Now have first one or two Y's */
/* ---------------------------------------------------------------------- */
	by[1] = ya;
	by[2] = ya1;
	if (ya1 == zero) {
	    *ncalc = 1;
	} else {
	    aye = one + *alpha;
	    twobyx = two / ex;
	    *ncalc = 2;
	    i__1 = *nb;
	    for (i__ = 3; i__ <= i__1; ++i__) {
		if (twobyx < one) {
		    if ((d__1 = by[i__ - 1], abs(d__1)) * twobyx >= xinf /
			    aye) {
			goto L450;
		    }
		} else {
		    if ((d__1 = by[i__ - 1], abs(d__1)) >= xinf / aye /
			    twobyx) {
			goto L450;
		    }
		}
		by[i__] = twobyx * aye * by[i__ - 1] - by[i__ - 2];
		aye += one;
		++(*ncalc);
/* L400: */
	    }
	}
L450:
	i__1 = *nb;
	for (i__ = *ncalc + 1; i__ <= i__1; ++i__) {
	    by[i__] = zero;
/* L460: */
	}
    } else {
	by[1] = zero;
	*ncalc = min(*nb,0) - 1;
    }
/* L900: */
    return 0;
/* ---------- Last line of RYBESL ---------- */
} /* rybesl_ */

/// @endcond
/** @}*/