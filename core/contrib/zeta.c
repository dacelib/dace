/*							zeta.c
 *
 *	Riemann zeta function of two arguments
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, q, y, zeta();
 *
 * y = zeta( x, q );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *
 *                 inf.
 *                  -        -x
 *   zeta(x,q)  =   >   (k+q)
 *                  -
 *                 k=0
 *
 * where x > 1 and q is not a negative integer or zero.
 * The Euler-Maclaurin summation formula is used to obtain
 * the expansion
 *
 *                n
 *                -       -x
 * zeta(x,q)  =   >  (k+q)
 *                -
 *               k=1
 *
 *           1-x                 inf.  B   x(x+1)...(x+2j)
 *      (n+q)           1         -     2j
 *  +  ---------  -  -------  +   >    --------------------
 *        x-1              x      -                   x+2j+1
 *                   2(n+q)      j=1       (2j)! (n+q)
 *
 * where the B2j are Bernoulli numbers.  Note that (see zetac.c)
 * zeta(x,1) = zetac(x) + 1.
 *
 *
 *
 * ACCURACY:
 *
 *
 *
 * REFERENCE:
 *
 * Gradshteyn, I. S., and I. M. Ryzhik, Tables of Integrals,
 * Series, and Products, p. 1073; Academic Press, 1980.
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier

Modified interface by Alexander Wittig, November 2019.
*/

#include <math.h>

/* Expansion coefficients
 * for Euler-Maclaurin summation formula
 * (2k)! / B2k
 * where B2k are Bernoulli numbers
 */
static double A[] = {
12.0,
-720.0,
30240.0,
-1209600.0,
47900160.0,
-1.8924375803183791606e9, /*1.307674368e12/691*/
7.47242496e10,
-2.950130727918164224e12, /*1.067062284288e16/3617*/
1.1646782814350067249e14, /*5.109094217170944e18/43867*/
-4.5979787224074726105e15, /*8.028576626982912e20/174611*/
1.8152105401943546773e17, /*1.5511210043330985984e23/854513*/
-7.1661652561756670113e18 /*1.6938241367317436694528e27/236364091*/
};
/* 30 Nov 86 -- error in third coefficient fixed */

static const double MACHEP = 1.11022302462515654042E-16;

double zeta_(const double x, const double q, unsigned int *err)
{
	double a, b, k, s, t, w;

	if( err ) *err = 0;

	if( x == 1.0 )
		goto retinf;

	if( x < 1.0 )
	{
domerr:
		if( err ) *err = 1;
		return(NAN);
	}

	if( q <= 0.0 )
	{
		if(q == floor(q))
		{
retinf:
			if( err ) *err = 2;
			return(INFINITY);
		}
		if( x != floor(x) )
			goto domerr; /* because q^-x not defined */
	}

	/* Euler-Maclaurin summation formula */
	/* Permit negative q but continue sum until n+q > +9 .
	* This case should be handled by a reflection formula.
	* If q<0 and x is an integer, there is a relation to
	* the polygamma function.
	*/
	s = pow( q, -x );
	a = q;
	b = 0.0;
	for( unsigned int i = 0; (i < 9) || (a <= 9.0); i++ )
	{
		a += 1.0;
		b = pow( a, -x );
		s += b;
		if( fabs(b/s) < MACHEP )
			goto done;
	}

	w = a;
	s += b*w/(x-1.0);
	s -= 0.5 * b;
	a = 1.0;
	k = 0.0;
	for( unsigned int i=0; i<12; i++ )
	{
		a *= x + k;
		b /= w;
		t = a*b/A[i];
		s = s + t;
		if( fabs(t/s) < MACHEP )
			goto done;
		k += 1.0;
		a *= x + k;
		b /= w;
		k += 1.0;
	}
done:
	return(s);
}
