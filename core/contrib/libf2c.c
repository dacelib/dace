/* Extract from libf2c (11 March 2019) containing only those
   routines used by the contributed code.

   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/** \addtogroup DACE Contrib
 *  @{
 */

/// @cond

#include "f2c.h"

#undef abs
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif
double d_int(const doublereal *x)
{
return( (*x>0) ? floor(*x) : -floor(- *x) );
}

double pow_dd(const doublereal *ap, const doublereal *bp)
{
return(pow(*ap, *bp) );
}

double pow_di(const doublereal *ap, const integer *bp)
{
double pow, x;
integer n;
unsigned long u;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for(u = n; ; )
		{
		if(u & 01)
			pow *= x;
		if(u >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}
#ifdef __cplusplus
}
#endif

/// @endcond
/** @}*/
