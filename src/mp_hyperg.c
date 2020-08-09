/*							mp_hyperg.c
 *
 *	Confluent hypergeometric function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, hyperg();
 *
 * y = hyperg( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the confluent hypergeometric function
 *
 *                          1           2
 *                       a x    a(a+1) x
 *   F ( a,b;x )  =  1 + ---- + --------- + ...
 *  1 1                  b 1!   b(b+1) 2!
 *
 * Many higher transcendental functions are special cases of
 * this power series.
 *
 * As is evident from the formula, b must not be a negative
 * integer or zero unless a is an integer with 0 >= a > b.
 *
 * The routine attempts both a direct summation of the series
 * and an asymptotic expansion.  In each case error due to
 * roundoff, cancellation, and nonconvergence is estimated.
 * The result with smaller estimated error is returned.
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random points (a, b, x), all three variables
 * ranging from 0 to 30.
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,30         2000       1.2e-15     1.3e-16
 qtst1:
 21800   max =  1.4200E-14   rms =  1.0841E-15  ave = -5.3640E-17 
 ltstd:
 25500   max = 1.2759e-14   rms = 3.7155e-16  ave = 1.5384e-18 
 *    IEEE      0,30        30000       1.8e-14     1.1e-15
 *
 * Larger errors can be observed when b is near a negative
 * integer or zero.  Certain combinations of arguments yield
 * serious cancellation error in the power series summation
 * and also are not in the region of near convergence of the
 * asymptotic series.  An error message is printed if the
 * self-estimated relative error is greater than 1.0e-12.
 *
 */

/*							hyperg.c */


/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier
*/

#include "mconf.h"
#include <Rmath.h>
#include <mpfr.h>

#ifdef ANSIPROT
extern double exp ( double );
extern double log ( double );
extern double fabs ( double );
double gam(double), lgam(double);
double hyp2f0 ( double, double, double, int, double * );
static double mphy1f1p(double, double, double, double *);
static double mphy1f1a(double, double, double, double *);
double hyperg (double, double, double);
#else
double exp(), log(), gammafn(), lgammafn(),fabs(), 
double hyp2f0();
double gam(), lgam();
static double mphy1f1p();
static double mphy1f1a();
double hyperg();
#endif
extern double MAXNUM, MACHEP;

double mphyperg( a, b, x)
double a, b, x;
{
double asum, psum, acanc, pcanc, temp;
    
/* See if a Kummer transformation will help */
temp = b - a;
if( fabs(temp) < 0.001 * fabs(a) )
	return( exp(x) * hyperg( temp, b, -x )  );


psum = mphy1f1p( a, b, x, &pcanc );
if( pcanc < 1.0e-15 )
	goto done;


/* try asymptotic series */

asum = mphy1f1a( a, b, x, &acanc );


/* Pick the result with less estimated error */

if( acanc < pcanc )
	{
	pcanc = acanc;
	psum = asum;
	}
// For gglm R package, I change the precision here:
// edited by Yingbo Li, 01/04/2013  
done:
//if( pcanc > 1.0e-12 )
  if( pcanc > 1.0e-12 )
	mtherr( "hyperg", PLOSS );

return( psum );
}




/* Power series summation for confluent hypergeometric function		*/

static double mphy1f1p( a, b, x, err )
  double a, b, x;
double *err;
{
  
  
  double n, t, temp, an, bn, pcanc, u_d;
  mpfr_t u, a0, sum, maxt, res; 
  
  an = a;
  bn = b;
  n = 1.0;
  t = 1.0;
  
  mpfr_init2 (u, 200);
  mpfr_init_set_d (a0, 1.0, 200);
  mpfr_init_set_d (sum, 1.0, 200);
  mpfr_init_set_d (maxt, 0.0, 200);
  mpfr_init2 (res, 200);
  
  while( t > MACHEP )
  {
    if( bn == 0 )			
    {
      mtherr( "hyperg", SING );
      return( MAXNUM );	
    }
    if( an == 0 )			
      return( mpfr_get_d(sum, MPFR_RNDN) );
    if( n > 200 ) 
      goto pdone;
    
    u_d = x * (an / (bn * n));
    mpfr_set_d (u, u_d, MPFR_RNDN);
      
    /* check for blowup */
    temp = fabs(mpfr_get_d(u, MPFR_RNDN));
    if ( (temp > 1.0) && (mpfr_get_d(maxt, MPFR_RNDN) > (MAXNUM/temp)) )
    {
      pcanc = 1.0;  /* estimate 100% error */
      goto blowup;
    } 
    
    mpfr_mul (a0, a0, u, MPFR_RNDN); 
    mpfr_add (sum, sum, a0, MPFR_RNDN);
    
    t = fabs(mpfr_get_d(a0, MPFR_RNDN));
    if (t > mpfr_get_d(maxt, MPFR_RNDN)) {
      mpfr_set_d(maxt, t, MPFR_RNDN);
    }
    
    an += 1.0;
    bn += 1.0;
    n += 1.0;
  }
  
  pdone:
    
    /* estimate error due to roundoff and cancellation */
    if( mpfr_get_d(sum, MPFR_RNDN) != 0.0 )
      mpfr_div_d (maxt, maxt, fabs(mpfr_get_d(sum, MPFR_RNDN)), MPFR_RNDN);
    
    mpfr_mul_d (maxt, maxt, MACHEP, MPFR_RNDN);   /* this way avoids multiply overflow */
    pcanc = fabs( MACHEP * n  +  mpfr_get_d(maxt, MPFR_RNDN) );
    
    blowup:
      
      *err = pcanc;
    
    mpfr_set (res, sum, MPFR_RNDN);
    
    mpfr_clear (u);
    mpfr_clear (a0);
    mpfr_clear (sum);
    mpfr_clear (maxt);
    
    return( mpfr_get_d(res, MPFR_RNDN) ); 
    
}


/*							hy1f1a()	*/
/* asymptotic formula for hypergeometric function:
 *
 *        (    -a                         
 *  --    ( |z|                           
 * |  (b) ( -------- 2f0( a, 1+a-b, -1/x )
 *        (  --                           
 *        ( |  (b-a)                      
 *
 *
 *                                x    a-b                     )
 *                               e  |x|                        )
 *                             + -------- 2f0( b-a, 1-a, 1/x ) )
 *                                --                           )
 *                               |  (a)                        )
 */

static double mphy1f1a( a, b, x, err )
double a, b, x;
double *err;
{
double h1, h2, t, u, temp, acanc, asum, err1, err2;

if( x == 0 )
	{
	acanc = 1.0;
	asum = MAXNUM;
	goto adone;
	}
temp = log( fabs(x) );
t = x + temp * (a-b);
u = -temp * a;

if( b > 0 )
	{
	temp = lgam(b);
	t += temp;
	u += temp;
	}

h1 = hyp2f0( a, a-b+1, -1.0/x, 1, &err1 );

temp = exp(u) / gam(b-a);
h1 *= temp;
err1 *= temp;

h2 = hyp2f0( b-a, 1.0-a, 1.0/x, 2, &err2 );

if( a < 0 )
	temp = exp(t) / gam(a);
else
	temp = exp( t - lgam(a) );

h2 *= temp;
err2 *= temp;

if( x < 0.0 )
	asum = h1;
else
	asum = h2;

acanc = fabs(err1) + fabs(err2);


if( b < 0 )
	{
	temp = gam(b);
	asum *= temp;
	acanc *= fabs(temp);
	}


if( asum != 0.0 )
	acanc /= fabs(asum);

acanc *= 30.0;	/* fudge factor, since error of asymptotic formula
		 * often seems this much larger than advertised */

adone:


*err = acanc;
return( asum );
}

/*							hyp2f0()	*/

double mphyp2f0( a, b, x, type, err )
double a, b, x;
int type;	/* determines what converging factor to use */
double *err;
{
double a0, alast, t, tlast, maxt;
double n, an, bn, u, sum, temp;

an = a;
bn = b;
a0 = 1.0e0;
alast = 1.0e0;
sum = 0.0;
n = 1.0e0;
t = 1.0e0;
tlast = 1.0e9;
maxt = 0.0;

do
	{
	if( an == 0 )
		goto pdone;
	if( bn == 0 )
		goto pdone;

	u = an * (bn * x / n);

	/* check for blowup */
	temp = fabs(u);
	if( (temp > 1.0 ) && (maxt > (MAXNUM/temp)) )
		goto error;

	a0 *= u;
	t = fabs(a0);

	/* terminating condition for asymptotic series */
	if( t > tlast )
		goto ndone;

	tlast = t;
	sum += alast;	/* the sum is one term behind */
	alast = a0;

	if( n > 200 )
		goto ndone;

	an += 1.0e0;
	bn += 1.0e0;
	n += 1.0e0;
	if( t > maxt )
		maxt = t;
	}
while( t > MACHEP );


pdone:	/* series converged! */

/* estimate error due to roundoff and cancellation */
*err = fabs(  MACHEP * (n + maxt)  );

alast = a0;
goto done;

ndone:	/* series did not converge */

/* The following "Converging factors" are supposed to improve accuracy,
 * but do not actually seem to accomplish very much. */

n -= 1.0;
x = 1.0/x;

switch( type )	/* "type" given as subroutine argument */
{
case 1:
	alast *= ( 0.5 + (0.125 + 0.25*b - 0.5*a + 0.25*x - 0.25*n)/x );
	break;

case 2:
	alast *= 2.0/3.0 - b + 2.0*a + x - n;
	break;

default:
	;
}

/* estimate error due to roundoff, cancellation, and nonconvergence */
*err = MACHEP * (n + maxt)  +  fabs ( a0 );


done:
sum += alast;
return( sum );

/* series blew up: */
error:
*err = MAXNUM;
mtherr( "hyperg", TLOSS );
return( sum );
}