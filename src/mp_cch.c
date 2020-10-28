/* search up segmentation error (Advanced R Programming Hadley Wickham)
 * lldb vs gdb will need to run on trig2 terminal via debugger
 * look at writing R extensions (specific to C part)
 * 
 * do compilation of mpfr on trig2 -lmpfr
 * the ./configure for gmp does not work on trig2 with is needed for mpft
 * ask zoyia to install things on trig2
 * 
 * PRIORITY: could run the code that throws error in the terminal
 * use a regular terminal window or rocker (docker) container (can also run interactively)
 * developers R-package-devel on CRAN (like a stack overflow) 
 */
/*--------------------------------------------------------------------------
 CONHYPTWO1(a,b,c,x,y)
 Phi1 hypergeometric function of two variables.
 Gradshteyn and Ryzhik (1965), equation 9.261.1
 This program assumes:  0<a<c, 0<b, y<=1.
 x is a column vector of real numbers, other parameters scalar.
 --------------------------------------------------------------------------*/
#include <math.h>
#include <Rmath.h>
#include <mpfr.h>
// #include "mex.h"
#define FACCURACY 1e-12

/*--------------------------------------------------------------------------
 HYPERG2F1
 Function for the Hypergeometric2F1(a,b,c,z).
 Covers cases:
 {0<b<c, 0<=z<=1} or {0<b, abs(c-a-b)<1, z<0}.
 --------------------------------------------------------------------------*/
double mphyperg2F1(double a, double b, double c, double z) {
  int j=0;
  mpfr_t w, F, res;
  
  mpfr_init2(F, 200);
  mpfr_set_d(F, 1.0, MPFR_RNDN);
  mpfr_init2(w, 200);
  mpfr_set_d(w, 1.0, MPFR_RNDN);
  mpfr_init2(res, 200);
  
  if (a<0) {
    /* Apply Abramowitz and Stegun 15.3.3 */
    mpfr_set_d(F, 1-z, MPFR_RNDN);
    mpfr_pow_si(F, F, c-a-b, MPFR_RNDN);
    mpfr_mul_d(F, F, mphyperg2F1(c-a,c-b,c,z), MPFR_RNDN);
  }
  else if (z<0) {
    /* Apply Abramowitz and Stegun 15.3.4 */
    mpfr_set_d(F, 1-z, MPFR_RNDN);
    mpfr_pow_si(F, F, -a, MPFR_RNDN);
    mpfr_mul_d(F, F, mphyperg2F1(a,c-b,c,z/(z-1)), MPFR_RNDN);
  }
  else if (z==1) {
    
    mpfr_set_d(F, lgammafn(c), MPFR_RNDN);
    mpfr_add_d(F, F, lgammafn(c-a-b)-lgammafn(c-a)-lgammafn(c-b), MPFR_RNDN);
    mpfr_exp(F, F, MPFR_RNDN);
  }
  else {
    while ((mpfr_get_d(w, MPFR_RNDN)/mpfr_get_d(F, MPFR_RNDN))>FACCURACY) {
      j++;
      
      mpfr_mul_d(w, w, ((a+j-1)*(b+j-1)/(c+j-1))*(z/j), MPFR_RNDN);
      mpfr_add(F, F, w, MPFR_RNDN);
    }
  }
  
  mpfr_set(res, F, MPFR_RNDN);
  
  mpfr_clear(w);
  mpfr_clear(F);
  
  return(mpfr_get_d(res, MPFR_RNDN));
}

/*--------------------------------------------------------------------------
 HyperTwo
 Function for Phi1(a,b,c,x,y)
 Assumes: 0<a<c, 0<b, y<1
 Use rule T5 if y<0.  Then use rule T1 if x>=0 or rule T2 if x<0.
 --------------------------------------------------------------------------*/
double mpHyperTwo(double a, double b, double c, double x, double y) {
  int m=0;
  mpfr_t F, zf, zfg, res;
  
  mpfr_init2(F, 200);
  mpfr_set_d(F, 0.0, MPFR_RNDN);
  mpfr_init2(zf, 200);
  mpfr_set_d(zf, 1.0, MPFR_RNDN);
  mpfr_init2(zfg, 200);
  mpfr_set_d(zfg, 0.0, MPFR_RNDN);
  mpfr_init2(res, 200);
  
  if (y<0) {
    mpfr_set_d(F, 1-y, MPFR_RNDN);
    mpfr_pow_si(F, F, -b, MPFR_RNDN);
    mpfr_mul_d(F, F, exp(x), MPFR_RNDN);
    mpfr_mul_d(F, F, mpHyperTwo(c-a,b,c,-x,y/(y-1)), MPFR_RNDN);
  }
  else {
    mpfr_set_d(zfg, mphyperg2F1(b,a,c,y), MPFR_RNDN);
    mpfr_set(F, zfg, MPFR_RNDN);
    if (x<0) {
      while (mpfr_get_d(zfg, MPFR_RNDN)/mpfr_get_d(F, MPFR_RNDN)>FACCURACY) {
        m++;
        mpfr_mul_d(zf, zf, ((c-a+m-1)/(c+m-1))*(-x/m), MPFR_RNDN);
        mpfr_mul_d(zfg, zf, mphyperg2F1(b,a,c+m,y), MPFR_RNDN);
        mpfr_add(F, F, zfg, MPFR_RNDN);
      }
      mpfr_mul_d(F, F, exp(x), MPFR_RNDN);
    }
    else {
      while (mpfr_get_d(zfg, MPFR_RNDN)/mpfr_get_d(F, MPFR_RNDN)>FACCURACY) {
        m++;
        mpfr_mul_d(zf, zf, ((a+m-1)/(c+m-1))*(x/m), MPFR_RNDN);
        mpfr_mul_d(zfg, zf, mphyperg2F1(b,a+m,c+m,y), MPFR_RNDN);
        mpfr_add(F, F, zfg, MPFR_RNDN);
      }
    }
  }
  
  mpfr_set(res, F, MPFR_RNDN);
  
  mpfr_clear(F);
  mpfr_clear(zf);
  mpfr_clear(zfg);
  
  return(mpfr_get_d(res, MPFR_RNDN));
}

/* check some syntax , also CLEAN AND REBUILD */

void mpphi1(double *a, double *b, double *c, double *x, double *y, double *phi, int *npara)
{
  int k;
  for (k = 0; k < *npara; k++) {
    //   if (x[k] <0) {
    /*  Since Linex system tends to report error for negative x, we
     use the following fomular to convert it to positive value
     1F1(a, b, x) = 1F1(b - a, b, -x) * exp(x) */
    /*      a[k] = b[k] - a[k];
     y[k] = hyperg(a[k], b[k], -x[k])*exp(x[k]);
     }
     else {
     y[k] = hyperg(a[k], b[k], x[k]);
     }*/
    phi[k] = mpHyperTwo(a[k], b[k], c[k], x[k], y[k]);
  }
}