#include <math.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <mpfr.h>

extern double mphyperg(double, double, double), lgammafn(double);
double mploghyperg1F1_laplace(double, double, double);

double mploghyperg1F1(double a, double b, double x, int laplace)
{
  double y;
  
  if (laplace == 0) {
    if (x <0) {
      /*  Since Linux system tends to report error for negative x, we
       use the following fomular to convert it to positive value
       1F1(a, b, x) = 1F1(b - a, b, -x) * exp(x) */
      y = log(mphyperg(b-a, b, -x)) + x;
    }
    else {
      y = log(mphyperg(a, b, x));
    }
    
    if (y <0) {
      y = mploghyperg1F1_laplace(a, b, x);
    }
  }
  else {
    /* Laplace approximation assumes -x  for x positive
     if ( x <= 0.0 ){y = loghyperg1F1_laplace(a, b, -x); }
     else { y = loghyperg1F1_laplace(b - a, b, x) + x;} */
    y = mploghyperg1F1_laplace(a, b, x);
  }
  
  //    Rprintf("LOG Cephes 1F1(%lf, %lf, %lf) = %lf (%lf)\n", a,b,x,log(y), y);
  
  
  //    ly = hyperg1F1_laplace(a,b,x);
  //    Rprintf("called from hyperg1F1: LOG Pos 1F1(%lf, %lf, %lf) = %lf (%lf)\n", a,b,x,ly, exp(ly));
  if (!R_FINITE(y)  &&  laplace == 0) {
    warning("Cephes 1F1 function returned NA, using Laplace approximation");
    y = mploghyperg1F1_laplace(a, b, x);  // try Laplace approximation
  }
  
  
  
  return(y);
}

double mploghyperg1F1_laplace(double a, double b, double x)
{
  double mode, mode1, mode2, lprec;
  mpfr_t prec, logy;
  mpfr_t gabtemp, powtemp, sqrtemp, diff, insqrt, res;
  
  /* int u^(a-1) (1-u)^(b-1) exp(-x u) du   assuming that x >= 0 */
  
  mpfr_init2(prec, 200);
  mpfr_set_d(prec, 0.0, MPFR_RNDN);
  mpfr_init2(logy, 200);
  mpfr_set_d(logy, 0.0, MPFR_RNDN);
  
  mpfr_init2(gabtemp, 200);
  mpfr_set_d(gabtemp, -lgammafn(b) - lgammafn(a), MPFR_RNDN);  
    
  mpfr_init2(powtemp, 200);
  mpfr_set_d(powtemp, 0.0, MPFR_RNDN);
  mpfr_init2(sqrtemp, 200);
  mpfr_set_d(sqrtemp, 0.0, MPFR_RNDN);
  mpfr_init2(diff, 200);
  mpfr_set_d(diff, a-b-x, MPFR_RNDN);
  mpfr_init2(insqrt, 200);
  mpfr_set_d(insqrt, 0, MPFR_RNDN);
  mpfr_init2(res, 200);
  mpfr_set_d(res, 0, MPFR_RNDN);
  
  if ( x <= 0.0) {
    if (x < 0.0) {
      x = -x;
      
      mpfr_add_d(logy, gabtemp, lgammafn(a+b), MPFR_RNDN);
      
      //	mode = (2.0 - 2.0* a + b - x - sqrt(pow(b, 2.0) - 2.0*b*x + x*(4.0*(a-1.0)+x)))/
      //  (2*(a - 1.0 - b));
      
      mpfr_sqr(powtemp, diff, MPFR_RNDN);
      mpfr_add_d(insqrt, powtemp, 4.*a*b, MPFR_RNDN);
      mpfr_sqrt(sqrtemp, insqrt, MPFR_RNDN);
      
      mode1 = .5*(-a + b + x - mpfr_get_d(sqrtemp, MPFR_RNDN))/a;
      mode1 =  1.0/(1.0 + mode1);
      mode2 = .5*(-a + b + x + mpfr_get_d(sqrtemp, MPFR_RNDN))/a;
      mode2 =  1.0/(1.0 + mode2);
      if (a*log(mode1) + b*log(1.0 - mode1) - x*mode1  >
            a*log(mode2) + b*log(1.0 - mode2) - x*mode2) mode = mode1;
      else mode = mode2;
      //       Rprintf("mode 1 %lf, mode %lf\n", mode1, mode2);
      if (mode < 0) {
        mode = 0.0;
        warning("1F1 Laplace approximation on boundary\n");
      }
      else{
        /*prec = a*mode*(1.0 - mode) + (1.0-mode)*(1.0 - mode)*b +
         x*pow(1.0-mode, 3.0) - x*mode*(1.0 -
         mode)*(1.0-mode); */
        
        mpfr_set_d(prec, (1.0-mode)*(a + b - x)*pow(mode,2), MPFR_RNDN);
        mpfr_add_d(prec, prec, (1.0-mode)*mode*(a + b  + x), MPFR_RNDN);             
        
        if (mpfr_get_d(prec, MPFR_RNDN) > 0)  {
          lprec = log(mpfr_get_d(prec, MPFR_RNDN));
          
          mpfr_add_d(logy, logy, a*log(mode) + b*log(1.0 - mode) - x*mode, MPFR_RNDN);
          mpfr_add_d(logy, logy, -0.5*lprec + M_LN_SQRT_2PI, MPFR_RNDN);
          
        }
        else {mpfr_set_d(prec, 0.0, MPFR_RNDN);}
      }
      
      //	Rprintf("mode %lf prec %lf, Lap 1F1(%lf, %lf, %lf) = %lf\n", mode, prec, a,b,x, logy);
    }
    else {mpfr_set_d(logy, 0.0, MPFR_RNDN);}
  }
  else {
    mpfr_set_d(logy, x, MPFR_RNDN);
    mpfr_add_d(logy, logy, mploghyperg1F1_laplace(b - a, a, -x), MPFR_RNDN);
  }
  
  mpfr_set(res, logy, MPFR_RNDN);
  
  mpfr_clear(prec);
  mpfr_clear(logy);
  mpfr_clear(gabtemp);
  mpfr_clear(powtemp);
  mpfr_clear(sqrtemp);
  
  return(mpfr_get_d(res, MPFR_RNDN));
}


void mphypergeometric1F1(double *a, double *b, double *x, double *y, int *npara, int *Method)
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
    y[k] = mploghyperg1F1(a[k], b[k], x[k], Method[k]);
  }
}
