#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <mpfr.h>
#include <assert.h>
#ifndef WITHOUT_OMP
#include <omp.h>
// #include "/usr/local/opt/libomp/include/omp.h"
#endif
#include <fenv.h>

#ifdef __APPLE__
/* Apple defines __exp10f, added before the name was standardized in C */
#define exp10f __exp10f
#define sincos __sincosf
#endif

/* redefine mpfr_lgamma to a function without the "int s" parameter,
   to match the lgamma function (thanks Vincent Lefèvre) */
static inline int
real_mpfr_lgamma(mpfr_t y, int *s, mpfr_t x, mpfr_rnd_t r)
{
  return mpfr_lgamma(y, s, x, r);
}

#undef mpfr_lgamma
#define mpfr_lgamma my_mpfr_lgamma

int my_mpfr_lgamma(mpfr_t y, mpfr_t x, mpfr_rnd_t r)
{
  int s;
  return real_mpfr_lgamma(y, &s, x, r);
}

int mpfr_sincos1(mpfr_t y, mpfr_t x, mpfr_rnd_t r)
{
  mpfr_t z;
  int inex;
  mpfr_init2(z, mpfr_get_prec(y));
  inex = mpfr_sin_cos(y, z, x, r);
  mpfr_clear(z);
  inex = inex % 4;
  return (inex == 0) ? 0 : (inex == 1) ? 1
                                       : -1;
}

int mpfr_sincos2(mpfr_t y, mpfr_t x, mpfr_rnd_t r)
{
  mpfr_t z;
  int inex;
  mpfr_init2(z, mpfr_get_prec(y));
  inex = mpfr_sin_cos(z, y, x, r);
  mpfr_clear(z);
  inex = inex / 4;
  return (inex == 0) ? 0 : (inex == 1) ? 1
                                       : -1;
}

/* https://stackoverflow.com/questions/1489932/how-to-concatenate-twice-with-the-c-preprocessor-and-expand-a-macro-as-in-arg */
#define FLOAT f
#define CAT1(X, Y) X##Y
#define CAT2(X, Y) CAT1(X, Y)
#if !defined(RLIBM) && !defined(RLIBMALL) && !defined(VDT)
#define FOO CAT2(STR, FLOAT)
#endif
#ifndef MPFR_FOO
#define MPFR_FOO CAT2(mpfr_, STR)
#endif
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define NAME TOSTRING(FOO)

int rnd1[] = {FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD};
mpfr_rnd_t rnd2[] = {MPFR_RNDN, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD};

mpfr_rnd_t rnd = MPFR_RNDN; /* default is to nearest */

int verbose = 0;

/* tgamma in C corresponds to mpfr_gamma */
int mpfr_tgamma(mpfr_t y, mpfr_t x, mpfr_rnd_t r)
{
  return mpfr_gamma(y, x, r);
}

/* return the ulp error between y and FOO(x), where FOO(x) is computed with
   MPFR with 100 bits of precision */
static double
ulp_error_double(float y, float x)
{
  mpfr_t yy, zz;
  mpfr_prec_t prec = 100;
  int ret;
  mpfr_exp_t e;
  double err;
  mpfr_exp_t emin = mpfr_get_emin();
  mpfr_exp_t emax = mpfr_get_emax();

  mpfr_set_emin(mpfr_get_emin_min());
  mpfr_set_emax(mpfr_get_emax_max());
  mpfr_init2(yy, 24);
  mpfr_init2(zz, prec);
  if (!isinf(y))
  {
    ret = mpfr_set_flt(yy, y, MPFR_RNDN);
    assert(ret == 0);
  }
  else
    mpfr_set_ui_2exp(yy, 1, 128, MPFR_RNDN);
  ret = mpfr_set_flt(zz, x, MPFR_RNDN);
  assert(ret == 0);
  MPFR_FOO(zz, zz, MPFR_RNDN);
  // printf("x:%f, y:", x);
  // printMPFR(yy);
  // printf("\n");

  e = mpfr_get_exp(zz);
  mpfr_sub(zz, zz, yy, MPFR_RNDA); // zz = zz - yy
  mpfr_abs(zz, zz, MPFR_RNDN);     // zz = abs(zz)
  /* we should add 2^(e - prec - 1) to |zz| */
  mpfr_set_ui_2exp(yy, 1, e - prec - 1, MPFR_RNDN); // yy = 1*2^(e-prec-1), RoundNearest
  mpfr_add(zz, zz, yy, MPFR_RNDA);                  // zz = zz + yy, Round away from zero.
  /* divide by ulp(y) */
  e = (e - 24 < -149) ? -149 : e - 24;
  mpfr_mul_2si(zz, zz, -e, MPFR_RNDN); // zz = zz*2^(-e*2)
  err = mpfr_get_d(zz, MPFR_RNDA);
  mpfr_set_emin(emin);
  mpfr_set_emax(emax);
  mpfr_clear(yy);
  mpfr_clear(zz);
  return err;
}

uint64_t errors = 0;
uint64_t errors2 = 0; /* errors with 2 ulps or more */
uint64_t maxerr_u = 0;
unsigned int nmax = 0;
double maxerr = 0;

typedef union // need
{
  uint32_t n;
  float x;
} union_t;

float asfloat(uint32_t n) // need
{
  union_t u;
  u.n = n;
  return u.x;
}

uint32_t
asuint(float x)
{
  union_t u;
  u.x = x;
  return u.n;
}

double check(unsigned int n, double (*WRAPPER)(const double))
{
  mpfr_set_emin(-148);
  mpfr_set_emax(128);
  float x, y;

  x = asfloat(n);

  assert(!isnan(x));
  assert(!isinf(x));

  fesetround(rnd1[rnd]);
  y = WRAPPER(x);
  // z = cr_foo(x);

  double err_double = ulp_error_double(y, x);
  return err_double;
}
