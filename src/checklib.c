#include <stdint.h>
#include <math.h>
#include <mpfr.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

typedef int (*mpfr_fun)(mpfr_t, const mpfr_t, mpfr_rnd_t);
mpfr_prec_t prec = 100;

double ulp_distance(double x, mpfr_t yy)
{
  assert(!isinf(x));
  assert(!mpfr_inf_p(yy));

  if(x==0.0){ 
    mpfr_exp_t e = 1074;
    mpfr_mul_2si(yy, yy, e, MPFR_RNDN);
    return mpfr_get_d(yy, MPFR_RNDA); 
  }

  int ret, sign;
  int prec_float = 53;
  double err;
  mpfr_t xx, dis;
  mpfr_exp_t e;

  mpfr_set_emin(mpfr_get_emin_min());
  mpfr_set_emax(mpfr_get_emax_max());
  mpfr_init2(dis, prec);
  mpfr_init2(xx, prec_float);

  ret = mpfr_set_d(xx, x, MPFR_RNDN);
  assert(ret == 0);

  if (mpfr_cmp(xx, yy) == 0)
  {
    mpfr_clear(xx);
    mpfr_clear(dis);
    return (double)0;
  }

  e = mpfr_get_exp(xx);
  mpfr_sub(dis, yy, xx, MPFR_RNDA); // dis = yy - xx, RoundAway

  sign = mpfr_sgn(dis);
  mpfr_abs(dis, dis, MPFR_RNDN); // dis = abs(dis), RoundNearest

  /* we should add 2^(e - prec - 1) to |dis| */
  mpfr_set_ui_2exp(yy, 1, e - prec - 1, MPFR_RNDN); // yy = 1*2^(e-prec-1), RoundNearest
  mpfr_add(dis, dis, yy, MPFR_RNDA);                // dis = dis + yy, RoundAway .

  /* divide by ulp(y) */
  if (isnormal(x))
    e = e - prec_float;
  else
    e = -1023 - prec_float + 2;          // + 1 + 1 (los subnormales tienen un bit menos de precisión)
  mpfr_mul_2si(dis, dis, -e, MPFR_RNDN); // dis = dis*2^(-e) RoundNearest = dis/2^e

  err = mpfr_get_d(dis, MPFR_RNDA); // err = (double)dis

  mpfr_clear(xx);
  mpfr_clear(dis);
  return err * sign;
}

double ulp_distance32(float x, mpfr_t yy)
{
  assert(!isinf(x));
  assert(!mpfr_inf_p(yy));
  
  if(x==0.0){ 
    mpfr_exp_t e = 149;
    mpfr_mul_2si(yy, yy, e, MPFR_RNDN);
    return mpfr_get_d(yy, MPFR_RNDA); 
  }

  int ret, sign;
  int prec_float = 24;
  double err;
  mpfr_t xx, dis;
  mpfr_exp_t e;

  mpfr_set_emin(mpfr_get_emin_min());
  mpfr_set_emax(mpfr_get_emax_max());
  mpfr_init2(dis, prec);
  mpfr_init2(xx, prec_float);

  ret = mpfr_set_flt(xx, x, MPFR_RNDN);
  assert(ret == 0);

  if (mpfr_cmp(xx, yy) == 0)
  {
    mpfr_clear(xx);
    mpfr_clear(dis);
    return (double)0;
  }

  e = mpfr_get_exp(xx);
  mpfr_sub(dis, yy, xx, MPFR_RNDA); // dis = yy - xx, RoundAway

  sign = mpfr_sgn(dis);
  mpfr_abs(dis, dis, MPFR_RNDN); // dis = abs(dis), RoundNearest

  /* we should add 2^(e - prec - 1) to |dis| */
  mpfr_set_ui_2exp(yy, 1, e - prec - 1, MPFR_RNDN); // yy = 1*2^(e-prec-1), RoundNearest
  mpfr_add(dis, dis, yy, MPFR_RNDA);                // dis = dis + yy, RoundAway .

  /* divide by ulp(y) */
  if (isnormal(x))
    e = e - prec_float;
  else
    e = -127 - prec_float + 2;           // + 1 + 1 (los subnormales tienen un bit menos de precision)
  mpfr_mul_2si(dis, dis, -e, MPFR_RNDN); // dis = dis*2^(-e) RoundNearest = dis/2^e
  err = mpfr_get_d(dis, MPFR_RNDA);      // err = (double)dis

  mpfr_clear(xx);
  mpfr_clear(dis);
  return err * sign;
}

double distance2inf64(double x)
{
  double steps = 0;
  x = fabs(x);
  while (!isinf(x))
  {
    steps++;
    x = nextafter(x, INFINITY);
  }
  return steps;
}

double distance2inf32(float x)
{
  double steps = 0;
  x = fabs(x);
  while (!isinf(x))
  {
    steps++;
    x = nextafterf(x, INFINITY);
  }
  return steps;
}

double ulp_error(double (*foo)(const double), mpfr_fun mpfr_foo, double x)
{
  int ret;
  double y;
  mpfr_t yy;

  y = foo(x);

  if (isinf(y))
  {
    mpfr_init2(yy, 53);
    mpfr_set_d(yy, x, MPFR_RNDN);
    mpfr_foo(yy, yy, MPFR_RNDN);
    y = mpfr_get_d(yy, MPFR_RNDN);
    mpfr_clear(yy);
    return distance2inf64(y);
  }

  mpfr_init2(yy, prec);
  ret = mpfr_set_d(yy, x, MPFR_RNDN);
  assert(ret == 0);
  mpfr_foo(yy, yy, MPFR_RNDN);
  y = ulp_distance(y, yy);
  mpfr_clear(yy);
  return y;
}

double ulp_error32(float (*foo)(const float), mpfr_fun mpfr_foo, float x)
{
  int ret;
  float y;
  double err;
  mpfr_t yy;

  y = foo(x);

  if (isinf(y))
  {
    mpfr_init2(yy, 24);
    mpfr_set_flt(yy, x, MPFR_RNDN);
    mpfr_foo(yy, yy, MPFR_RNDN);
    y = mpfr_get_flt(yy, MPFR_RNDN);
    mpfr_clear(yy);
    return distance2inf32(y);
  }

  mpfr_init2(yy, prec);
  ret = mpfr_set_flt(yy, x, MPFR_RNDN);
  assert(ret == 0);
  mpfr_foo(yy, yy, MPFR_RNDN);
  err = ulp_distance32(y, yy);
  mpfr_clear(yy);
  return err;
}

mpfr_fun get_mpfr_fun(char *strfoo)
{
  if (strcmp(strfoo, "cos") == 0)
  {
    return mpfr_cos;
  }
  else if (strcmp(strfoo, "sin") == 0)
  {
    return mpfr_sin;
  }
  else if (strcmp(strfoo, "tan") == 0)
  {
    return mpfr_tan;
  }
  else if (strcmp(strfoo, "cospi") == 0)
  {
    return mpfr_cospi;
  }
  else if (strcmp(strfoo, "sinpi") == 0)
  {
    return mpfr_sinpi;
  }
  else if (strcmp(strfoo, "acos") == 0)
  {
    return mpfr_acos;
  }
  else if (strcmp(strfoo, "asin") == 0)
  {
    return mpfr_asin;
  }
  else if (strcmp(strfoo, "atan") == 0)
  {
    return mpfr_atan;
  }
  else if (strcmp(strfoo, "csc") == 0)
  {
    return mpfr_csc;
  }
  else if (strcmp(strfoo, "sec") == 0)
  {
    return mpfr_sec;
  }
  else if (strcmp(strfoo, "cot") == 0)
  {
    return mpfr_cot;
  }
  else if (strcmp(strfoo, "cosh") == 0)
  {
    return mpfr_cosh;
  }
  else if (strcmp(strfoo, "sinh") == 0)
  {
    return mpfr_sinh;
  }
  else if (strcmp(strfoo, "tanh") == 0)
  {
    return mpfr_tanh;
  }
  else if (strcmp(strfoo, "acosh") == 0)
  {
    return mpfr_acosh;
  }
  else if (strcmp(strfoo, "asinh") == 0)
  {
    return mpfr_asinh;
  }
  else if (strcmp(strfoo, "atanh") == 0)
  {
    return mpfr_atanh;
  }
  else if (strcmp(strfoo, "exp") == 0)
  {
    return mpfr_exp;
  }
  else if (strcmp(strfoo, "expm1") == 0)
  {
    return mpfr_expm1;
  }
  else if (strcmp(strfoo, "exp2") == 0)
  {
    return mpfr_exp2;
  }
  else if (strcmp(strfoo, "exp10") == 0)
  {
    return mpfr_exp10;
  }
  else if (strcmp(strfoo, "log") == 0)
  {
    return mpfr_log;
  }
  else if (strcmp(strfoo, "log2") == 0)
  {
    return mpfr_log2;
  }
  else if (strcmp(strfoo, "log10") == 0)
  {
    return mpfr_log10;
  }
  else if (strcmp(strfoo, "log1p") == 0)
  {
    return mpfr_log1p;
  }
  printf("function %s not found\n", strfoo);
  return mpfr_log1p;
}

double check(double (*foo)(const double), char *mpfrStr, double x)
{
  mpfr_fun mpfr_foo = get_mpfr_fun(mpfrStr);
  return ulp_error(foo, mpfr_foo, x);
}
