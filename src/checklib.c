#include <stdint.h>
#include <math.h>
#include <mpfr.h>
#include <assert.h>
#include <string.h>

#include <stdio.h>

typedef int (*mpfr_fun)(mpfr_t, const mpfr_t, mpfr_rnd_t);

double ulp_distance(double x, mpfr_t yy)
{
  assert(!isinf(x));
  assert(!mpfr_inf_p(yy));
  int ret, sign;
  int prec_float = 53;
  double err;
  mpfr_t xx, dis;
  mpfr_prec_t prec = 256;
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
    e = -1023 - prec_float + 2;          // + 1 + 1 (los subnormales tienen un bit menos de precision)
  mpfr_mul_2si(dis, dis, -e, MPFR_RNDN); // dis = dis*2^(-e) RoundNearest = dis/2^e

  err = mpfr_get_d(dis, MPFR_RNDA); // err = (double)dis

  mpfr_clear(xx);
  mpfr_clear(dis);
  return err * sign;
}

double distance2inf64(double x)
{
  double steps = 0;
  x = fabs(x);
  printf("C] abs(x)=%a\n", x);
  while (!isinf(x))
  {
    x = nextafter(x, INFINITY);
    printf("C] x = step(x)=%a\n", x);

    steps++;
  }
  printf("steps= %f", steps);
  printf("\n");
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
    mpfr_init2(yy, 52);
    ret = mpfr_set_d(yy, x, MPFR_RNDN);
    assert(ret == 0);
    mpfr_foo(yy, yy, MPFR_RNDN);
    y = mpfr_get_d(yy, MPFR_RNDN);
    mpfr_clear(yy);
    return distance2inf64(y);
  }

  mpfr_init2(yy, 160);
  ret = mpfr_set_d(yy, x, MPFR_RNDN);
  assert(ret == 0);
  mpfr_foo(yy, yy, MPFR_RNDN);
  y = ulp_distance(y, yy);
  mpfr_clear(yy);
  return y;
}

// mpfr_foo get_mpfr_foo(char *strfoo)
// {
//   // if(strcmp(strfoo)==0){
//   //   return
//   // }
//   // else if(strcmp(strfoo)==0){

//   // }
//   return mpfr_cos;
// }
