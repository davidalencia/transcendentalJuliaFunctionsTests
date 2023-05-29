
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <mpfr.h>
#include <assert.h>
#include <fenv.h>

#define CAT1(X, Y) X##Y
#define CAT2(X, Y) CAT1(X, Y)
#ifndef MPFR_FOO
#define MPFR_FOO CAT2(mpfr_, STR)
#endif

// ----------------------------------------------
// global vars
const float Inf = 1.0 / 0.0; // positive infinte

// ----------------------------------------------
// type defs
typedef union // need
{
    uint64_t n;
    double x;
} union_t;

// ----------------------------------------------
// functions
double asdouble(uint64_t n) // need
{
    union_t u;
    u.n = n;

    return u.x;
}

static double
ulp_error_double(double y, double x)
{
    mpfr_t yy, zz;
    mpfr_prec_t prec = 160;
    int ret;
    mpfr_exp_t e;
    double err;
    mpfr_exp_t emin = mpfr_get_emin();
    mpfr_exp_t emax = mpfr_get_emax();

    mpfr_set_emin(mpfr_get_emin_min());
    mpfr_set_emax(mpfr_get_emax_max());
    mpfr_init2(yy, 53);
    mpfr_init2(zz, prec);
    if (!isinf(y))
    {
        ret = mpfr_set_d(yy, y, MPFR_RNDN);
        assert(ret == 0);
    }
    else
        mpfr_set_ui_2exp(yy, 1, 128, MPFR_RNDN);
    ret = mpfr_set_d(zz, x, MPFR_RNDN);
    assert(ret == 0);
    // mpfr_printf("xx =  %.80Rf \n ", zz);
    MPFR_FOO(zz, zz, MPFR_RNDN);
    // mpfr_printf("x = %f \n y = %f \n yy = %.80Rf \n zz = %.80Rf \n ", x, y, yy, zz);

    if (mpfr_cmp(yy, zz) == 0)
    {
        mpfr_clear(yy);
        mpfr_clear(zz);
        return (double)0;
    }

    e = mpfr_get_exp(zz);
    mpfr_sub(zz, zz, yy, MPFR_RNDA); // zz = zz - yy
    mpfr_abs(zz, zz, MPFR_RNDN);     // zz = abs(zz)
    // mpfr_printf("dif = %.80Rf \n ", zz);

    // /* we should add 2^(e - prec - 1) to |zz| */
    // printf("e: %ld\n", e);
    int eint = (int)e;
    // printf("e: %d\n", eint);

    mpfr_set_ui_2exp(yy, 1, e - prec - 1, MPFR_RNDN); // yy = 1*2^(e-prec-1), RoundNearest
    mpfr_add(zz, zz, yy, MPFR_RNDA);                  // zz = zz + yy, Round away from zero.
    // mpfr_printf("dif2 = %.80Rf \n ", zz);

    // /* divide by ulp(y) */
    // e_min = -1022 (f32 -126), p=52 (f32 p=24),          tabla 3.13 HFPA
    // e_min-p+1 = -126-24+1= -149
    // e -p+1 = e-24+1 = e-23
    eint = (eint - 52 < -1073) ? -1073 : eint - 52;

    mpfr_mul_2si(zz, zz, -eint, MPFR_RNDN); // zz = zz*2^(-e), esta usando la definición de Ulp de Goldberg
    // mpfr_printf("divided by ulps = %.80Rf \n ", zz);

    err = mpfr_get_d(zz, MPFR_RNDA);
    mpfr_set_emin(emin);
    mpfr_set_emax(emax);
    mpfr_clear(yy);
    mpfr_clear(zz);
    return err;
}

double distance2inf32(float x)
{
    mpfr_t zz;
    mpfr_init2(zz, 52); // 52 = precision of normal float64
    int ret = mpfr_set_d(zz, Inf, MPFR_RNDN);
    assert(ret == 0);

    MPFR_FOO(zz, zz, MPFR_RNDN); // mpfr trig function
    ret = mpfr_inf_p(zz);        // Return non-zero if op is an infinity
    double steps = 0;
    if (mpfr_sgn(zz)) // Return a positive value if op > 0, zero if op = 0, and a negative
        while (ret == 0)
        {
            mpfr_nextabove(zz);
            ret = mpfr_inf_p(zz);
            steps++;
        }
    else
        while (ret == 0)
        {
            mpfr_nextbelow(zz);
            ret = mpfr_inf_p(zz);
            steps++;
        }
    printf("steps= %f\n", steps);
    return (double)steps;
}

double check(long unsigned int n, double (*WRAPPER)(const double))
{
    mpfr_set_emin(-148);
    mpfr_set_emax(128);
    double x, y;

    x = asdouble(n);

    assert(!isnan(x));
    assert(!isinf(x));

    fesetround(FE_TONEAREST);
    y = WRAPPER(x);

    if (isinf(y))
        return distance2inf32(x);

    return ulp_error_double(y, x);
}

// int main(int argc, char const *argv[])
// {
//     // long unsigned int x = 0x400921fb54442c46; // pi-ish                                    i
//     long unsigned int x = 0x3ff921fb54442d16; // pi-ish /2
//     // long unsigned int x = 0x3ff0000000000000; // 1
//     // long unsigned int x = 0x0; // 0
//     // long unsigned int x = 0x3fe41b089a027525; // 0.6283
//     printf("%f", check(x, tan));
//     return 0;
// }
