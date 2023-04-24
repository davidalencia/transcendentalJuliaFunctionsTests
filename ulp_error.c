/* Search worst cases of a univariate binary32 function, by exhaustive search.

   This program is open-source software distributed under the terms
   of the GNU General Public License <http://www.fsf.org/copyleft/gpl.html>.

   Compile with:

   $ gcc -DSTR=acos -O3 check_exhaustive.c -lmpfr -lgmp -lm -fopenmp
   $ icc -DSTR=acos -no-ftz -O3 check_exhaustive.c -lmpfr -lgmp -fopenmp

   Example of compilation command under Windows with MSVC:

   /permissive- /MP /ifcOutput "x64\Release\" /GS /GL /W3 /Gy /Zc:wchar_t
   /I"..\..\..\mpir\lib\x64\Release\" /I"..\..\..\mpfr\src" /I"..\"
   /Zi /Gm- /O2 /sdl /Fd"x64\Release\vc143.pdb" /Zc:inline /fp:strict
   /D "NDEBUG" /D "_CONSOLE" /D "NO_FLOAT128" /D "WORST" /D "USE_FLOAT" /D "UCRT" /D "DSTR=acos"
   /D "FOO=rndn" /D "_UNICODE" /D "UNICODE"
   /fp:except /errorReport:prompt /WX- /Zc:forScope /std:c17 /arch:SSE
   /Gd /Oi /MT /openmp- /FC /Fa"x64\Release\" /EHsc /nologo /Fo"x64\Release\"
   /Fp"x64\Release\binary32_exhaustive.pch" /diagnostics:column /openmp:llvm

   By default it uses all threads available. To use for example 32 threads:

   $ OMP_NUM_THREADS=32 ./a.out

   Options:
   -rndn: check rounding to nearest (default)
   -rndz: check rounding towards zero
   -rndu: check rounding towards +Inf
   -rndd: check rounding towards -Inf
   -v   : program is more verbose
   -nthreads nnn: uses nnn threads

   For NEWLIB: add -DNEWLIB (to avoid compilation error with __errno).


   clang -DSTR=acos -O3 check_exhaustive.c -lmpfr -lgmp -lm -Xpreprocesssor -I/opt/homebrew/Cellar/mpfr/4.2.0/include -L/opt/homebrew/Cellar/mpfr/4.2.0/lib  -I/opt/homebrew/Cellar/gmp/6.2.1_1/include -L/opt/homebrew/Cellar/gmp/6.2.1_1/lib
   pkg-config  --cflags --libs mpfr
*/

// #if !defined(__INTEL_COMPILER) && !defined(_GNU_SOURCE)
// #define _GNU_SOURCE
// #endif

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#ifdef _MSC_VER
/* The Microsoft library does not define j0f but _j0f,
   you should compile with -DSTR=_j0 */
#define mpfr__j0 mpfr_j0
#define _j0f _j0
#define mpfr__j1 mpfr_j1
#define _j1f _j1
#define mpfr__y0 mpfr_y0
#define _y0f _y0
#define mpfr__y1 mpfr_y1
#define _y1f _y1
#endif
#include <mpfr.h>
#include <assert.h>
#ifdef GLIBC
#include <gnu/libc-version.h>
#endif
#include <fenv.h>
#ifdef RLIBM
#include "rlibm.h"
#endif
#ifdef RLIBMALL
#include "float_rno_lib.h"
#endif
#ifdef VDT
#include "vdtMath.h"
#endif

#ifdef NEWLIB
/* Newlib defines different values for FE_DOWNWARD... */
#undef FE_DOWNWARD
#define FE_DOWNWARD 1
#undef FE_TONEAREST
#define FE_TONEAREST 0
#undef FE_TOWARDZERO
#define FE_TOWARDZERO 3
#undef FE_UPWARD
#define FE_UPWARD 2
/* RedHat's libm claims:
   undefined reference to `__errno' in j1f/y1f */
int errno;
int *__errno() { return &errno; }
#endif

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
#ifdef RLIBM
#ifndef FOO /* use -DFOO=rlibm_log2_8 or -DFOO=rlibm_log10_8 */
#define TMP CAT2(Rlibm_, STR)
#define FOO CAT2(TMP, FLOAT)
#endif
#endif
#ifdef RLIBMALL
#ifndef FOO
#define FOO CAT2(rlibm_all_fast_, STR)
#endif
#endif
#ifdef VDT
#define FOO_AUX CAT2(STR, FLOAT)
#define FOO CAT2(vdt::fast_, FOO_AUX)
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

/* Return the error in ulps between y and z,
   where y is the result computed by libm (for input x),
   and z is the result computed by MPFR.
   Both y and z should not be NaN.
   Only one of y and z is allowed to be infinite. */
static float
ulp_error(float y, float z, float x)
{
    float err, ulp;
    if (y == z)
        return 0;
    if (isinf(z))
    {
        mpfr_t zz;
        if (isinf(y)) /* then y and z are of different signs */
        {
            assert(y * z < 0);
            return UINTMAX_MAX;
        }
        /* we divide everything by 2, taking as reference the MPFR function */
        y = y / 2;
        mpfr_init2(zz, 24);
        mpfr_set_flt(zz, x, MPFR_RNDN);
        MPFR_FOO(zz, zz, rnd2[rnd]);
        // z = (float) (STR ((double) x) / 2);
        mpfr_div_2ui(zz, zz, 1, MPFR_RNDN);
        z = mpfr_get_flt(zz, MPFR_RNDN);
        /* If z is +Inf or -Inf, set it to +/-2^127 (since we divided y by 2) */
        if (isinf(z))
            z = (z > 0) ? 0x1p127 : -0x1p127;
        mpfr_clear(zz);
    }
    if (isinf(y))
    {
        assert(isinf(z) == 0);
        /* If the library gives +/-Inf but the correct rounding is in the
           binary32 range, assume the library gives +/-2^128. */
        z = z / 2; /* scale y and z by 1/2 */
        y = (y > 0) ? 0x1p127f : -0x1p127f;
    }
    err = y - z;
    ulp = nextafterf(z, y) - z;
    err = fabsf(err / ulp);
    return err;
    // return (err >= (float)UINTMAX_MAX) ? UINTMAX_MAX : (uint64_t)err;
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
    e = mpfr_get_exp(zz);
    mpfr_sub(zz, zz, yy, MPFR_RNDA);
    mpfr_abs(zz, zz, MPFR_RNDN);
    /* we should add 2^(e - prec - 1) to |zz| */
    mpfr_set_ui_2exp(yy, 1, e - prec - 1, MPFR_RNDN);
    mpfr_add(zz, zz, yy, MPFR_RNDA);
    /* divide by ulp(y) */
    e = (e - 24 < -149) ? -149 : e - 24;
    mpfr_mul_2si(zz, zz, -e, MPFR_RNDN);
    err = mpfr_get_d(zz, MPFR_RNDA);
    mpfr_set_emin(emin);
    mpfr_set_emax(emax);
    mpfr_clear(yy);
    mpfr_clear(zz);
    return err;
}
