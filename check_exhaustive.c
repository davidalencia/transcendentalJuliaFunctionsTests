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
#ifndef WITHOUT_OMP
#include <omp.h>
// #include "/usr/local/opt/libomp/include/omp.h"
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

// #ifndef HAVE_NO_SINCOS
// float sincos1f(float x)
// {
//   float s, c;
//   sincosf(x, &s, &c);
//   return s;
// }

// float sincos2f(float x)
// {
//   float s, c;
//   sincosf(x, &s, &c);
//   return c;
// }
// #endif

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

static float
cr_foo(float x)
{
  int inex;
  mpfr_t yy;
  float ret;
  mpfr_init2(yy, 24);
  mpfr_set_flt(yy, x, MPFR_RNDN);
  inex = MPFR_FOO(yy, yy, rnd2[rnd]);
  mpfr_subnormalize(yy, inex, rnd2[rnd]);
  ret = mpfr_get_flt(yy, MPFR_RNDN);
  mpfr_clear(yy);
  return ret;
}

/* Return the error in ulps between y and z,
   where y is the result computed by libm (for input x),
   and z is the result computed by MPFR.
   Both y and z should not be NaN.
   Only one of y and z is allowed to be infinite. */
static uint64_t
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
  return (err >= (float)UINTMAX_MAX) ? UINTMAX_MAX : (uint64_t)err;
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

#ifdef DEBUG
/* check if z is the correct rounding of x by computing another value
   with more precision and rounding it back */
static void
check_mpfr(float x, float z)
{
  mpfr_t zz;
  mpfr_prec_t prec2 = 100;
  int inex;
  mpfr_init2(zz, prec2);
  mpfr_set_flt(zz, x, MPFR_RNDN);
  MPFR_FOO(zz, zz, MPFR_RNDN);
  inex = mpfr_prec_round(zz, 24, rnd2[rnd]);
  mpfr_subnormalize(zz, inex, rnd2[rnd]);
  /* if inex=0, we can't conclude */
  if (inex != 0 && mpfr_get_flt(zz, MPFR_RNDN) != z)
  {
    fprintf(stderr, "Possible error in MPFR for x=%a\n", x);
    fprintf(stderr, "mpfr_%s (x) gives %a at precision 24\n", NAME, z);
    fprintf(stderr, "mpfr_%s (x) gives %a at precision 100\n", NAME,
            mpfr_get_flt(zz, MPFR_RNDN));
    fflush(stdout);
    exit(1);
  }
  mpfr_clear(zz);
}
#endif

uint64_t errors = 0;
uint64_t errors2 = 0; /* errors with 2 ulps or more */
uint64_t maxerr_u = 0;
unsigned int nmax = 0;
double maxerr = 0;

typedef union
{
  uint32_t n;
  float x;
} union_t;

float asfloat(uint32_t n)
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

#ifdef LIBMVEC
/* LIBMVEC should be 128 (sse4.2, default), or 256 (avx2) or 512 (avx512f) */
#if LIBMVEC == 256
#define LIBMVEC_N 8
#elif LIBMVEC == 512
#define LIBMVEC_N 16
#else /* default */
#define LIBMVEC_N 4
#endif
#endif

/* return FOO(x) */
static float
WRAPPER(float x)
{
#ifdef LIBMVEC
  float xx[LIBMVEC_N] = {
      x,
  },
        yy[LIBMVEC_N];
  for (int i = 0; i < LIBMVEC_N; i++)
    yy[i] = FOO(xx[i]);
  return yy[0];
#else
  return FOO(x);
#endif
}

/* e <- d with double exponent range */
static void
mpfr_set_d_safe(mpfr_t e, double d)
{
  mpfr_exp_t emax = mpfr_get_emax();
  mpfr_set_emax(1024);
  int ret = mpfr_set_d(e, d, MPFR_RNDN);
  mpfr_set_emax(emax);
  assert(ret == 0);
}

static void
print_maximal_error(unsigned int n)
{
  float x, y, z;
  mpfr_t e;
  uint64_t err;
  double err_double;

  x = asfloat(n);

  fesetround(rnd1[rnd]);
  y = WRAPPER(x); // resultado de la función que estamos probando
  z = cr_foo(x);
  err = ulp_error(y, z, x);
  err_double = ulp_error_double(y, x);

  mpfr_init2(e, 53);
  mpfr_set_d_safe(e, err_double);
  mpfr_printf("libm wrong by up to %.2RUe ulp(s) [%lu] for x=%a\n",
              e, err, x);
  printf("%s      gives %a\n", TOSTRING(STR), y);
  printf("mpfr_%s gives %a\n", TOSTRING(STR), z);
  fflush(stdout);
  mpfr_clear(e);
}

static void
check(unsigned int n)
{
  float x, y, z;

  x = asfloat(n);

  assert(!isnan(x));
  assert(!isinf(x));

#ifdef EXCLUDE
  if (exclude(x))
    return;
#endif

  fesetround(rnd1[rnd]);
  y = WRAPPER(x);
  z = cr_foo(x);

  if (y != z && !(isnan(y) && isnan(z)))
  {
    uint64_t err;
    double err_double;
#ifdef DEBUG
    if (!isinf(z))
      check_mpfr(x, z);
#endif
#ifndef WITHOUT_OMP
#pragma omp atomic update
#endif
    errors++;
    err = ulp_error(y, z, x);
    if (err > 1)
#ifndef WITHOUT_OMP
#pragma omp atomic update
#endif
      errors2++;
    err_double = ulp_error_double(y, x);
#ifndef WITHOUT_OMP
#pragma omp critical
#endif
    if (err > maxerr_u || (err == maxerr_u && err_double > maxerr))
    {
      maxerr_u = err;
      maxerr = err_double;
      nmax = n;
      if (verbose)
        print_maximal_error(nmax);
    }
  }
}

int main(int argc, char *argv[])
{
  unsigned int n;
  int nthreads = 0;
  while (argc >= 2)
  {
    if (strcmp(argv[1], "-v") == 0)
    {
      verbose++;
      argc--;
      argv++;
    }
    else if (strcmp(argv[1], "-rndn") == 0)
    {
      rnd = (mpfr_rnd_t)0;
      argc--;
      argv++;
    }
    else if (strcmp(argv[1], "-rndz") == 0)
    {
      rnd = (mpfr_rnd_t)1;
      argc--;
      argv++;
    }
    else if (strcmp(argv[1], "-rndu") == 0)
    {
      rnd = (mpfr_rnd_t)2;
      argc--;
      argv++;
    }
    else if (strcmp(argv[1], "-rndd") == 0)
    {
      rnd = (mpfr_rnd_t)3;
      argc--;
      argv++;
    }
    else if (argc >= 3 && strcmp(argv[1], "-nthreads") == 0)
    {
      nthreads = atoi(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else
    {
      fprintf(stderr, "Error, unknown option %s\n", argv[1]);
      exit(1);
    }
  }

#ifdef GLIBC
  printf("GNU libc version: %s\n", gnu_get_libc_version());
  printf("GNU libc release: %s\n", gnu_get_libc_release());
#endif
#ifdef __INTEL_COMPILER
  printf("Using Intel Math Library\n");
#endif
#ifdef AMD
  printf("Using AMD's library\n");
#endif
#ifdef NEWLIB
  printf("Using RedHat newlib\n");
  // __fdlib_version = -1; /* __fdlibm_ieee */
#endif
#ifdef OPENLIBM
  printf("Using OpenLibm\n");
#endif
#ifdef MUSL
  printf("Using Musl\n");
#endif
#ifdef LLVM
  printf("Using llvm-libc\n");
#endif
#ifdef _MSC_VER
  printf("Using Microsoft math library %d\n", _MSC_VER);
#endif
  printf("MPFR library: %-12s\nMPFR header:  %s (based on %d.%d.%d)\n",
         mpfr_get_version(), MPFR_VERSION_STRING, MPFR_VERSION_MAJOR,
         MPFR_VERSION_MINOR, MPFR_VERSION_PATCHLEVEL);
  printf("Checking function %s with %s\n", NAME,
         mpfr_print_rnd_mode(rnd2[rnd]));
  fflush(stdout);

#define MAXN 2139095040U
#ifndef WITHOUT_OMP
  /* Apparently Visual Studio does not properly set the number of threads. */
#pragma omp parallel
  if (nthreads <= 0)
    nthreads = omp_get_num_threads();
  if (verbose)
    printf("Using %d threads\n", nthreads);
  omp_set_num_threads(nthreads);
/* a dynamic schedule is better than 'guided' for example, especially for some
   functions like exp which is much faster for inputs yielding an overflow */
#pragma omp parallel for schedule(dynamic, 1024)
#endif
  for (n = 0; n < MAXN; n++)
  {
    /* we have to set emin/emax here, so that it is thread-local */
    mpfr_set_emin(-148);
    mpfr_set_emax(128);

    check(n);
    check(0x80000000 + n); /* negative values */
  }

  /* if verbose > 1, the maximal error was already printed */
  if (verbose == 0 && errors > 0)
    print_maximal_error(nmax);

  /* reset the rounding mode to nearest to print the %age below */
  fesetround(FE_TONEAREST);

  mpfr_t e;
  mpfr_init2(e, 53);
  mpfr_set_d_safe(e, maxerr);
  mpfr_printf("Total: errors=%lu (%.2f%%) errors2=%lu maxerr=%.2RUe ulp(s)\n",
              errors, 100.0 * (double)errors / (double)(2 * MAXN),
              errors2, e);
  mpfr_clear(e);
  mpfr_free_cache();
  fflush(stdout);
  return 0;
}
