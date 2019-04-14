/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
  
#include "arb_poly.h"
#include "acb_poly.h"
#include "arb_hypgeom.h"
#include "flint/profiler.h"
 
typedef arb_poly_struct *arb_poly_ptr;
typedef const arb_poly_struct *arb_poly_srcptr;

//procedure to calculate the absolute value of an arb_poly function (of length 1)
void
arb_poly_abs_series(arb_poly_t res, const arb_poly_t p, slong n, slong prec)
{
    arb_t a, one;
    arb_init(a);
    arb_init(one);

    arb_set_si(one, 1);

    arb_poly_evaluate(a, p, one, prec);

    arb_abs(a, a);

    arb_poly_set_arb(res, a);

    arb_clear(a);
    arb_clear(one);
}

//procedure to estabish the max value of the real parts of two arb_poly functions (of length 1)
void
arb_poly_max_series(arb_poly_t res, const arb_poly_t p, const arb_poly_t q, slong n, slong prec)
{
    arb_t a, b, one;
 
    arb_init(a);
    arb_init(b);
    arb_init(one);
    arb_set_si(one, 1);

    arb_poly_evaluate(a, p, one, prec);

    arb_poly_evaluate(b, q, one, prec);

    arb_max(a, a, b, prec);

    arb_poly_set_arb(res, a);

    arb_clear(a);
    arb_clear(b);
    arb_clear(one);
}

//procedure to calculate the modgamma-value
static void
_arb_poly_modgamma_series(arb_poly_t res, const arb_poly_t xN, const arb_poly_t y, slong n, slong prec)
{
    arb_poly_t a;
    arb_t ar;
 
    arb_poly_init(a);

    arb_init(ar);

    arb_set_ui(ar, 1005);
    arb_div_ui(ar, ar, 1000, prec);

    arb_poly_neg(a, y);
    arb_poly_pow_series(a, xN, a, n, prec);
	
    arb_poly_scalar_mul(res, a, ar, prec);
 
    arb_poly_clear(a);
    arb_clear(ar);
}
 
//procedure to calculate the sig exponent-value 
static void
_arb_poly_sig_series(arb_poly_t res, const arb_poly_t xN, const arb_poly_t t, const arb_poly_t y, slong n, slong prec)
{
    arb_poly_t a, b, c, d;

    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(d);
 
    arb_poly_add_si(a, y, 1, prec);
    arb_poly_scalar_mul_2exp_si(a, a, -1);

    arb_poly_log_series(b, xN, n, prec);
    arb_poly_mullow(b, b, t, n, prec);
    arb_poly_scalar_mul_2exp_si(b, b, -1);

    arb_poly_add_series(res, a, b, n, prec);

    arb_poly_set_si(c, 1);
    arb_poly_set_si(d, 1000);
    arb_poly_div_series(c, c, d, n, prec);

    arb_poly_sub_series(res, res, c, n, prec);

    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(d);
}

//procedure to establish some max values  
static void
lbound_numerator(arb_poly_t summandnd, const arb_poly_t ndsum_b, const arb_poly_t ndsum_a, arb_poly_t modgam, slong n, slong prec)
{
    arb_poly_t a, b, c, M1, M2;
    arb_poly_init(a); 
    arb_poly_init(b); 
    arb_poly_init(c);
    arb_poly_init(M1);
    arb_poly_init(M2);

    arb_poly_one(a);
    arb_poly_sub_series(b, a, modgam, n, prec);
    arb_poly_add_series(c, a, modgam, n, prec);
    arb_poly_div_series(a, b, c, n, prec);
    arb_poly_add_series(b, ndsum_b, ndsum_a, n, prec);
    arb_poly_abs_series(b, b, n, prec);
    arb_poly_mullow(M1, b, a, n, prec);

    arb_poly_sub_series(a, ndsum_b, ndsum_a, n, prec);
    arb_poly_abs_series(M2, a, n, prec);

    arb_poly_max_series(summandnd, M1, M2, n, prec);

    arb_poly_clear(a); 
    arb_poly_clear(b); 
    arb_poly_clear(c); 
    arb_poly_clear(M1);
    arb_poly_clear(M2);
}

//procedure to calculate the bt-value 
static void
_arb_poly_bt_series(arb_poly_t res, const arb_t logk, const arb_poly_t t, slong n, slong prec)
{
    arb_poly_scalar_mul(res, t, logk, prec);
    arb_poly_scalar_mul(res, res, logk, prec);
    arb_poly_scalar_mul_2exp_si(res, res, -2);
    arb_poly_exp_series(res, res, n, prec);
}

//procedure to calculate the at-value  
static void
_arb_poly_at_series(arb_poly_t res, const arb_t logk, const arb_poly_t t, const arb_poly_t y, slong n, slong prec)
{
    arb_poly_t p;
    arb_poly_init(p);

    arb_poly_scalar_mul(p, t, logk, prec);
    arb_poly_scalar_mul_2exp_si(p, p, -2);
    arb_poly_add_series(p, p, y, n, prec);
 
    arb_poly_scalar_mul(res, p, logk, prec);
    arb_poly_exp_series(res, res, n, prec);
 
    arb_poly_clear(p);
}

//procedure to establish the first n primes  
static void
_first_n_primes(ulong *primes, ulong n)
{
    ulong p = 2;
    ulong k = 0;
    while (k < n)
    {
        primes[k++] = p;
        p = n_nextprime(p, 1);
    }
}

//main procedure to establish the Lemma-bound  
void
arb_poly_abbeff_lbound(arb_poly_t lbound, const arb_poly_t t, const arb_poly_t y,
        arb_poly_srcptr moll, slong nprimes, const ulong *primes, const slong *d,
        slong D, const slong N1, const slong N2, slong n, slong prec)
{
    slong divs;

    arb_poly_t a, ndsum_a, ndsummand_a, ndsum_b, ndsummand_b, sumnd, summandnd;
    arb_poly_t atpoly, btpoly, Npoly;
    arb_poly_t kterm, modgam, sig;

    arb_poly_init(a);

    arb_poly_init(ndsum_a);
    arb_poly_init(ndsummand_a);
    arb_poly_init(ndsum_b);
    arb_poly_init(ndsummand_b);
    arb_poly_init(sumnd);
    arb_poly_init(summandnd);

    arb_poly_init(btpoly);
    arb_poly_init(atpoly);
    arb_poly_init(Npoly);

    arb_poly_init(kterm);
    arb_poly_init(modgam);
    arb_poly_init(sig);

    slong j, k, nd;
 
    divs = 1 << nprimes;

    arb_t logk, logj;
    arb_init(logk);
    arb_init(logj);

    //define common
    arb_poly_set_si(Npoly, N1);

    _arb_poly_sig_series(sig, Npoly, t, y, n, prec);
    _arb_poly_modgamma_series(modgam, Npoly, y, n, prec);

    //main loop from n=1 to D*N2
    arb_poly_zero(sumnd);
    for (k = 1; k <= D*N2; k++)
    {
        arb_poly_zero(ndsum_a);
        arb_poly_zero(ndsum_b);
        //sub loop from nd=1 to ndivs
        for (nd = 1; nd <= divs; nd++)
        {
            if (((k+1) % d[nd] == 0) && ((k+1) <= d[nd]*N2))
            {
                j = (k+1) / d[nd];
 
                arb_set_si(logj, j);
                arb_log(logj, logj, prec);
 
                _arb_poly_bt_series(btpoly, logj, t, n, prec);
                arb_poly_mullow(ndsummand_b, moll + nd, btpoly, n, prec);
                arb_poly_add_series(ndsum_b, ndsum_b, ndsummand_b, n, prec);

                _arb_poly_at_series(atpoly, logj, t, y, n, prec);
                arb_poly_mullow(ndsummand_a, moll + nd, atpoly, n, prec);
                arb_poly_mullow(ndsummand_a, ndsummand_a, modgam, n, prec);
                arb_poly_add_series(ndsum_a, ndsum_a, ndsummand_a, n, prec);
            }
        }

        if (!arb_poly_is_zero(ndsum_b) || !arb_poly_is_zero(ndsum_a))
        {
            arb_set_si(logk, k+1);
            arb_log(logk, logk, prec);
            arb_poly_scalar_mul(kterm, sig, logk, prec);
            arb_poly_neg(kterm, kterm);
            arb_poly_exp_series(kterm, kterm, n, prec);

            lbound_numerator(summandnd, ndsum_b, ndsum_a, modgam, n, prec);

            arb_poly_mullow(summandnd, kterm, summandnd, n, prec);
            arb_poly_add_series(sumnd, sumnd, summandnd, n, prec);
        }
    }

    //perform final calculations to complete lbound
    arb_poly_one(a);
    arb_poly_sub_series(a, a, modgam, n, prec);
    arb_poly_sub_series(lbound, a, sumnd, n, prec);

    arb_clear(logk);
    arb_clear(logj);
 
    arb_poly_clear(a);
    arb_poly_clear(ndsum_a);
    arb_poly_clear(ndsummand_a);
    arb_poly_clear(ndsum_b);
    arb_poly_clear(ndsummand_b);
    arb_poly_clear(sumnd);
    arb_poly_clear(summandnd);
    arb_poly_clear(btpoly);
    arb_poly_clear(atpoly);
    arb_poly_clear(Npoly);
    arb_poly_clear(kterm);
    arb_poly_clear(modgam);
    arb_poly_clear(sig);
}

//initialise and fill arrays for bt, at, moll. Process main loop  
void static
set_bound_constants(slong res, const arb_t t, const arb_t y, slong N1, slong N2, slong nprimes, slong digits, slong n, slong prec)
{

    slong D;

    slong i, j;

    arb_poly_ptr moll;
    slong divs;
    ulong *primes;
    slong *d;

    divs=0; D=0; 
 
    arb_t logp, coeff;
    arb_init(logp);
    arb_init(coeff);

    arb_poly_t tpoly, ypoly, btpoly, bound;
    arb_poly_init(tpoly);
    arb_poly_init(ypoly);
    arb_poly_init(btpoly);
    arb_poly_init(bound);
 
    arb_poly_set_arb(tpoly, t);
    arb_poly_set_arb(ypoly, y);
 
    divs = 1 << nprimes;
    primes = flint_malloc(nprimes * sizeof(*primes));
    _first_n_primes(primes, nprimes);
 
    D = 1;
    for (i = 0; i < nprimes; i++)
    {
        D *= primes[i];
    }
 
    d = flint_malloc((divs + 1) * sizeof(*d));
    moll = flint_malloc((divs + 1) * sizeof(*moll));
    for (i = 0; i < divs; i++)
    {
        arb_poly_init(moll + i + 1);
        arb_poly_one(moll + i + 1);
        d[i + 1] = 1;
        for (j = 0; j < nprimes; j++)
        {
            if (i & (1 << j))
            {
                d[i + 1] *= primes[j];
                arb_set_si(logp, primes[j]);
                arb_log(logp, logp, prec);
 
                _arb_poly_bt_series(btpoly, logp, tpoly, n, prec);
                arb_poly_mullow(
                        moll + i + 1, moll + i + 1, btpoly, n, prec);
                arb_poly_neg(moll + i + 1, moll + i + 1);
            }
        }
    }

    //this script was mainly written in arb_poly language, however with all polynomial lengths = 1
    arb_poly_abbeff_lbound(bound, tpoly, ypoly, moll, nprimes, primes, d, D, N1, N2, 1, prec);

    //print output
    printf("\n");
    printf("Lower Lemma bound ");
    printf("for t=");
    arb_printn(t, 1, ARB_STR_NO_RADIUS);
    printf(", y=");
    arb_printn(y, 1, ARB_STR_NO_RADIUS);
    printf(" in the range [N1=%ld,N2=%ld] ", N1, N2);
    printf("using %ld primes:", nprimes);
    printf("\n");
    printf("\n");
    arb_poly_get_coeff_arb(coeff, bound, 0);
    arf_printd(arb_midref(coeff), digits);
    printf("\n");
    printf("\n");

    //clear all variables/arrays/etc. 
    for (i = 0; i < divs; i++)
    {
        arb_poly_clear(moll + i + 1);
    }
    flint_free(d);
    flint_free(moll);
    flint_free(primes);

    arb_clear(logp);
    arb_clear(coeff);

    arb_poly_clear(tpoly);
    arb_poly_clear(ypoly);
    arb_poly_clear(btpoly);
    arb_poly_clear(bound);
}
 
//main initialisation - processing and validating input values 
int main(int argc, char *argv[])
{
    const char *t_str, *y_str, *N1_str, *N2_str, *m_str, *d_str;

    slong N1, N2, m, d, prec, res;

    int result = EXIT_SUCCESS;
    int usage_error = 0;
    res = 0;
 
    arb_t t, y;
    arb_init(t);
    arb_init(y);
 
    if (argc != 7)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }
 
    t_str  = argv[1];
    y_str  = argv[2];
    N1_str = argv[3];
    N2_str = argv[4];
    m_str  = argv[5];
    d_str  = argv[6];
 
    N1 = atol(N1_str);
    N2 = atol(N2_str);
    m  = atol(m_str);
    d  = atol(d_str);

    if (N1 > N2 || m < 0 || m > 20 || d < 1)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }
 
    prec = d * 3.32192809488736 + 30;

    arb_set_str(t, t_str, prec);
    arb_set_str(y, y_str, prec);

TIMEIT_ONCE_START

        set_bound_constants(res, t, y, N1, N2, m, d, 1, prec);

TIMEIT_ONCE_STOP 
 
finish:
 
    if (usage_error && result == EXIT_FAILURE)
    {
        flint_printf("Usage:\n");
        flint_printf("%s t y N1 N2 m d\n\n", argv[0]);
        flint_printf(
    "The script computes the lower bound for N in the range [N1,N2]\n"
    "for t,y. It uses an Euler mollification that includes m primes,\n"
    "so that for example m=0 corresponds to no mollification\n"
    "and m=3 corresponds to an Euler mollification that uses\n"
    "the first three primes {2, 3, 5}.\n"
    "The number of decimal digits to be displayed can be set by d."
    "\n");
    }
 
    arb_clear(t);
    arb_clear(y);;
 
    flint_cleanup();
    return result;
}
