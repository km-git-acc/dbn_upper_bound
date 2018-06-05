/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
#include "acb_poly.h"
#include "acb_mat.h"

//procedure to estabish the max value of the real parts of two acb_poly functions
void
acb_poly_max_series(acb_poly_t res,
        const acb_poly_t p, const acb_poly_t q, slong n, slong prec)
{
    acb_t a, b, one;
    arb_t ar, br;
 
    arb_init(ar);
    arb_init(br);
    acb_init(a);
    acb_init(b);
    acb_init(one);
    acb_set_si(one, 1);

    acb_poly_evaluate(a, p, one, prec);
    acb_get_real(ar, a);

    acb_poly_evaluate(b, q, one, prec);
    acb_get_real(br, b);

    if (arb_gt(ar, br) > 0)
            {
                acb_poly_set(res, p);
            } else {
                acb_poly_set(res, q);
            }

    arb_clear(ar);
    arb_clear(br);
    acb_clear(a);
    acb_clear(b);
    acb_clear(one);
}

//procedure to establish all primes < n  
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

//procedure to calculate the absolute value of an acb_poly function
void
acb_poly_abs_series(acb_poly_t res,
        const acb_poly_t p, slong n, slong prec)
{
	arb_t abs;
    acb_t one, a; 
 
    arb_init(abs);
    acb_init(a);
    acb_init(one);

 
	acb_set_si(one, 1);
    acb_poly_evaluate(a, p, one, prec);
    acb_abs(abs, a, prec);

    acb_set_arb(a,abs);

    acb_poly_set_acb(res,a);
 
    arb_clear(abs);
    acb_clear(one);
    acb_clear(a);
}

void acb_poly_onei(acb_poly_t res)
{
    acb_poly_one(res);
    acb_onei(acb_poly_get_coeff_ptr(res, 0));
}
 
void
acb_poly_add_acb(acb_poly_t res,
        const acb_poly_t x, const acb_t y, slong prec)
{
    slong len = x->length;
 
    if (len == 0)
    {
        acb_poly_set_acb(res, y);
    }
    else
    {
        acb_poly_fit_length(res, len);
 
        acb_add(res->coeffs, x->coeffs, y, prec);
 
        if (res != x)
            _acb_vec_set(res->coeffs + 1, x->coeffs + 1, len - 1);
 
        _acb_poly_set_length(res, len);
        _acb_poly_normalise(res);
    }
}
 
void
acb_poly_add_arb(acb_poly_t res,
        const acb_poly_t x, const arb_t y, slong prec)
{
    acb_t z;
    acb_init(z);
    acb_set_arb(z, y);
    acb_poly_add_acb(res, x, z, prec);
    acb_clear(z);
}
 
 
void _acb_poly_alpha1_series(acb_poly_t z, const acb_poly_t s,
        slong n, slong prec)
{
    acb_poly_t a;
    acb_poly_init(a);
 
    acb_poly_zero(z);
 
    acb_poly_inv_series(a, s, n, prec);
    acb_poly_scalar_mul_2exp_si(a, a, -1);
    acb_poly_add_series(z, z, a, n, prec);
 
    acb_poly_add_si(a, s, -1, prec);
    acb_poly_inv_series(a, a, n, prec);
    acb_poly_add_series(z, z, a, n, prec);
 
    acb_poly_log_series(a, s, n, prec);
    acb_poly_scalar_mul_2exp_si(a, a, -1);
    acb_poly_add_series(z, z, a, n, prec);
 
    {
        arb_t c;
        arb_ptr u;
        arb_init(c);
        arb_const_log_sqrt2pi(c, prec);
        u = acb_realref(acb_poly_get_coeff_ptr(z, 0));
        arb_sub(u, u, c, prec);
        arb_clear(c);
    }
 
    acb_poly_clear(a);
}
 
void acb_poly_alpha1_series(acb_poly_t z, const acb_poly_t s,
        slong n, slong prec)
{
    if (z == s)
    {
        acb_poly_t u;
        acb_poly_init(u);
        _acb_poly_alpha1_series(u, s, n, prec);
        acb_poly_swap(z, u);
        acb_poly_clear(u);
    }
    else
    {
        _acb_poly_alpha1_series(z, s, n, prec);
    }
}
 
void _acb_poly_H01_series(acb_poly_t z, const acb_poly_t s,
        slong n, slong prec)
{
    arb_t c;
    acb_poly_t a, b;
    acb_t u;
    arb_ptr x;
 
    acb_poly_init(a);
    acb_poly_init(b);
 
    acb_init(u);
 
    arb_init(c);
    arb_const_log_sqrt2pi(c, prec);
 
    acb_poly_zero(z);
 
    acb_poly_add_si(a, s, -1, prec);
    acb_poly_log_series(a, a, n, prec);
    acb_poly_add_series(z, z, a, n, prec);
 
    acb_zero(u);
    x = acb_realref(u);
    arb_log_ui(x, 2, prec);
    arb_neg(x, x);
    arb_add_ui(x, x, 1, prec);
    arb_mul_2exp_si(x, x, -1);
    arb_add(x, x, c, prec);
    acb_poly_scalar_mul(a, s, u, prec);
    acb_poly_sub_series(z, z, a, n, prec);
 
    x = acb_realref(acb_poly_get_coeff_ptr(z, 0));
    arb_add(x, x, c, prec);
 
    acb_poly_add_si(b, s, 1, prec);
    acb_poly_scalar_mul_2exp_si(a, s, -1);
    acb_poly_log_series(a, a, n, prec);
    acb_poly_scalar_mul_2exp_si(a, a, -1);
    acb_poly_mullow(a, a, b, n, prec);
    acb_poly_add_series(z, z, a, n, prec);
 
    acb_poly_exp_series(z, z, n, prec);
 
    acb_poly_clear(a);
    acb_poly_clear(b);
    acb_clear(u);
    arb_clear(c);
}
 
void acb_poly_H01_series(acb_poly_t z, const acb_poly_t s,
        slong n, slong prec)
{
    if (z == s)
    {
        acb_poly_t u;
        acb_poly_init(u);
        _acb_poly_H01_series(u, s, n, prec);
        acb_poly_swap(z, u);
        acb_poly_clear(u);
    }
    else
    {
        _acb_poly_H01_series(z, s, n, prec);
    }
}
 
static void
_acb_poly_s_series(acb_poly_t res,
        const acb_poly_t x, const acb_poly_t y,
        slong n, slong prec)
{
    acb_poly_t ix;
 
    acb_poly_init(ix);
    acb_poly_onei(ix);
 
    acb_poly_mullow(ix, ix, x, n, prec);
 
    acb_poly_sub_series(res, ix, y, n, prec);
    acb_poly_add_si(res, res, 1, prec);
    acb_poly_scalar_mul_2exp_si(res, res, -1);
 
    acb_poly_clear(ix);
}
 
static void
_acb_poly_abbeff_summand_series(acb_poly_t res,
        const acb_t logk, const acb_poly_t a, const acb_poly_t b,
        slong n, slong prec)
{
    acb_poly_scalar_mul(res, a, logk, prec);
    acb_poly_sub_series(res, res, b, n, prec);
    acb_poly_scalar_mul(res, res, logk, prec);
    acb_poly_exp_series(res, res, n, prec);
}

void
abbeff_series(acb_t res, const arb_t rx, const arb_t ry, const arb_t rt, slong N, slong n, slong prec)
{
    acb_t one, logk;
    acb_poly_t x, y, t;
    acb_poly_t s, sc;
    acb_poly_t alph1, alph2;
    acb_poly_t hs, hsc;
    acb_poly_t ta1, ta2;
    acb_poly_t A0, B0;
    acb_poly_t a, b1, b2;
    acb_poly_t summand;
    acb_poly_t Asum, Bsum;
    acb_poly_t A, B;

    acb_init(one); 
    acb_poly_init(x); 
    acb_poly_init(y); 
    acb_poly_init(t); 
    //convert arb-inputs to acb_poly variables
    acb_poly_one(x);
    arb_set(acb_realref(acb_poly_get_coeff_ptr(x, 0)), rx);
    acb_poly_one(y);
    arb_set(acb_realref(acb_poly_get_coeff_ptr(y, 0)), ry);
    acb_poly_one(t);
    arb_set(acb_realref(acb_poly_get_coeff_ptr(t, 0)), rt);
 
    acb_poly_init(s);
    _acb_poly_s_series(s, x, y, n, prec);
 
    acb_poly_init(sc);
    acb_poly_add_si(sc, s, -1, prec);
    acb_poly_neg(sc, sc);
 
    acb_poly_init(alph1);
    acb_poly_alpha1_series(alph1, s, n, prec);
 
    acb_poly_init(alph2);
    acb_poly_alpha1_series(alph2, sc, n, prec);
 
    acb_poly_init(hs);
    acb_poly_H01_series(hs, s, n, prec);
 
    acb_poly_init(hsc);
    acb_poly_H01_series(hsc, sc, n, prec);
 
    acb_poly_init(ta1);
    acb_poly_mullow(ta1, t, alph1, n, prec);
 
    acb_poly_init(ta2);
    acb_poly_mullow(ta2, t, alph2, n, prec);
 
    acb_poly_init(A0);
    acb_poly_mullow(A0, ta1, alph1, n, prec);
    acb_poly_scalar_mul_2exp_si(A0, A0, -2);
    acb_poly_exp_series(A0, A0, n, prec);
    acb_poly_mullow(A0, A0, hs, n, prec);
    acb_poly_scalar_mul_2exp_si(A0, A0, -3);
 
    acb_poly_init(B0);
    acb_poly_mullow(B0, ta2, alph2, n, prec);
    acb_poly_scalar_mul_2exp_si(B0, B0, -2);
    acb_poly_exp_series(B0, B0, n, prec);
    acb_poly_mullow(B0, B0, hsc, n, prec);
    acb_poly_scalar_mul_2exp_si(B0, B0, -3);
 
    acb_poly_init(a);
    acb_poly_scalar_mul_2exp_si(a, t, -2);
 
    acb_poly_init(b1);
    acb_poly_scalar_mul_2exp_si(b1, ta1, -1);
    acb_poly_add_series(b1, b1, s, n, prec);
 
    acb_poly_init(b2);
    acb_poly_scalar_mul_2exp_si(b2, ta2, -1);
    acb_poly_add_series(b2, b2, sc, n, prec);
 
    acb_init(logk);
    acb_poly_init(summand);
    acb_poly_init(Asum);
    acb_poly_init(Bsum);
    {
        slong k;
        for (k = N; k >= 2; k--)
        {
            acb_set_si(logk, k);
            acb_log(logk, logk, prec);
 
            _acb_poly_abbeff_summand_series(summand, logk, a, b1, n, prec);
            acb_poly_add_series(Asum, Asum, summand, n, prec);
 
            _acb_poly_abbeff_summand_series(summand, logk, a, b2, n, prec);
            acb_poly_add_series(Bsum, Bsum, summand, n, prec);
        }
    }
    acb_poly_add_si(Asum, Asum, 1, prec);
    acb_poly_add_si(Bsum, Bsum, 1, prec);
 
    acb_poly_init(A);
    acb_poly_mullow(A, A0, Asum, n, prec);
 
    acb_poly_init(B);
    acb_poly_mullow(B, B0, Bsum, n, prec);
 
    acb_poly_div_series(a, A, B0, n, prec);
    acb_poly_add_series(a, a, Bsum, n, prec);

    acb_set_si(one, 1);
    acb_poly_evaluate(res, a, one, prec);
 
    acb_clear(one);
    acb_clear(logk);
    acb_poly_clear(s);
    acb_poly_clear(x);
    acb_poly_clear(y);
    acb_poly_clear(t);
    acb_poly_clear(sc);
    acb_poly_clear(alph1);
    acb_poly_clear(alph2);
    acb_poly_clear(hs);
    acb_poly_clear(hsc);
    acb_poly_clear(ta1);
    acb_poly_clear(ta2);
    acb_poly_clear(A0);
    acb_poly_clear(B0);
    acb_poly_clear(a);
    acb_poly_clear(b1);
    acb_poly_clear(b2);
    acb_poly_clear(summand);
    acb_poly_clear(Asum);
    acb_poly_clear(Bsum);
    acb_poly_clear(A);
    acb_poly_clear(B);
}

slong get_N(const arb_t t, const arb_t x, slong prec)
{
    arb_t pi, u;
    slong N;
    slong result;
 
    arb_init(pi);
    arb_init(u);
 
    arb_const_pi(pi, prec);
    arb_mul(u, pi, t, prec);
    arb_mul_2exp_si(u, u, -2);
    arb_add(u, u, x, prec);
    arb_div(u, u, pi, prec);
    arb_mul_2exp_si(u, u, -2);
    arb_sqrt(u, u, prec);
    arb_floor(u, u, prec);
 
    N = (slong) arf_get_d(arb_midref(u), ARF_RND_DOWN);
 
    if (arb_contains_si(u, N) &&
        !arb_contains_si(u, N-1) &&
        !arb_contains_si(u, N+1))
    {
        result = N;
    }
    else
    {
        fprintf(stderr, "Unexpected error: could not compute N\n");
        flint_abort();
    }
   
    arb_clear(pi);
    arb_clear(u);
 
    return result;
}
 
int main(int argc, char *argv[])
{

    arb_t x, y, t, a, b, X, eulabs, thres, one, d, re, im, out, y0, t0;
    acb_t cb, cc, cp, cone, eulprod, out1;
    const char *x_str, *numx_str, *nprimes_str, *thres_str, *sw_str, *y0_str, *t0_str ;
    slong N, prec, tht, i, v, nprimes, numx, sw;
    int result = EXIT_SUCCESS;
    int usage_error = 0;
 
    arb_mat_t offset;
    arb_mat_t output;
 
    arb_init(x);
    arb_init(y);
    arb_init(t);
    arb_init(a);
    arb_init(b);
    arb_init(X);
    arb_init(eulabs);
    arb_init(thres);
    arb_init(d);
    arb_init(re);
    arb_init(im);
    arb_init(out);
    arb_init(y0);
    arb_init(t0);
    arb_init(one);
    arb_one(one);

    acb_init(cb);
    acb_init(cc);
    acb_init(cp);
    acb_init(eulprod);
    acb_init(out1);
    acb_init(cone);
    acb_one(cone);

    arb_mat_init(offset, 3,1);
    arb_mat_init(output, 6,1);
 
    if (argc != 8)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }

    x_str = argv[1];
    numx_str = argv[2];
    nprimes_str = argv[3];
    thres_str = argv[4];
    y0_str = argv[5];
    t0_str = argv[6];
    sw_str = argv[7];

    prec = 128;

    arb_set_str(x, x_str, prec);
    numx = atol(numx_str);
    nprimes = atol(nprimes_str);
    arb_set_str(thres, thres_str, prec);
    sw = atol(sw_str);
    arb_set_str(y0, y0_str, prec);
    arb_set_str(t0, t0_str, prec);
	

    N = get_N(t, x, prec);

    ulong *primes;
    primes = flint_malloc(nprimes * sizeof(*primes));
    _first_n_primes(primes, nprimes);

    printf("\n");
    printf("Searching for Barrier location candidates around N = %ld ...\n", N);
    printf("\n");

    arb_zero(a);
    arb_set(arb_mat_entry(offset, 1, 0), a );
    arb_mul_2exp_si(a, one, -1);
    arb_set(arb_mat_entry(offset, 2, 0), a );
    arb_neg(a, a);
    arb_set(arb_mat_entry(offset, 0, 0), a );

    for (tht=0; tht <= numx; tht++)
    {
        arb_mul_2exp_si(a, one, -1);
        arb_set_si(b, tht);
        arb_add(b, b, a, prec);
        arb_add(X, x, b, prec);
        for (i=0; i<3; i++)
        {
            arb_mul_2exp_si(a, one, -1);
            arb_neg(a, a);
            arb_add(b, X, arb_mat_entry(offset, i, 0), prec);
            arb_mul_2exp_si(b, b, -1);
            arb_neg(b, b);
            acb_set_arb_arb(cb, a, b);
            acb_one(eulprod);
            for (v = 0; v < nprimes; v++)
            {
                acb_set_si(cp, primes[v]);
                acb_pow(cc, cp, cb, prec);
                acb_sub(cc, cone, cc, prec);
                acb_div(cc, cone, cc, prec);
                acb_mul(eulprod, eulprod, cc, prec);
            }
            acb_abs(eulabs, eulprod, prec);
            arb_set(arb_mat_entry(output, i, 0), eulabs);
        }

        for (i=0; i<3; i++)
        {
            arb_neg(a, one);
            arb_add(b, X, arb_mat_entry(offset, i, 0), prec);
            arb_mul_2exp_si(b, b, -1);
            arb_neg(b, b);
            acb_set_arb_arb(cb, a, b);
            acb_one(eulprod);
            for (v = 0; v < nprimes; v++)
            {
                acb_set_si(cp, primes[v]);
                acb_pow(cc, cp, cb, prec);
                acb_sub(cc, cone, cc, prec);
                acb_div(cc, cone, cc, prec);
                acb_mul(eulprod, eulprod, cc, prec);
            }
            acb_abs(eulabs, eulprod, prec);
            arb_set(arb_mat_entry(output, i+3, 0), eulabs);
        }
        if (arb_gt(arb_mat_entry(output, 3, 0), thres) && arb_gt(arb_mat_entry(output, 4, 0), thres)
		    && arb_gt(arb_mat_entry(output, 5, 0), thres))
        { 
            arb_printn(X, 30, ARB_STR_NO_RADIUS);
            printf(", ");
            for (i=0; i<6; i++)
            {
               arb_printn(arb_mat_entry(output, i, 0), 10, ARB_STR_NO_RADIUS);
               printf(", ");
            }
            printf("\n");
            if (sw == 1)
            {
                for (i=0; i<3; i++)
                {
                   arb_zero(t);
                   arb_set(a, y0);
                   arb_add(b, X, arb_mat_entry(offset, i, 0), prec);
                   printf("-> ");
                   arb_printn(b, 30, ARB_STR_NO_RADIUS); 
                   printf(", ");
                   arb_printn(a, 10, ARB_STR_NO_RADIUS);  
                   printf(", ");
                   arb_printn(t, 10, ARB_STR_NO_RADIUS);  
                   printf(", ");                   
                   abbeff_series(out1, b, a, t, N, 1, prec);
                   acb_abs(out,out1,prec);
                   arb_printn(out, 10, ARB_STR_NO_RADIUS);
                   printf("\n");
                }
                for (i=0; i<3; i++)
                {
                   arb_one(a);
                   arb_set(t, t0);
                   arb_add(b, X, arb_mat_entry(offset, i, 0), prec);
                   printf("-> ");
                   arb_printn(b, 30, ARB_STR_NO_RADIUS); 
                   printf(", ");
                   arb_printn(a, 10, ARB_STR_NO_RADIUS);  
                   printf(", ");
                   arb_printn(t, 10, ARB_STR_NO_RADIUS);  
                   printf(", ");                   
                   abbeff_series(out1, b, a, t, N, 1, prec);
                   acb_abs(out,out1,prec);
                   arb_printn(out, 10, ARB_STR_NO_RADIUS);
                   printf("\n");
                }

            } 
		}
    }

    flint_free(primes);

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Usage:\n");
        flint_printf("%s x xnum nprimes thres y0 t0 sw d\n\n", argv[0]);
        flint_printf(
    "This script assists in finding an optimal locations to place\n"
    "a Barrier. Input values are x for the search start, numx for\n"
    "the search length, nprimes for the number of primes in the Euler \n"
    "product. Thres sets the threshold for last three (y=1) outputs. \n"
    "When sw=0 it just prints 6 fields (y=0 and y=1), when sw=1 it also\n"
    "outputs the associated abbeff values for y=y0, y=1 and t=0, t=t0.\n"
    "                         \n");
    }
 
    arb_mat_clear(offset);
    arb_mat_clear(output);
 
    arb_clear(x);
    arb_clear(y);
    arb_clear(t);
    arb_clear(a);
    arb_clear(b);
    arb_clear(X);
    arb_clear(eulabs);
    arb_clear(one);
    arb_clear(thres);
    arb_clear(d);
    arb_clear(re);
    arb_clear(im);
    arb_clear(y0);
    arb_clear(t0);
    arb_clear(out);

    acb_clear(cb);
    acb_clear(cc);
    acb_clear(cp);
    acb_clear(cone);
    acb_clear(eulprod);
    acb_clear(out1);
 
    flint_cleanup();

    return result;
}
