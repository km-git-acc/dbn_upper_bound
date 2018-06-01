/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
#include "acb_poly.h"
#include "acb_mat.h"
#include "flint/profiler.h"

//procedure to fill the matet matrix
static void
matet(arb_mat_t resarb, const arb_t xval, const slong expterms, const slong taylorterms, 
           arb_mat_t expvec, arb_mat_t tayvec, slong prec)
{
    arb_t a, b, eacb, tacb;

    slong e, t;

    arb_init(a);
    arb_init(b);
    arb_init(eacb);
    arb_init(tacb);

    arb_set_si(arb_mat_entry(expvec, 0, 0), 1);
    for (e = 0; e < expterms-1; e++)
    {
        arb_set_ui(eacb, e);
        arb_add_si(b, eacb, 1, prec);
        arb_div(a, xval, b, prec);
        arb_mul(b, a, arb_mat_entry(expvec, e, 0), prec);
        arb_set(arb_mat_entry(expvec, e+1, 0), b);
    }

    arb_set_si(arb_mat_entry(tayvec, 0, 0), 1);
    arb_mul(a, xval, xval, prec);
    arb_mul_2exp_si(a, a, -2);
    for (t = 0; t < taylorterms-1; t++)
    {
        arb_set_ui(tacb, t);
        arb_add_si(b, tacb, 1, prec);
        arb_div(b, a, b, prec);
        arb_mul(b, b, arb_mat_entry(tayvec, 0, t), prec);
        arb_set(arb_mat_entry(tayvec, 0, t+1), b);
    }

    arb_mat_mul(resarb, expvec, tayvec, prec);

    arb_clear(a);
    arb_clear(b);
    arb_clear(eacb);
    arb_clear(tacb);
}

void 
storedsums(slong res, const acb_t cX, const acb_t cy, 
           const slong H, const slong taylorterms, const slong expterms, slong n, slong prec)

{
    acb_mat_t finalmatb, finalmata, resacb;
    acb_mat_init(finalmatb, expterms, taylorterms);
    acb_mat_init(finalmata, expterms, taylorterms);
    acb_mat_init(resacb, expterms, taylorterms);

    arb_mat_t expvec, tayvec, resarb;
    arb_mat_init(expvec, expterms, 1);
    arb_mat_init(tayvec, 1, taylorterms);
    arb_mat_init(resarb, expterms, taylorterms);

    arb_t ab, harb, Harb, n0;
    arb_init(ab);
    arb_init(harb);
    arb_init(Harb);
    arb_init(n0);

    acb_t a, b, X, y, one, expob, expoa, hexpob, hexpoa, hacb;
    acb_init(a);
    acb_init(b);
    acb_init(X);
    acb_init(y);
    acb_init(expob);
    acb_init(expoa);
    acb_init(hexpob);
    acb_init(hexpoa);
    acb_init(hacb);
    acb_init(one);
    acb_one(one);

    acb_poly_t finpolyb, finpolya;
    acb_poly_init(finpolyb);
    acb_poly_init(finpolya);
   
    arb_set_si(Harb, H);

    printf("\n");
    printf("Calculating the stored sums for N= %ld.\n", H);
    printf("\n");

    slong h;

    //change X and y to midpoints
    acb_mul_2exp_si(a, one, -1);  
    acb_add(X, cX, a, prec); 
    acb_add_si(a, cy, 1, prec);
    acb_mul_2exp_si(y, a, -1);

    //n0
    arb_mul_2exp_si(n0, Harb, -1);

    //prepare the expob and expoa powers
    acb_mul_onei(b, X);
    acb_sub_si(b, b, 1, prec);
    acb_sub(b, b, y, prec);
    acb_mul_2exp_si(expob, b, -1);

    acb_mul_onei(b, X);
    acb_neg(b, b);
    acb_sub_si(b, b, 1, prec);
    acb_add(b, b, y, prec);
    acb_mul_2exp_si(expoa, b, -1);

    //produce the stored sums matrix(expterms, taylorterms)
    for (h = 0; h < H; h++)
    {
        arb_set_ui(harb, h);
        arb_add_si(harb, harb, 1, prec);
        acb_set_arb(hacb, harb);

        arb_div(ab, harb, n0, prec);
        arb_log(ab, ab, prec);

        matet(resarb, ab, expterms, taylorterms, expvec, tayvec, prec);

        acb_mat_set_arb_mat(resacb, resarb);

        acb_pow(hexpob, hacb, expob, prec);
        acb_pow(hexpoa, hacb, expoa, prec);
        acb_mat_scalar_addmul_acb(finalmatb, resacb, hexpob, prec);
        acb_mat_scalar_addmul_acb(finalmata, resacb, hexpoa, prec);
    }

    acb_mat_printd(finalmatb, 100);
    printf("\n");
    printf("\n");
    printf("\n");
    acb_mat_printd(finalmata, 100);

    arb_mat_clear(expvec);
    arb_mat_clear(tayvec);
    acb_mat_clear(resacb);
    arb_mat_clear(resarb);

    arb_clear(ab);
    arb_clear(harb);
    arb_clear(Harb);
    arb_clear(n0);

    acb_clear(a);
    acb_clear(b);
    acb_clear(X);
    acb_clear(y);
    acb_clear(one);
    acb_clear(hacb);
    acb_clear(expob);
    acb_clear(expoa);
    acb_clear(hexpob);
    acb_clear(hexpoa);
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
    arb_t X, y0, t0;
    acb_t cX, cy0;

    const char *X_str, *y0_str, *expterms_str, *taylorterms_str;
    slong N, prec, expterms, taylorterms, res;
    int result = EXIT_SUCCESS;
    res=0;

    arb_init(X);
    arb_init(y0);
    arb_init(t0);

    acb_init(cX);
    acb_init(cy0);

    if (argc != 5)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    X_str = argv[1];
    y0_str = argv[2];
    expterms_str = argv[3];
    taylorterms_str = argv[4];

    expterms = atol(expterms_str);
	taylorterms = atol(taylorterms_str);

    //precision set after calibration with other software to produce at least 12 digits accuracy for all X <= 10^15
    prec = 400;

    arb_set_str(X, X_str, prec);
    arb_set_str(y0, y0_str, prec);
    acb_set_arb(cX, X);
    acb_set_arb(cy0, y0);
    arb_zero(t0);

    N = get_N(t0, X, prec);

TIMEIT_ONCE_START

    storedsums(res, cX, cy0, N, taylorterms, expterms, 1, prec);
	acb_printd(cy0,10);

TIMEIT_ONCE_STOP

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Required inputs:\n");
        flint_printf("%s X y0 expterms taylorterms \n\n", argv[0]);
        flint_printf(
    "This script computes the stored sums as input '3D-Barrier',\n"
    "calculations. Precision has been set at 100 digits accuracy.\n"
    "It prints the finalmatb and finalmata sums separated by 3 lines.\n"
    "t0 has been fixed at 0.\n");
    }

    arb_clear(X);
    arb_clear(y0);
    arb_clear(t0);

    acb_clear(cX);
    acb_clear(cy0);
 
    flint_cleanup();

    return result;
}
