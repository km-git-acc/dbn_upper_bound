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

//procedure to establish required number of exp and taylorterms
static void
ttermmagnitude(arb_t res, const ulong expterms, const ulong taylorterms, arb_t x1, arb_t x2, arb_t preterm, slong prec)
{
    arb_t a, b, c, d, f, iarb, one, half, zerodot2, borderisum, borderjsum, bordermax;

    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(f);
    arb_init(iarb);
    arb_init(one);
    arb_init(zerodot2);
    arb_init(half);
    arb_init(borderisum);
    arb_init(borderjsum);
    arb_init(bordermax);
    arb_one(one);

    ulong i, j;

    arb_mul_2exp_si(half, one, -1);
    arb_set_si(a, 2);
    arb_set_si(b, 10);
    arb_div(zerodot2, a, b, prec);
    arb_mul_2exp_si(half, one, -1);

    arb_zero(borderjsum);
    for (j = 0; j < taylorterms; j++)
    {
        arb_pow_ui(a, x1, expterms, prec);
        arb_fac_ui(f, expterms, prec);
        arb_div(a, a, f, prec);
        arb_pow_ui(b, x2, j, prec);
        arb_fac_ui(f, j, prec);
        arb_div(b, b, f, prec);
        arb_pow_ui(c, half, expterms, prec);
        arb_pow_ui(d, zerodot2, j, prec);
        arb_mul(b, a, b, prec);
        arb_mul(b, c, b, prec);
        arb_mul(b, d, b, prec);
        arb_add(borderjsum, borderjsum, b, prec);
    }

    arb_zero(borderisum);
    for (i = 0; i < expterms; i++)
    {
        arb_pow_ui(a, x1, i, prec);
        arb_fac_ui(f, i, prec);
        arb_div(a, a, f, prec);
        arb_pow_ui(b, x2, taylorterms, prec);
        arb_fac_ui(f, taylorterms, prec);
        arb_div(b, b, f, prec);
        arb_pow_ui(c, half, i, prec);
        arb_pow_ui(d, zerodot2, taylorterms, prec);
        arb_mul(b, a, b, prec);
        arb_mul(b, c, b, prec);
        arb_mul(b, d, b, prec);
        arb_add(borderisum, borderisum, b, prec);
    }

    arb_add(bordermax, borderisum, borderjsum, prec);
    arb_mul(res, bordermax, preterm, prec);

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(f);
    arb_clear(iarb);
    arb_clear(one);
    arb_clear(half);
    arb_clear(borderisum);
    arb_clear(borderjsum);
    arb_init(bordermax);
}


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
storedsums(slong res, const arb_t X, const arb_t y0, 
           const slong H, slong digits, slong n, slong prec)

{
    acb_mat_t finalmatb, finalmata, resacb;

    arb_mat_t expvec, tayvec, resarb;

    arb_t ab, bb, harb, Harb, x, y, x1, x2, preterm, n0, re, im, zerodot25, res1, accreq;
    arb_init(ab);
    arb_init(bb);
    arb_init(harb);
    arb_init(Harb);
    arb_init(x);
    arb_init(y);
    arb_init(x1);
    arb_init(x2);
    arb_init(preterm);
    arb_init(n0);
    arb_init(re);
    arb_init(im);
    arb_init(zerodot25);
    arb_init(res1);
    arb_init(accreq);

    acb_t a, b, one, expob, expoa, hexpob, hexpoa, hacb;
    acb_init(a);
    acb_init(b);
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
    arb_set_si(zerodot25, 25);
    arb_set_si(ab, 100);
    arb_div(zerodot25, zerodot25, ab, prec);

    slong e, h, t, expterms, taylorterms;

    //change X and y to midpoints
    arb_one(ab);
    arb_mul_2exp_si(ab, ab, -1);  
    arb_add(x, X, ab, prec); 
    arb_add_si(ab, y0, 1, prec);
    arb_mul_2exp_si(y, ab, -1);

    //n0
    arb_mul_2exp_si(n0, Harb, -1);

    //establish required expterms and taylorterms
    arb_log(x1, n0, prec);
    arb_pow_ui(x2, x1, 2, prec);
    arb_mul_2exp_si(x2, x2, -2);
    arb_pow(preterm, n0, zerodot25, prec);
    arb_mul_ui(preterm, preterm, H, prec);

    arb_set_si(ab, digits);
    arb_neg(ab, ab);
    arb_set_si(bb, 10);
    arb_pow(accreq, bb, ab, prec);

    expterms = 3;
    taylorterms = 3;
    arb_set_si(res1, 1000);

    while (arb_ge(res1, accreq))
    {
          expterms = expterms + 1;
          taylorterms = taylorterms + 1;  
          ttermmagnitude(res1, expterms, taylorterms, x1, x2, preterm, prec);
    }
    expterms = expterms + 1;
    taylorterms = taylorterms + 1;

    acb_mat_init(finalmatb, expterms, taylorterms);
    acb_mat_init(finalmata, expterms, taylorterms);
    acb_mat_init(resacb, expterms, taylorterms);

    arb_mat_init(expvec, expterms, 1);
    arb_mat_init(tayvec, 1, taylorterms);
    arb_mat_init(resarb, expterms, taylorterms);

    //prepare the expob and expoa powers
    acb_set_arb(a, x);
    acb_mul_onei(b, a);
    acb_sub_si(b, b, 1, prec);
    acb_sub_arb(b, b, y, prec);
    acb_mul_2exp_si(expob, b, -1);

    acb_mul_onei(b, a);
    acb_neg(b, b);
    acb_sub_si(b, b, 1, prec);
    acb_add_arb(b, b, y, prec);
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

    //print the stored sums for storage in file
    arb_printn(X, 30, ARB_STR_NO_RADIUS);
    printf(", ");
    arb_printn(y0, 10, ARB_STR_NO_RADIUS);
    printf(", %ld, %ld, %ld\n", expterms, taylorterms, digits);

    for (e = 0; e < expterms; e++)
    {
         for (t = 0; t < taylorterms; t++)
         {
              acb_get_real(re, acb_mat_entry(finalmatb, e, t));
              acb_get_imag(im, acb_mat_entry(finalmatb, e, t));
              arb_printn(re, digits + 10, ARB_STR_NO_RADIUS);
              printf(", ");
              arb_printn(im, digits + 10, ARB_STR_NO_RADIUS);
              if (t < taylorterms - 1) 
                  printf(", ");
              else
                  printf("\n");
         }
    }
    printf("***\n");
    for (e = 0; e < expterms; e++)
    {
         for (t = 0; t < taylorterms; t++)
         {
              acb_get_real(re, acb_mat_entry(finalmata, e, t));
              acb_get_imag(im, acb_mat_entry(finalmata, e, t));
              arb_printn(re, digits+10, ARB_STR_NO_RADIUS);
              printf(", ");
              arb_printn(im, digits+10, ARB_STR_NO_RADIUS);

              if (t < taylorterms - 1) 
                  printf(", ");
              else
                  printf("\n");
         }
    }

    arb_mat_clear(expvec);
    arb_mat_clear(tayvec);
    acb_mat_clear(resacb);
    arb_mat_clear(resarb);

    arb_clear(ab);
    arb_clear(bb);
    arb_clear(harb);
    arb_clear(Harb);
    arb_clear(x);
    arb_clear(y);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(preterm);
    arb_clear(n0);
    arb_clear(re);
    arb_clear(im);
    arb_clear(zerodot25);
    arb_clear(res1);
    arb_clear(accreq);

    acb_mat_clear(finalmatb);
    acb_mat_clear(finalmata);

    acb_clear(a);
    acb_clear(b);
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
    arb_t a, b, X, y0, t0;

    const char *X_str, *y0_str, *digits_str;
    slong N, prec, digits, res;
    int result = EXIT_SUCCESS;
    res=0;

    arb_init(a);
    arb_init(b);
    arb_init(X);
    arb_init(y0);
    arb_init(t0);

    if (argc != 4)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    X_str = argv[1];
    y0_str = argv[2];
    digits_str = argv[3];

    //precision 
    digits = atol(digits_str);
    prec = digits * 3.32192809488736 + 60;

    arb_set_str(X, X_str, prec);
    arb_set_str(y0, y0_str, prec);
    arb_zero(t0);

    N = get_N(t0, X, prec);

    storedsums(res, X, y0, N, digits, 1, prec);

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Required inputs:\n");
        flint_printf("%s X y0 digits \n\n", argv[0]);
        flint_printf(
    "This script computes the stored sums as input for 'Barrier-t-loop',\n"
    "calculations. Precision will be set by a fixed formula based on digits\n"
    "and automatic selection of the required number of Taylor expansion terms." 
    "It first outputs X, y0, Expterms, Taylorterms, digits and then" 
    "the sums stored in the finalmatb and finalmata matrices separated by ***.\n"
    "t0 has been fixed at 0.\n");
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(X);
    arb_clear(y0);
    arb_clear(t0);
 
    flint_cleanup();

    return result;
}
