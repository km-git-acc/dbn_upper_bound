/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
#include "acb_poly.h"
#include "flint/profiler.h"

//procedure to calculate a single stored sum for the stored sums matrix
static void
storedsumsinstance(slong res, arb_t X, slong N, slong taylorterms, slong index, slong prec)
{
    arb_t ab, bb, cb, x, n0, Narb, narb;

    arb_init(ab);
    arb_init(bb);
    arb_init(cb);
    arb_init(x);
    arb_init(n0);
    arb_init(Narb);
    arb_init(narb);

    acb_t a, b, nacb, expo, nexpo, sum;
    acb_init(a);
    acb_init(b);
    acb_init(nacb);
    acb_init(expo);
    acb_init(nexpo);
    acb_init(sum);

    slong d, i, j, n;

    i = index / taylorterms;
    j = index % taylorterms;

    //the exponent of i and j together
    d = i + 2*j;    

    //shift X
    arb_one(ab);
    arb_mul_2exp_si(ab, ab, -1);  
    arb_add(x, X, ab, prec); 

    //n0
    arb_set_si(Narb, N);
    arb_mul_2exp_si(n0, Narb, -1);

    //prepare the expo power
    acb_set_arb(a, x);
    acb_mul_onei(b, a);
    acb_sub_si(b, b, 1, prec);
    acb_mul_2exp_si(expo, b, -1);

    //produce the stored sum instance 'cell' of the matrix
    acb_zero(sum);
    for (n = 0; n < N; n++)
    {
        arb_set_ui(narb, n);
        arb_add_si(narb, narb, 1, prec);
        acb_set_arb(nacb, narb);

        arb_div(ab, narb, n0, prec);
        arb_log(ab, ab, prec);
        arb_pow_ui(ab, ab, d, prec);

        acb_pow(nexpo, nacb, expo, prec);
        acb_mul_arb(a, nexpo, ab, prec);

        acb_add(sum, sum, a, prec);
    }

    arb_fac_ui(ab, i, prec);
    arb_fac_ui(bb, j, prec);
    arb_mul(ab, ab, bb, prec);
    arb_set_si(cb, 4);
    arb_pow_ui(bb, cb, j, prec);
    arb_mul(ab, ab, bb, prec);
    acb_div_arb(a, sum, ab, prec);

    arb_printn(X, 30, ARB_STR_NO_RADIUS);
    printf(", %ld, %ld, %ld, %ld, ", N, taylorterms, i, j);
    acb_printn(a, 30, ARB_STR_NO_RADIUS);
    printf("\n");

    arb_clear(ab);
    arb_clear(bb);
    arb_clear(cb);
    arb_clear(x);
    arb_clear(n0);
    arb_clear(narb);
    arb_clear(Narb);

    acb_clear(a);
    acb_clear(b);
    acb_clear(nacb);
    acb_clear(expo);
    acb_clear(nexpo);
    acb_clear(sum);
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
    arb_t X, t0;
    arb_init(X);
    arb_init(t0);

    const char *index_str;
    slong index, size, taylorterms, N, prec, res;
    int result = EXIT_SUCCESS;
    res=0;

    if (argc != 2)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    index_str = argv[1];
    index = atol(index_str);

    //precision 
    prec = 20 * 3.32192809488736 + 60;

    //*** here are all fixed parameters set ***
    arb_set_str(X, "100000000000000037221", prec);
    taylorterms = 112;
    //*** here are all fixed parameters set ***
    
    size = taylorterms * taylorterms;

    if ((index < 0) || (index >= size))
        goto finish;

    arb_zero(t0);
    N = get_N(t0, X, prec);

TIMEIT_ONCE_START

    storedsumsinstance(res, X, N, taylorterms, index, prec);

TIMEIT_ONCE_STOP

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Required inputs:\n");
        flint_printf("%s index \n\n", argv[0]);
        flint_printf(
    "This script computes a single cell of the stored sums matrix that serves as \n"
    "input for 'Barrier-t-loop' calculations.\n" 
	"X = 10^20 + 37221; Taylor expansion = 112 x 112 terms (i.e. 20 digit accuracy).\n"
    "The only input required is a matrix indexnumber from 0 .. 12543 (112 x 112 - 1). \n"
    "Output: X, N, taylorterms, matrix row, matrix column, calculated value \n");
    }

    arb_clear(X);
    arb_clear(t0);
 
    flint_cleanup();

    return result;
}
