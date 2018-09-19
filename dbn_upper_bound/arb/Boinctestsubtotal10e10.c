/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
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
storedsums(slong res, const arb_t X, const slong N, const slong expterms, const slong taylorterms,
           const slong indexsubtotal, const slong chunks, const slong prec)

{
    acb_mat_t finalmat, resacb;

    arb_mat_t expvec, tayvec, resarb;

    arb_t ab, bb, x, narb, Narb, accreq, n0, re, im, res1;
    arb_init(ab);
    arb_init(bb);
    arb_init(x);
    arb_init(narb);
    arb_init(Narb);
    arb_init(accreq);
    arb_init(n0);
    arb_init(re);
    arb_init(im);
    arb_init(res1);

    acb_t a, b, expo, nexpo, nacb;
    acb_init(a);
    acb_init(b);
    acb_init(expo);
    acb_init(nexpo);
    acb_init(nacb);
   
    arb_set_si(Narb, N);

    slong e, n, t;
    slong tasksize, startloop, endloop;

    //shift X
    arb_one(ab);
    arb_mul_2exp_si(ab, ab, -1);  
    arb_add(x, X, ab, prec); 

    //n0
    arb_mul_2exp_si(n0, Narb, -1);

    acb_mat_init(finalmat, expterms, taylorterms);
    acb_mat_init(resacb, expterms, taylorterms);

    arb_mat_init(expvec, expterms, 1);
    arb_mat_init(tayvec, 1, taylorterms);
    arb_mat_init(resarb, expterms, taylorterms);

    //prepare the expo power
    acb_set_arb(a, x);
    acb_mul_onei(b, a);
    acb_sub_si(b, b, 1, prec);
    acb_mul_2exp_si(expo, b, -1);

    //establish the start and endvalue of the subloop based on indexsubtotal
    tasksize = N / chunks;
    startloop = indexsubtotal * tasksize;
    
    if (indexsubtotal == chunks - 1)
    {
       endloop = N;
    }
    else
    {
       endloop = (indexsubtotal + 1) * tasksize;
    }

    printf("start: %ld end: %ld \n", startloop, endloop);

    //produce the stored sums matrix(expterms, taylorterms)
    for (n = startloop; n < endloop; n++)
    {
        arb_set_ui(narb, n);
        arb_add_si(narb, narb, 1, prec);
        acb_set_arb(nacb, narb);

        arb_div(ab, narb, n0, prec);
        arb_log(ab, ab, prec);

        matet(resarb, ab, expterms, taylorterms, expvec, tayvec, prec);

        acb_mat_set_arb_mat(resacb, resarb);

        acb_pow(nexpo, nacb, expo, prec);
        acb_mat_scalar_addmul_acb(finalmat, resacb, nexpo, prec);
    }

    //print the stored sums for storage in file
    arb_printn(X, 30, ARB_STR_NO_RADIUS);
    printf(", %ld, %ld, %d\n", expterms, taylorterms, 20);

    for (e = 0; e < expterms; e++)
    {
         for (t = 0; t < taylorterms; t++)
         {
              acb_get_real(re, acb_mat_entry(finalmat, e, t));
              acb_get_imag(im, acb_mat_entry(finalmat, e, t));
              arb_printn(re, 30, ARB_STR_NO_RADIUS);
              printf(", ");
              arb_printn(im, 30, ARB_STR_NO_RADIUS);
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
    arb_clear(x);
    arb_clear(narb);
    arb_clear(Narb);
    arb_clear(accreq);
    arb_clear(n0);
    arb_clear(re);
    arb_clear(im);
    arb_clear(res1);

    acb_mat_clear(finalmat);

    acb_clear(a);
    acb_clear(b);
    acb_clear(nacb);
    acb_clear(expo);
    acb_clear(nexpo);

}
 
int main(int argc, char *argv[])
{
    arb_t X;

    const char *indexsubtotal_str;

    slong N, prec, res, expterms, taylorterms, chunks, indexsubtotal;
    int result = EXIT_SUCCESS;
    res=0;

    arb_init(X);

    if (argc != 2)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    indexsubtotal_str = argv[1];
    indexsubtotal = atol(indexsubtotal_str);

    //****** all fixed parameters are set in this section ******
    prec = 20 * 3.32192809488736 + 100;

    arb_set_str(X, "60000155019", prec);
    N = 69098;

    expterms = 54;
    taylorterms = 54;

    chunks = 9;
    //***********************************************************

    if ((indexsubtotal < 0) || (indexsubtotal > chunks - 1))
        goto finish;

TIMEIT_ONCE_START
    storedsums(res, X, N, expterms, taylorterms, indexsubtotal, chunks, prec);
TIMEIT_ONCE_STOP

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Required inputs:\n");
        flint_printf("%s subtotalindex \n\n", argv[0]);
        flint_printf(
    "This script computes the subtotal of a single stored sums matrix that serves as input\n"
    "for 'Barrier-t-loop' calculations. All parameters are 'hard coded' in the relevant \n"
    "section in the script. The main loop of this script is split up into a number of 'chunks'\n"
    "(sub-loops) and the 'subtotalindex' (0..chunks-1) determines which 'chunk' will be computed.  \n");
    }

    arb_clear(X);
 
    flint_cleanup();

    return result;

}
