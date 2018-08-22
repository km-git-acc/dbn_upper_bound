/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_calc.h"
#include "flint/profiler.h"
#include "pthread.h"

//define some global parameters to be accessed by individual threads
acb_mat_t finmat1;
acb_mat_t finmat2;
acb_mat_t finmat3;
acb_mat_t finmat4;
acb_mat_t finmat5;
acb_mat_t finmat6;
acb_mat_t finmat7;
acb_mat_t finmat8;
acb_mat_t finmat9;
acb_mat_t finmat10;

arb_t XX;

struct ThreadData {
slong expterms, taylorterms, prec, n00, start, stop, id;
}; 


// est1 vecsum
static void
est1_vecsum(arb_t estsum1, arb_t var, const arb_t n0, slong expterms, slong taylorterms, slong prec)
{
    arb_t a, b, c, d, abslog, zerodot2, zerodot66;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(abslog);
    arb_init(zerodot2);
    arb_init(zerodot66);   

    ulong j;

    //common values
    arb_set_si(a, 2);
    arb_set_si(b, 10);
    arb_div(zerodot2, a, b, prec);
    arb_set_si(a, 66);
    arb_set_si(b, 100);
    arb_div(zerodot66, a, b, prec);

    arb_div(a, var, n0, prec);
    arb_log(a, a, prec);
    arb_abs(abslog, a);

    arb_pow_ui(a, abslog, expterms, prec);
    arb_fac_ui(b, expterms, prec);   
    arb_div(d, a, b, prec);

    arb_zero(estsum1);
    for (j = 0; j <= taylorterms; j++)
    {
        arb_pow_ui(b, abslog, 2, prec);
        arb_mul_2exp_si(b, b, -2);
        arb_pow_ui(b, b, j, prec);
        arb_fac_ui(c, j, prec);   
        arb_div(b, b, c, prec);
        arb_mul(b, b, d, prec);

        arb_pow_ui(a, zerodot66, expterms, prec);
        arb_mul(b, b, a, prec);
        arb_pow_ui(c, zerodot2, j, prec);
        arb_mul(b, b, c, prec);
 
        arb_add(estsum1, estsum1, b, prec);
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(abslog);
    arb_clear(zerodot2);
    arb_clear(zerodot66);
}

// est2 vecsum
static void
est2_vecsum(arb_t estsum2, arb_t var, const arb_t n0, slong expterms, slong taylorterms, slong prec)
{
    arb_t a, b, c, d, abslog, zerodot2, zerodot66;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(abslog);
    arb_init(zerodot2);
    arb_init(zerodot66);   

    ulong i;

    //common values
    arb_set_si(a, 2);
    arb_set_si(b, 10);
    arb_div(zerodot2, a, b, prec);
    arb_set_si(a, 66);
    arb_set_si(b, 100);
    arb_div(zerodot66, a, b, prec);

    arb_div(a, var, n0, prec);
    arb_log(a, a, prec);
    arb_abs(abslog, a);

    arb_pow_ui(b, abslog, 2, prec);
    arb_mul_2exp_si(b, b, -2);
    arb_pow_ui(a, b, taylorterms, prec);
    arb_fac_ui(c, taylorterms, prec);   
    arb_div(d, a, c, prec);

    arb_zero(estsum2);
    for (i = 0; i <= expterms; i++)
    {
        arb_pow_ui(b, abslog, i, prec);
        arb_fac_ui(c, i, prec);   
        arb_div(b, b, c, prec);
        arb_mul(b, b, d, prec);

        arb_pow_ui(a, zerodot66, i, prec);
        arb_mul(b, b, a, prec);
        arb_pow_ui(c, zerodot2, taylorterms, prec);
        arb_mul(b, b, c, prec);
 
        arb_add(estsum2, estsum2, b, prec);
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(abslog);
    arb_clear(zerodot2);
    arb_clear(zerodot66);
}

// prepare the integrand for est1
int
f_est1_integrand(acb_ptr estint1, const acb_t z, void *param, slong order, slong prec)
{
    arb_t a, u, Nu, estsum1;
    arb_init(a);
    arb_init(u);
    arb_init(Nu);
    arb_init(estsum1);

    arb_ptr(Narb);
    arb_ptr(n0);

    slong expterms, taylorterms;

    Narb = ((arb_ptr *)(param))[0];
    n0 = ((arb_ptr *)(param))[1];
    expterms = ((slong *)(param))[2];
    taylorterms = ((slong *)(param))[3];

    acb_get_real(u, z);
    arb_mul(Nu, u, Narb, prec);

    est1_vecsum(estsum1, Nu, n0, expterms, taylorterms, prec);

    arb_one(a);
    arb_mul_2exp_si(a, a, -1);
    arb_neg(a, a);
    arb_pow(a, Nu, a, prec);
    arb_mul(estsum1, a, estsum1, prec);

    acb_set_arb(estint1, estsum1);

    arb_clear(a);
    arb_clear(u);
    arb_clear(Nu);
    arb_clear(estsum1);

    return 0;
}

// prepare the integrand for est1
int
f_est2_integrand(acb_ptr estint2, const acb_t z, void *param, slong order, slong prec)
{
    arb_t a, u, Nu, estsum2;
    arb_init(a);
    arb_init(u);
    arb_init(Nu);
    arb_init(estsum2);

    arb_ptr(Narb);
    arb_ptr(n0);

    slong expterms, taylorterms;

    Narb = ((arb_ptr *)(param))[0];
    n0 = ((arb_ptr *)(param))[1];
    expterms = ((slong *)(param))[2];
    taylorterms = ((slong *)(param))[3];

    acb_get_real(u, z);
    arb_mul(Nu, u, Narb, prec);

    est2_vecsum(estsum2, Nu, n0, expterms, taylorterms, prec);

    arb_one(a);
    arb_mul_2exp_si(a, a, -1);
    arb_neg(a, a);
    arb_pow(a, Nu, a, prec);
    arb_mul(estsum2, a, estsum2, prec);

    acb_set_arb(estint2, estsum2);

    arb_clear(a);
    arb_clear(u);
    arb_clear(Nu);
    arb_clear(estsum2);

    return 0;
}

//establish the required number of Taylorterms
void
est_taylorterms(arb_t taylorest, const arb_t Narb, const arb_t n0, slong expterms, slong taylorterms, slong prec)

{  
    acb_t tmp, ai, bi, intest1, intest2;
    acb_init(ai);
    acb_init(bi);
    acb_init(tmp);
    acb_init(intest1);
    acb_init(intest2);

    arb_t a, b, esta, estb, estinit1, estinit2, preterm;
    arb_init(a);
    arb_init(b);
    arb_init(esta);
    arb_init(estb);
    arb_init(estinit1);
    arb_init(estinit2);
    arb_init(preterm);

    slong goal;

    //preterm
    arb_one(a);
    arb_mul_2exp_si(a, a, -1);
    arb_pow(preterm, n0, a, prec);

    //estinit1 and estinit2
    arb_one(a);
    est1_vecsum(estinit1, a, n0, expterms, taylorterms, prec);
    est2_vecsum(estinit2, a, n0, expterms, taylorterms, prec);

    //evaluate the integrals
    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

    mag_t tol;
    mag_init(tol);

    goal = prec;
    mag_set_ui_2exp_si(tol, 1, -prec);

    arb_inv(a, Narb, prec);
    arb_zero(b);
    acb_set_arb_arb(ai, a, b);
    arb_one(a);
    acb_set_arb_arb(bi, a, b);

    void *param[4];

    param[0] = (void *) Narb;
    param[1] = (void *) n0;
    param[2] = (void *) expterms;
    param[3] = (void *) taylorterms;

    acb_calc_integrate(intest1, f_est1_integrand, param, ai, bi, goal, tol, NULL, prec);
    acb_abs(esta, intest1, prec);
    arb_mul(esta, esta, Narb, prec);
    arb_add(esta, esta, estinit1, prec);

    acb_calc_integrate(intest2, f_est2_integrand, param, ai, bi, goal, tol, NULL, prec);
    acb_abs(estb, intest2, prec);
    arb_mul(estb, estb, Narb, prec);
    arb_add(estb, estb, estinit2, prec);

    arb_add(a, esta, estb, prec);
    arb_mul(taylorest, preterm, a, prec);

    mag_clear(tol);

    acb_clear(tmp);
    acb_clear(ai);
    acb_clear(bi);
    acb_clear(intest1);
    acb_clear(intest2);

    arb_clear(a);
    arb_clear(b);
    arb_clear(esta);
    arb_clear(estb);
    arb_clear(estinit1);
    arb_clear(estinit2);
    arb_clear(preterm);

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

void* storedsums_calc(void *voidData)
{
    struct ThreadData* data=voidData;

    acb_mat_t resacb;
    arb_mat_t expvec, tayvec, resarb;

    acb_t a, b, nacb, expo, nexpo;
    acb_init(a);
    acb_init(b);
    acb_init(nacb);
    acb_init(expo);
    acb_init(nexpo);

    arb_t ab, narb, n0;
    arb_init(ab);
    arb_init(narb);
    arb_init(n0);

    slong n;

    //recover the data passed to this specific thread
    slong expterms=data->expterms;
    slong taylorterms=data->taylorterms;
    slong prec=data->prec;
    slong n00=data->n00;
    slong start=data->start;
    slong stop=data->stop;
    slong id=data->id;

    arb_set_si(n0, n00);

    //initialise the local matrices for the matet function
    acb_mat_init(resacb, expterms, taylorterms);
    arb_mat_init(expvec, expterms, 1);
    arb_mat_init(tayvec, 1, taylorterms);
    arb_mat_init(resarb, expterms, taylorterms);

    //prepare the expo power
    acb_set_arb(a, XX);
    acb_mul_onei(b, a);
    acb_sub_si(b, b, 1, prec);
    acb_mul_2exp_si(expo, b, -1);

    //produce the stored sums matrix(expterms, taylorterms)
    for (n = start; n < stop; n++)
    {
        arb_set_ui(narb, n);
        arb_add_si(narb, narb, 1, prec);
        acb_set_arb(nacb, narb);

        arb_div(ab, narb, n0, prec);
        arb_log(ab, ab, prec);

        matet(resarb, ab, expterms, taylorterms, expvec, tayvec, prec);

        acb_mat_set_arb_mat(resacb, resarb);

        acb_pow(nexpo, nacb, expo, prec);
        if (id == 0)
            acb_mat_scalar_addmul_acb(finmat1, resacb, nexpo, prec);
        if (id == 1)
            acb_mat_scalar_addmul_acb(finmat2, resacb, nexpo, prec);
        if (id == 2)
            acb_mat_scalar_addmul_acb(finmat3, resacb, nexpo, prec);
        if (id == 3)
            acb_mat_scalar_addmul_acb(finmat4, resacb, nexpo, prec);
        if (id == 4)
            acb_mat_scalar_addmul_acb(finmat5, resacb, nexpo, prec);
        if (id == 5)
            acb_mat_scalar_addmul_acb(finmat6, resacb, nexpo, prec);
        if (id == 6)
            acb_mat_scalar_addmul_acb(finmat7, resacb, nexpo, prec);
        if (id == 7)
            acb_mat_scalar_addmul_acb(finmat8, resacb, nexpo, prec);
        if (id == 8)
            acb_mat_scalar_addmul_acb(finmat9, resacb, nexpo, prec);
        if (id == 9)
            acb_mat_scalar_addmul_acb(finmat10, resacb, nexpo, prec);
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(nacb);
    acb_clear(expo);
    acb_clear(nexpo);

    arb_clear(ab);
    arb_clear(narb);
    arb_clear(n0);

    arb_mat_clear(expvec);
    arb_mat_clear(tayvec);
    acb_mat_clear(resacb);
    arb_mat_clear(resarb);

    return(NULL);
}


void
storedsums(slong res, const arb_t X, const slong N, slong numthreads, slong digits, slong prec)

{
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

    slong e, t, expterms, taylorterms, i;

    arb_init(XX);

    //shift X
    arb_one(ab);
    arb_mul_2exp_si(ab, ab, -1);  
    arb_add(XX, X, ab, prec); 

    //n0
    arb_mul_2exp_si(n0, Narb, -1);

    //establish required expterms and taylorterms
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
          est_taylorterms(res1, Narb, n0, expterms, taylorterms, prec);
    }
    expterms = expterms + 1;
    taylorterms = taylorterms + 1;

    //define output matrices per thread (defined gloabbly)
    acb_mat_t finalmat;
    acb_mat_init(finalmat, expterms, taylorterms);
    acb_mat_init(finmat1, expterms, taylorterms);
    acb_mat_init(finmat2, expterms, taylorterms);
    acb_mat_init(finmat3, expterms, taylorterms);
    acb_mat_init(finmat4, expterms, taylorterms);
    acb_mat_init(finmat5, expterms, taylorterms);
    acb_mat_init(finmat6, expterms, taylorterms);
    acb_mat_init(finmat7, expterms, taylorterms);
    acb_mat_init(finmat8, expterms, taylorterms);
    acb_mat_init(finmat9, expterms, taylorterms);
    acb_mat_init(finmat10, expterms, taylorterms);

    //prep the threads
    pthread_t thread[numthreads];

    struct ThreadData data[numthreads];

    //prep all the thread data (divide up ranges in the overall loop with start and stop of range)
    slong threadtasks = (N+numthreads-1)/numthreads;
    for (i = 0; i < numthreads; i++)
    {
    data[i].expterms=expterms;
    data[i].taylorterms=taylorterms;
    data[i].prec=prec;
    data[i].n00= N/2;
    data[i].start= i*threadtasks;
    data[i].stop= (i+1)*threadtasks;
    data[i].id= i;
    }

    data[numthreads-1].stop = N;
 
    //start the threads with an indexed (array) of a data-structure (with pointers to it from the threads).
    for (i = 0; i < numthreads; i++)
    {
    pthread_create(&thread[i], 0, storedsums_calc, &data[i]);
    }

    //wait for all threads to complete
    for (i = 0; i < numthreads; i++)
    {
    pthread_join(thread[i], NULL);
    }

    //sum all thread outputs (stored in global variables) together
    acb_mat_add(finalmat, finmat1, finmat2, prec);
    acb_mat_add(finalmat, finalmat, finmat3, prec);
    acb_mat_add(finalmat, finalmat, finmat4, prec);
    acb_mat_add(finalmat, finalmat, finmat5, prec);
    acb_mat_add(finalmat, finalmat, finmat6, prec);
    acb_mat_add(finalmat, finalmat, finmat7, prec);
    acb_mat_add(finalmat, finalmat, finmat8, prec);
    acb_mat_add(finalmat, finalmat, finmat9, prec);
    acb_mat_add(finalmat, finalmat, finmat10, prec);

    //print the stored sums for storage in file
    arb_printn(X, 30, ARB_STR_NO_RADIUS);
    printf(", %ld, %ld, %ld\n", expterms, taylorterms, digits);

    for (e = 0; e < expterms; e++)
    {
         for (t = 0; t < taylorterms; t++)
         {
              acb_get_real(re, acb_mat_entry(finalmat, e, t));
              acb_get_imag(im, acb_mat_entry(finalmat, e, t));
              arb_printn(re, digits + 10, ARB_STR_NO_RADIUS);
              printf(", ");
              arb_printn(im, digits + 10, ARB_STR_NO_RADIUS);
              if (t < taylorterms - 1) 
                  printf(", ");
              else
                  printf("\n");
         }
    }

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
    arb_clear(XX);

    acb_mat_clear(finalmat);
    acb_mat_clear(finmat1);
    acb_mat_clear(finmat2);
    acb_mat_clear(finmat3);
    acb_mat_clear(finmat4);
    acb_mat_clear(finmat5);
    acb_mat_clear(finmat6);

    acb_clear(a);
    acb_clear(b);
    acb_clear(nacb);
    acb_clear(expo);
    acb_clear(nexpo);

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

    const char *X_str, *thread_str, *digits_str;
	
    slong N, prec, numthreads, digits, res;
    res=0;

    int result = EXIT_SUCCESS;

    if (argc != 4)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    X_str = argv[1];
    thread_str = argv[2];
    digits_str = argv[3];

    numthreads = atol(thread_str);

    //precision 
    digits = atol(digits_str);
    prec = digits * 3.32192809488736 + 60;

    if (numthreads < 1 || numthreads > 10 || digits < 1)
    {
        result = EXIT_FAILURE;
        goto finish;
    }
    
    arb_set_str(X, X_str, prec);

    arb_zero(t0);
    N = get_N(t0, X, prec);

TIMEIT_ONCE_START
    storedsums(res, X, N, numthreads, digits, prec);
TIMEIT_ONCE_STOP

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Required inputs:\n");
        flint_printf("%s X numthreads digits \n\n", argv[0]);
        flint_printf(
    "This script computes the single stored sums matrix as input for 'Barrier-t-loop',\n"
    "calculations. Precision will be set by a fixed formula based on digits\n"
    "and automatic selection of the required number of Taylor expansion terms.\n" 
    "It outputs X, Expterms, Taylorterms, digits and then the finalmat sums.\n"
    "The number of threads (1 to 10) can be selected by numthreads.\n"
    "t0 has been fixed at 0.\n");
    }

    arb_clear(X);
    arb_clear(t0);
 
    flint_cleanup();

    return result;
}

