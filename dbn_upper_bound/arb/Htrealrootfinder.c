/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
#include "acb_mat.h"
#include "acb_calc.h"
#include "acb_dirichlet.h"
#include "flint/profiler.h"
#include "pthread.h"

arb_t t;
arb_mat_t Threadwork;

acb_t z;
acb_mat_t Rsums;

struct ThreadData {
    slong id, prec;
}; 

// Prepare the 'integrand' ofr the Riemann summation
void
xi_integrand(acb_t result, const arb_t u, slong prec)
{
    arb_t a;
    arb_init(a);

    acb_t ac, bc, xi;
    acb_init(ac);
    acb_init(bc);
    acb_init(xi);

    //(1+z*I)/2
    acb_mul_onei(ac, z);
    acb_add_si(ac, ac, 1, prec);
    acb_mul_2exp_si(ac, ac, -1);

    //sqrt(t)*u (t can be negative, hence acb)
    acb_set_arb(bc, t);
    acb_sqrt(bc, bc, prec);	
    acb_mul_arb(bc, bc, u, prec);

    acb_add(ac, ac, bc, prec);

    //xi((1+z*I)/2 + sqrt(t)*u)
    acb_dirichlet_xi(xi, ac, prec);

    //exp(-u^2) (could be 'memoized' later)
    arb_mul(a, u, u, prec);
    arb_neg(a, a);
    arb_exp(a, a, prec);

    acb_mul_arb(result, xi, a, prec);

    arb_clear(a);

    acb_clear(ac);
    acb_clear(bc);
    acb_clear(xi);
}

//Evaluate the Riemann sum (thread).
void* Ht_Riemann_sum_thread(void *voidData)
{
    //recover the data passed to this specific thread
    struct ThreadData* data=voidData;
    slong id=data->id;
    slong prec=data->prec;

    arb_t a, i, start, stop, step;
    arb_init(a);
    arb_init(i);
    arb_init(start);
    arb_init(stop);
    arb_init(step);

    acb_t res, resd, rsum, rsumd;
    acb_init(res);
    acb_init(resd);
    acb_init(rsum);
    acb_init(rsumd);

    arb_set(start, arb_mat_entry(Threadwork, id, 0));
    arb_set(stop, arb_mat_entry(Threadwork, id, 1));
    arb_set(step, arb_mat_entry(Threadwork, id, 2));

    acb_zero(rsum); acb_zero(rsumd); arb_set(i, start);

    while(arb_lt(i, stop))
    {
        arb_get_mid_arb(i, i);
        xi_integrand(res, i, prec);
        acb_add(rsum, rsum, res, prec);
        acb_mul_arb(resd, res, i, prec);
        acb_add(rsumd, rsumd, resd, prec);
        arb_add(i, i, step, prec);
    }

    acb_set(acb_mat_entry(Rsums, id, 0), rsum);
    acb_set(acb_mat_entry(Rsums, id, 1), rsumd);

    arb_clear(a);
    arb_clear(i);
    arb_clear(start);
    arb_clear(stop);
    arb_clear(step);

    acb_clear(res);
    acb_clear(resd);
    acb_clear(rsum);
    acb_clear(rsumd);

    flint_cleanup();

    return(NULL);
}

//Evaluate the Riemann sum. The loop is split up in line with the number of threads chosen.
void
Ht_Riemann_sum(acb_t rs, arb_t step, slong numthreads, slong prec)
{
    arb_t a;
    arb_init(a);

    acb_t ac;
    acb_init(ac);
	
    slong i;
    
    //prep the threads
    pthread_t thread[numthreads];

    struct ThreadData data[numthreads];

    //prep all the thread data
    for (i = 0; i < numthreads; i++)
    {
        data[i].id = i;
        data[i].prec = prec;
    }

    //start the threads with an indexed (array) of a data-structure (with pointers to it from the threads).
    for (i = 0; i < numthreads; i++)
    {
        pthread_create(&thread[i], 0, Ht_Riemann_sum_thread, &data[i]);
    }

    //wait for all threads to complete
    for (i = 0; i < numthreads; i++)
    {
        pthread_join(thread[i], NULL);
    }

    acb_zero(rs); //acb_zero(rs + 1);
    for (i = 0; i < numthreads; i++)
    {
        acb_add(rs, rs, acb_mat_entry(Rsums, i, 0), prec);
    }

    //rs/(8*sqrt(pi)*step)
    acb_mul_2exp_si(rs, rs, -3);
    arb_const_pi(a, prec);
    arb_sqrt(a, a, prec);
    acb_div_arb(rs, rs, a, prec);
    acb_mul_arb(rs, rs, step, prec);

    arb_clear(a);

    acb_clear(ac);
}

//Find an exact root at 5 decimal places of Ht(t) in the interval begm to endm
void
FindExactRoot(arb_t res, arb_t begm, arb_t endm, arb_t step, slong numthreads, slong prec)
{
    arb_t a, b, c, d, e, diffba, accrt;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(e);
    arb_init(diffba);
    arb_init(accrt);

    acb_struct outcome[2];
    acb_init(outcome);
    acb_init(outcome + 1);

    arb_set_si(a, 100000);
    arb_inv(accrt, a, prec);

    arb_set(a, begm);
    arb_set(b, endm);
    arb_set(res, a);
    arb_sub(diffba, b, a, prec);

    while(arb_ge(diffba, accrt))
    {
        arb_add(c, a, b, prec);
        arb_mul_2exp_si(res, c, -1);

        arb_zero(e);
        acb_set_arb_arb(z, res, e);
        Ht_Riemann_sum(outcome, step, numthreads, prec);
        acb_get_real(c, outcome);

        acb_set_arb_arb(z, a, e);
        Ht_Riemann_sum(outcome, step, numthreads, prec);
        acb_get_real(d, outcome);

        arb_mul(d, c, d, prec);
        if (arb_is_negative(d))
           arb_set(b, res);
        else
           arb_set(a, res);

        arb_sub(diffba, b, a, prec);
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(e);
    arb_clear(diffba);
    arb_clear(accrt);

    acb_clear(outcome);
    acb_clear(outcome + 1);
}

//Allocate the work across threads (stored in global matrix Threadwork)
void
Allocate_worktothreads(slong nix, arb_t intlim, arb_t step, slong numthreads, slong prec)
{
    arb_t a, workload;
    arb_init(a);
    arb_init(workload);

    slong i;

    arb_mul_2exp_si(workload, intlim, 1);
    arb_div(workload, workload, step, prec);
    arb_add_si(workload, workload, 1, prec);
    arb_div_si(workload, workload, numthreads, prec);
    arb_floor(workload, workload, prec);

    for(i = 0; i < numthreads; i++)
    {
        arb_mul_si(a, workload, i, prec);
        arb_mul(a, a, step, prec);
        arb_sub(a, a, intlim, prec);
        arb_set(arb_mat_entry(Threadwork, i, 0), a);

        arb_mul_si(a, workload, i + 1, prec);
        arb_mul(a, a, step, prec);
        arb_sub(a, a, intlim, prec);
        arb_set(arb_mat_entry(Threadwork, i, 1), a);

        arb_set(arb_mat_entry(Threadwork, i, 2), step);
    }
    arb_set(arb_mat_entry(Threadwork, numthreads - 1, 1), intlim);

    arb_clear(a);
    arb_clear(workload);
}

//Establish the required 'stepsize' to achieve 6 digits accuracy of the Riemann sum
void
Establish_step(arb_t step, arb_t intlim, arb_t stepprev, arb_t acc, slong numthreads, slong prec)
{
    arb_t prevs, currs, diffs;
    arb_init(prevs);
    arb_init(currs);
    arb_init(diffs);

    acb_t s;
    acb_init(s);

    slong res; res = 0;

    arb_set(step, stepprev);

    arb_one(diffs);

    Allocate_worktothreads(res, intlim, step, numthreads, prec);
    Ht_Riemann_sum(s, step, numthreads, prec);
    acb_get_real(prevs, s);

    while (arb_gt(diffs, acc))
    {
        arb_mul_2exp_si(step, step, -1);

        Allocate_worktothreads(res, intlim, step, numthreads, prec);
        Ht_Riemann_sum(s, step, numthreads, prec);
        acb_get_real(currs, s);

        arb_div(diffs, currs, prevs, prec);
        arb_sub_si(diffs, diffs, 1, prec);
        arb_abs(diffs, diffs);
        arb_get_mid_arb(diffs, diffs);

        arb_set(prevs, currs);
    }

    arb_clear(prevs);
    arb_clear(currs);
    arb_clear(diffs);

    acb_clear(s);
}


//Establish the required 'intlim' to achieve 6 digits accuracy of the Riemann sum
void
Establish_intlim(arb_t intlim, arb_t step, arb_t intlimprev, arb_t acc, slong numthreads, slong prec)
{
    arb_t prevs, currs, diffs;
    arb_init(prevs);
    arb_init(currs);
    arb_init(diffs);

    acb_t s;
    acb_init(s);

    slong res; res = 0;

    arb_set(intlim, intlimprev);

    arb_one(diffs);

    Allocate_worktothreads(res, intlim, step, numthreads, prec);
    Ht_Riemann_sum(s, step, numthreads, prec);
    acb_get_real(prevs, s);

    while (arb_gt(diffs, acc))
    {
        arb_add_si(intlim, intlim, 1, prec);

        Allocate_worktothreads(res, intlim, step, numthreads, prec);
        Ht_Riemann_sum(s, step, numthreads, prec);
        acb_get_real(currs, s);

        arb_div(diffs, currs, prevs, prec);
        arb_sub_si(diffs, diffs, 1, prec);
        arb_abs(diffs, diffs);
        arb_get_mid_arb(diffs, diffs);

        arb_set(prevs, currs);
    }

    arb_clear(prevs);
    arb_clear(currs);
    arb_clear(diffs);

    acb_clear(s);
}

//Main t and x-loops to find the real roots in the rectangle xs..xe (in xsteps), ts..te (in tsteps) 
void
Find_roots_loops(slong res, arb_t xs, arb_t xe, arb_t xstep, arb_t ts, arb_t te, arb_t tstep, slong numthreads, slong prec)
{
    arb_t a, acc, y, xcnt, tcnt, currz, prevz, step, stepstrt, intlim, intlimprev, root;
    arb_init(a);
    arb_init(acc);
    arb_init(y);
    arb_init(xcnt);
    arb_init(tcnt);
    arb_init(currz);
    arb_init(prevz);
    arb_init(step);
    arb_init(stepstrt);
    arb_init(intlim);
    arb_init(intlimprev);
    arb_init(root);

    acb_t s;
    acb_init(s);

    arf_t u;
    arf_init(u);

    slong rest; rest = 0;
    double tmp;

    arb_zero(y);

    arb_set_d(acc, 0.000001);

    //clean the radius of the variables to ensure proper comparison occurs in while-loops
    arb_get_mid_arb(xs, xs);
    arb_get_mid_arb(xe, xe);
    arb_get_mid_arb(ts, ts);
    arb_get_mid_arb(te, te);

    arb_set(tcnt, ts);
    while(arb_le(tcnt, te))
    {
        arb_set(t, tcnt);

        //for each incremental t step, automatically establish the right intlim and sum stepsize at the highest x (xe)
        acb_set_arb_arb(z, xe, y);
        arb_set_d(intlimprev, 3); arb_one(stepstrt);

        Establish_intlim(intlim, stepstrt, intlimprev, acc, numthreads, prec);

        Establish_step(step, intlim, stepstrt, acc, numthreads, prec);

        arb_set(xcnt, xs);
        acb_set_arb_arb(z, xcnt, y);
        Ht_Riemann_sum(s, step, numthreads, prec);
        acb_get_real(prevz, s);
        while(arb_lt(xcnt, xe))
        {
            arb_add(xcnt, xcnt, xstep, prec);
            arb_get_mid_arb(xcnt, xcnt);

            acb_set_arb_arb(z, xcnt, y);
            Ht_Riemann_sum(s, step, numthreads, prec);
            acb_get_real(currz, s);

            arb_mul(a, currz, prevz, prec);
            if (arb_is_negative(a))
            {
               arb_get_ubound_arf(u, t, prec);
               tmp=arf_get_d(u, ARF_RND_NEAR);
               printf("%f, ", tmp);
               arb_sub(a, xcnt, xstep, prec);
               FindExactRoot(root, a, xcnt, step, numthreads, prec);
               arb_get_ubound_arf(u, root, prec);
               tmp=arf_get_d(u, ARF_RND_NEAR);
               printf("%3.6f\n", tmp);
            }

            arb_set(prevz, currz);
        }

        arb_add(tcnt, tcnt, tstep, prec);
        arb_get_mid_arb(tcnt, tcnt);
    }

    arb_clear(a);
    arb_clear(acc);
    arb_clear(y);
    arb_clear(xcnt);
    arb_clear(tcnt);
    arb_clear(currz);
    arb_clear(prevz);
    arb_clear(step);
    arb_clear(stepstrt);
    arb_clear(intlim);
    arb_clear(intlimprev);
    arb_clear(root);

    acb_clear(s);

    arf_clear(u);
}

int main(int argc, char *argv[])
{
    arb_t a, xs, xe, xstmp, xetmp, xwork, xstep, ts, te, tstep;
    arb_init(a);
    arb_init(xs);
    arb_init(xe);
    arb_init(xstmp);
    arb_init(xetmp);
    arb_init(xwork);
    arb_init(xstep);
    arb_init(ts);
    arb_init(te);
    arb_init(tstep);
    arb_init(t);

    acb_init(z);

    const char *xs_str, *xe_str, *xstep_str, *ts_str, *te_str, *tstep_str, *numclients_str, *clientid_str, *numthreads_str;
	
    slong numclients, clientid, numthreads, prec, res; res = 0;

    int result = EXIT_SUCCESS;

    if (argc != 10)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    xs_str = argv[1];
    xe_str = argv[2];
    xstep_str = argv[3];
    ts_str = argv[4];
    te_str = argv[5];
    tstep_str = argv[6];
    numclients_str = argv[7];
    clientid_str = argv[8];
    numthreads_str = argv[9];

    //Working precision, sufficient for x < 10^14 and t > - 10^5 
    prec = 20 * 3.32192809488736 + 10;

    arb_set_str(xs, xs_str, prec);
    arb_set_str(xe, xe_str, prec);
    arb_set_str(xstep, xstep_str, prec);
    arb_set_str(ts, ts_str, prec);
    arb_set_str(te, te_str, prec);
    arb_set_str(tstep, tstep_str, prec);
    numclients = atol(numclients_str);
    clientid = atol(clientid_str);

    numthreads = atol(numthreads_str);

    if(clientid < 0 || clientid > numclients - 1 || numclients < 1 || numthreads < 1 || numthreads > 254)
    {
        result = EXIT_FAILURE; 
        goto finish;
    }

    if(arb_lt(xe, xs) || arb_lt(te, ts) || arb_is_nonpositive(xstep) || arb_is_nonpositive(tstep))
    {
        result = EXIT_FAILURE; 
        goto finish;
    }

    arb_mat_init(Threadwork, numthreads, 3);
    acb_mat_init(Rsums, numthreads, 2);

TIMEIT_ONCE_START

    //Establish the x-range workload for the selected client.
    arb_sub(xwork, xe, xs, prec);
    arb_div_si(xwork, xwork, numclients, prec);
    arb_floor(xwork, xwork, prec);
    //Set xs
    arb_mul_si(xstmp, xwork, clientid, prec);
    arb_add(xs, xs, xstmp, prec);
    //Set xe
    if (clientid < numclients - 1)
        arb_add(xe, xs, xwork, prec);

    Find_roots_loops(res, xs, xe, xstep, ts, te, tstep, numthreads, prec);

TIMEIT_ONCE_STOP

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Required inputs:\n");
        flint_printf("%s xs, xe, xstep, ts, te, tstep, numclients, clientid, numthreads \n\n", argv[0]);
        flint_printf(
    "This script aims to compute all real zeros of H_t(x)\n"
    "in the rectangle xs <= x <= xe, ts <= t <= te, given\n"
    "the 'stepsizes' xstep and tstep that should be > 0.\n"
    "The work can be parallellized through grid computing by chosing\n" 
    "a number of clients (1..N) and a clientid (0..N-1) to split the x-range.\n"
    "When numclients = 1 and clientid = 0 the xrange is not split.\n"
    "With numthreads the number of threads (< 250) can be selected.\n"
    "The accuracy of the zeros has been fixed at 5 decimal places.\n"
    "NOTE: this script has been designed to explore the behavior of the zeros\n"
    "of Ht at larger negative t. It does not intend to find all real zeros,\n"
    "although it will when xstep is made arbitrary small.\n");
    }

    arb_clear(a);
    arb_clear(xs);
    arb_clear(xe);
    arb_clear(xstmp);
    arb_clear(xetmp);
    arb_clear(xwork);
    arb_clear(xstep);
    arb_clear(ts);
    arb_clear(te);
    arb_clear(tstep);

    arb_mat_clear(Threadwork);
    acb_mat_clear(Rsums);
 
    flint_cleanup();

    return result;
}