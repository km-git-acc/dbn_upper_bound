#include <stdio.h>
/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_calc.h"
#include "acb_dirichlet.h"
#include "flint/profiler.h"

//count all zeros on the critical line (for test purposes)
void
N_count(arb_t rcount, arb_t xa, slong prec)
{  
    arb_t a, b, c, pi;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(pi);

    acb_t ac, bc, cc, x;
    acb_init(ac);
    acb_init(bc);
    acb_init(cc);
    acb_init(x);

    arb_const_pi(pi, prec);

    acb_set_arb(x, xa);

    arb_printd(xa, 20);
    printf("\n");

    acb_mul_onei(ac, x);
    acb_mul_2exp_si(ac, ac, -1);
    acb_one(bc);
    acb_mul_2exp_si(bc, bc, -2);
    acb_add(ac, ac, bc, prec);
    acb_lgamma(ac, ac, prec);
    acb_get_imag(a, ac);
    arb_div(a, a, pi, prec);

    arb_mul_2exp_si(b, pi, 1);
    arb_div(b, xa, b, prec);  
    arb_log(c, pi, prec);
    arb_mul(b, b, c, prec);

    arb_sub(a, a, b, prec);

    acb_mul_onei(ac, x);
    acb_one(bc);
    acb_mul_2exp_si(bc, bc, -1);
    acb_add(bc, ac, bc, prec);
    acb_zeta(cc, bc, prec);
    acb_log(cc, cc, prec);
    acb_get_imag(c, cc);
    arb_div(c, c, pi, prec);

    arb_add(a, a, c, prec);
    arb_add_si(rcount, a, 1, prec);

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(pi);

    acb_clear(ac);
    acb_clear(bc);
    acb_clear(cc);
    acb_clear(x);
}

//f(z) = zeta'(z) / zeta(z)
int
f_zeta_frac(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_struct t[2];

    acb_init(t);
    acb_init(t + 1);

    acb_dirichlet_zeta_jet(t, z, 0, 2, prec);
    acb_div(res, t + 1, t, prec);

    acb_clear(t);
    acb_clear(t + 1);

    return 0;
}

//evaluate the contour integral
void
countour_integral(arb_t rcount, const arb_t Te, const slong prec)
{  
    arb_t ar, br;
    arb_init(ar);
    arb_init(br);

    acb_t s, t, a, b;
    acb_init(a);
    acb_init(b);
    acb_init(s);
    acb_init(t);

    mag_t tol;
    mag_init(tol);

    slong goal;

    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

    goal = prec;
    mag_set_ui_2exp_si(tol, 1, -prec);
 
    acb_zero(s);
 
    arb_set_si(ar, 100);
    arb_zero(br);
    acb_set_arb_arb(a, ar, br);
    acb_set_arb_arb(b, ar, Te);
    acb_calc_integrate(t, f_zeta_frac, NULL, a, b, goal, tol, options, prec);
    acb_add(s, s, t, prec);

    arb_one(br);
    arb_mul_2exp_si(br, br, -1);
    acb_set_arb_arb(a, ar, Te);
    acb_set_arb_arb(b, br, Te);
    acb_calc_integrate(t, f_zeta_frac, NULL, a, b, goal, tol, options, prec);
    acb_add(s, s, t, prec);
    acb_printd(t, 20);
    printf("\n");

    acb_div_onei(s, s);
    arb_zero(acb_imagref(s));
 
    acb_set_arb(t, Te);
    acb_dirichlet_hardy_theta(t, t, NULL, NULL, 1, prec);
    acb_add(s, s, t, prec);
 
    acb_const_pi(t, prec);
    acb_div(s, s, t, prec);
    acb_add_ui(s, s, 1, prec);
    acb_abs(rcount, s, prec);

    arb_clear(ar);
    arb_clear(br);

    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
    acb_clear(t);

    mag_clear(tol);
}

//test procedure for gram function abs(tn-tn1) > gram target precision. 0=false, 1 = true.
int
tntest(slong test, arb_t tn, arb_t tn1, arb_t gramprec, slong prec)
{
    arb_t a;
	arb_init(a);

    test = 0;
    arb_sub(a, tn, tn1, prec);
    arb_abs(a, a);
    if (arb_gt(a, gramprec))
        test = 1;

    arb_clear(a);

    return(test);
}

//establish good approximation of the n-th gram point (gn)
void
gram(arb_t tn1, arb_t gn, slong prec)

{  
    arb_t tn, a, b, c, d, gramprec, logtndiv2pi, pi;
    arb_init(tn);
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(gramprec);
    arb_init(logtndiv2pi);
    arb_init(pi);

    slong res; res=0;

    arb_const_pi(pi, prec);

    arb_zero(tn);
    arb_mul_2exp_si(a, gn, -1);
    arb_set_si(b, 20);
    arb_add(tn1, a, b, prec);

    arb_set_si(a, 1);
    arb_set_si(b, 10000000000000000);
    arb_div(gramprec, a, b, prec);

    while(tntest(res, tn, tn1, gramprec, prec))
    {
        arb_set(tn, tn1);
        arb_mul_2exp_si(a, pi, 1);
        arb_div(a, tn, a, prec);
        arb_log(a , a, prec);
        arb_mul_2exp_si(logtndiv2pi, a, -1);

        arb_mul(a, tn, logtndiv2pi, prec);
        arb_mul_2exp_si(b, tn, -1);
        arb_sub(a, a, b, prec);

        arb_mul_2exp_si(b, pi, -3);
        arb_sub(a, a, b, prec);

        arb_mul_si(b, tn, 48, prec);
        arb_inv(b, b, prec);
        arb_add(a, a, b, prec);

        arb_pow_ui(b, tn, 3, prec);
        arb_mul_si(b, b, 5760, prec);
        arb_set_si(c, 7);
        arb_div(b, c, b, prec);
        arb_add(a, a, b, prec);

        arb_mul(b, gn, pi, prec);
        arb_sub(a, a, b, prec);	

        arb_pow_ui(b, tn, 2, prec);
        arb_mul_si(b, b, 48, prec);
        arb_inv(b, b, prec);
        arb_sub(c, logtndiv2pi, b, prec);

        arb_pow_ui(b, tn, 4, prec);
        arb_mul_si(b, b, 1920, prec);
        arb_set_si(d, 7);
        arb_div(b, d, b, prec);
        arb_sub(c, c, b, prec);

        arb_div(a, a, c, prec);
        arb_sub(tn1, tn, a, prec);
        arb_get_mid_arb(tn, tn);
    }

    arb_clear(tn);
    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(gramprec);
    arb_clear(logtndiv2pi);
    arb_clear(pi);
}

//test procedure for gramblock function num < (gm-gn)  target precision. 0=false, 1 = true.
int
gntest(slong test, arb_t num, arb_t gn, arb_t gm, slong prec)
{
    arb_t a;
	arb_init(a);

    test = 0;
    arb_sub(a, gm, gn, prec);
    if (arb_lt(num, a))
        test = 1;

    arb_clear(a);

    return(test);
}

//test procedure to establish whether a number is even. -1=odd, 1 = even.
void
even(arb_t even, arb_t v, slong prec)
{
   arb_one(even); 
   if (arb_is_int_2exp_si(v, 1)==0)
       arb_set_si(even, -1);
}

//test procedure to establish whether a gram point is good or bad. 0=bad, 1 = good.
int
gramqual(slong test, arb_t eveninp, arb_t zinp, slong prec)
{

    arb_t a, b;
    arb_init(a);
    arb_init(b);

    acb_t s, t;
    acb_init(s);
    acb_init(t);

    test = 0;

    even(a, eveninp, prec);

    acb_set_arb(t, zinp);
    acb_dirichlet_hardy_z(s, t, NULL, NULL, 1, prec);
    acb_get_real(b, s);
    arb_mul(a, a, b, prec);
    arb_zero(b);
    if (arb_gt(a, b))
        test = 1;

    arb_clear(a);
    arb_clear(b);

    acb_clear(s);
    acb_clear(t);

    return(test);
}

//establish upper bound on S(gm)
void
adjustments(arb_t sbound, arb_t gm, slong samp, slong prec)

{  
    arb_t a, b, c, d, e, f, n, htot, step, pi, n0d128, n2d30, samparb;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(e);
    arb_init(f);
    arb_init(n);
    arb_init(htot);
    arb_init(step);
    arb_init(pi);
    arb_init(n0d128);
    arb_init(n2d30);
    arb_init(samparb);

    acb_t s, t;
    acb_init(s);
    acb_init(t);

    arf_t u;
    arf_init(u);

    slong i, k, m, rosserfail, res; 

    res=0; rosserfail = 0;

    arb_const_pi(pi, prec);

    arb_set_si(a, 1);
    arb_set_si(b, 10);
    arb_div(step, a, b, prec);

    arb_set_si(a, 128);
    arb_set_si(b, 1000);
    arb_div(n0d128, a, b, prec);

    arb_set_si(a, 230);
    arb_set_si(b, 100);
    arb_div(n2d30, a, b, prec);

    arb_set_si(samparb, samp);

    arb_zero(htot);

    arb_mat_t gmat, hmat;
    arb_mat_init(gmat, 1, samp);
    arb_mat_init(hmat, 1, samp);

    for (i = 0; i < samp; i++)
    {
       arb_add_si(a, gm, i, prec);
       gram(b, a, prec);  
       arb_set(arb_mat_entry(gmat, 0, i), b);
       arb_zero(b); 
       arb_set(arb_mat_entry(hmat, 0, i), b); 
    }

    k = 0;
    while (k < samp)
    {
        arb_zero(n);
        arb_zero(c);
        arb_add_si(a, gm, k, prec);
        arb_set(b, arb_mat_entry(gmat, 0, k));
        arb_set(d, b);
        arb_set(e, b);
        while(gramqual(res, a, d, prec) == 0 && gramqual(res, a, e, prec) == 0)
        {
           arb_mul(c, n, step, prec);
           arb_sub(d, b, c, prec);
           arb_add(e, b, c, prec);
           arb_add_si(n, n, 1, prec);

           if (arb_gt(n, samparb))
           {
              arb_add_si(f, gm, k, prec);
              arb_get_ubound_arf(u, f, prec);
              rosserfail = arf_get_d(u, ARF_RND_NEAR) - 1;
              break;
           }
        }

        if (gramqual(res, a, d, prec) > 0)
        {
           arb_neg(e, c);
           arb_set(arb_mat_entry(hmat, 0, k), e);
        }
        else
        {
           arb_set(arb_mat_entry(hmat, 0, k), c);
        }

        arb_add(htot, htot, arb_mat_entry(hmat, 0, k), prec);

        if (k > 0)
        {
            arb_add(a, arb_mat_entry(gmat, 0, k), arb_mat_entry(hmat, 0, k), prec);
            arb_add(b, arb_mat_entry(gmat, 0, k-1), arb_mat_entry(hmat, 0, k-1), prec);
            if (arb_le(a, b))
            {
                arb_zero(htot);
                k = -1;
                arb_mul_2exp_si(step, step, -1);
            }
        }
        k = k + 1; 
    }

    printf("Table headers:  | m | g_m | h_m | g_m + h_m | Z(g_m + h_m) | \n\n");
    m = 0;
    while (m < samp)
    {
        arb_add_si(a, gm, m, prec);
        arb_set(b, arb_mat_entry(gmat, 0, m));
        arb_set(c, arb_mat_entry(hmat, 0, m));
        arb_add(d, arb_mat_entry(gmat, 0, m), arb_mat_entry(hmat, 0, m), prec);
        acb_set_arb(t, d);
        acb_dirichlet_hardy_z(s, t, NULL, NULL, 1, prec);
        acb_get_real(e, s);
        arb_get_mid_arb(e, e);
        arb_get_ubound_arf(u, a, prec);
        slong a1=arf_get_si(u, ARF_RND_NEAR);
        arb_get_ubound_arf(u, b, prec);
        double b1=arf_get_d(u, ARF_RND_NEAR);
        arb_get_ubound_arf(u, c, prec);
        double c1=arf_get_d(u, ARF_RND_NEAR);
        arb_get_ubound_arf(u, d, prec);
        double d1=arf_get_d(u, ARF_RND_NEAR);
        arb_get_ubound_arf(u, e, prec);
        double e1=arf_get_d(u, ARF_RND_NEAR);
        printf("| %ld | % 8.3f | % 03.5f | % 8.3f | % 8.3f | ", a1, b1, c1, d1, e1);
        if (a1 == rosserfail)
           printf(" <- failure of Rosser's rule");
        printf("\n");

        m = m + 1;
    }

    arb_set(a, arb_mat_entry(gmat, 0, 0));
    arb_set(b, arb_mat_entry(gmat, 0, samp-1));
    arb_sub(a, b, a, prec);

    arb_mul_2exp_si(d, pi, 1);
    arb_div(e, b, d, prec);
    arb_log(e, e, prec);

    arb_mul(e, e, n0d128, prec);
    arb_add(e, e, n2d30, prec);
    arb_add(e, e, htot, prec);
    arb_div(e, e, a, prec);
    arb_add_si(sbound, e, 1, prec);

    arb_get_ubound_arf(u, sbound, prec);
    double t1a=arf_get_d(u, ARF_RND_NEAR);
    arb_get_ubound_arf(u, gm, prec);
    slong t1b=arf_get_d(u, ARF_RND_NEAR);
    arb_get_ubound_arf(u, arb_mat_entry(gmat, 0, 0), prec);
    double t1c=arf_get_d(u, ARF_RND_NEAR);
    printf("\nUpper bound of S(g_m) : %f \n\n", t1a);

    if (rosserfail == 0)
    {
       if (t1a < 2)
          printf("Number of non-trivial zeros in the strip <= Gram(%ld)=%f : %ld \n\n", t1b, t1c, t1b +1);
       else
          printf("! Turing test failed, upper bound for S(t) > 2 ! \n\n");
    }
    else
          printf("! Turing test failed, a violation of Rosser's rule occurred in Gramblock %ld : !\n\n", rosserfail);

    arb_mat_clear(gmat);
    arb_mat_clear(hmat);

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(e);
    arb_clear(f);
    arb_clear(n);
    arb_clear(htot);
    arb_clear(step);
    arb_clear(pi);
    arb_clear(n0d128);
    arb_clear(n2d30);
    arb_clear(samparb);

    acb_clear(s);
    acb_clear(t);

    arf_clear(u);
}

//establish the required accuracy for Z(t)
slong
z_accuracy(slong zpreccnt, arb_t xa, slong targetprec)
{
    acb_t a, s;
    acb_init(a);
    acb_init(s);

    arb_t ar;
    arb_init(ar);

    acb_set_arb(a, xa);

    zpreccnt = targetprec; 
    acb_dirichlet_hardy_z(s, a, NULL, NULL, 1, zpreccnt);

    while(arb_rel_accuracy_bits(acb_realref(s)) < targetprec)
    {
        zpreccnt = zpreccnt + 1;
        acb_dirichlet_hardy_z(s, a, NULL, NULL, 1, zpreccnt);
    }

    acb_clear(a);
    acb_clear(s);

    arb_clear(ar);

    return(zpreccnt);
}

int main(int argc, char *argv[])
{
    arb_t a, b, xa, res;
    arb_init(a);
    arb_init(b);
    arb_init(xa);
    arb_init(res);

    acb_t s, t;
    acb_init(s);
    acb_init(t);

    arf_t u;
    arf_init(u);

    const char *xa_str, *mode_str;
    int result = EXIT_SUCCESS;
    slong prec, ztargetprec, samp, out; out = 0;
    int mode;

    if (argc != 3)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    //precision 
    prec = 120;

    xa_str = argv[1];
    arb_set_str(xa, xa_str, prec);
    mode_str = argv[2];
    mode=atoi(mode_str);
   
TIMEIT_ONCE_START

    switch (mode)
    {
    case 1:
        gram(a, xa, prec);
        arb_get_mid_arb(a,a);
        if (gramqual(out, xa, a, prec) == 0)
        {
           printf("\n! This is a bad Gram point and the Turing test requires a good Gram point ! \n\n");
           break;
        }

        arb_log(b, a, prec); 
        arb_add_si(b, b, 10, prec);
        arb_get_ubound_arf(u, b, prec);
        samp=arf_get_si(u, ARF_RND_NEAR);
        printf("\nNumber of samples: %ld\n\n", samp);

        //prec = z_accuracy(prec, a, 40, 16);

        adjustments(res, xa, samp, prec);

        break;

     case 2:
        gram(a, xa, prec);
        arb_get_mid_arb(a,a);
        if (gramqual(out, xa, a, prec) == 0)
        {
           printf("\n! This is a bad Gram point and the N(t)-formula then does become unreliable ! \n\n");
           break;
        }

        printf("\nProcessing... \n\n");

        N_count(res, a, prec);

        arb_get_ubound_arf(u, xa, prec);
        slong t2a=arf_get_d(u, ARF_RND_NEAR);
        arb_get_ubound_arf(u, a, prec);
        double t2b=arf_get_d(u, ARF_RND_NEAR);
        arb_get_ubound_arf(u, res, prec);
        long t2c=arf_get_d(u, ARF_RND_NEAR);
        printf("\nNumber of non-trivial zeros in the strip <= Gram(%ld)=%8.8f : %ld \n\n", t2a, t2b, t2c);
        break;

     case 3:
        gram(a, xa, prec);
        arb_get_mid_arb(a, a);

        prec = 32;
        arb_set_si(b, 10000);
        if (arb_gt(xa, b))
           prec = 48;

        arb_set_si(b, 1000000000);
        if (arb_gt(xa, b))
        {
           printf("\n! Please pick a value < 10^9, otherwise the integral takes forever ! \n\n");
           break;
        }

        printf("\nProcessing... \n\n");

        countour_integral(res, a, prec);

        prec = 120;

        arb_get_ubound_arf(u, xa, prec);
        slong t3a=arf_get_d(u, ARF_RND_NEAR);
        arb_get_ubound_arf(u, a, prec);
        double t3b=arf_get_d(u, ARF_RND_NEAR);
        arb_get_ubound_arf(u, res, prec);
        long t3c=arf_get_d(u, ARF_RND_NEAR);
        printf("\nNumber of non-trivial zeros in the strip <= Gram(%ld)=%8.8f : %ld \n\n", t3a, t3b, t3c);
        break;

     case 4:
        gram(a, xa, prec);

        arb_get_ubound_arf(u, xa, prec);
        slong t4a=arf_get_d(u, ARF_RND_NEAR);
        printf("\nGram point(%ld) = ", t4a);
        arb_printn(a, 20, ARB_STR_NO_RADIUS);
        if (gramqual(out, xa, a, prec))
           printf(" and is Good.");
        else
           printf(" and is Bad.");
        printf(" (15 digits accurate) \n\n");
        break;

     case 5:
        //determine the required decimal precision (15 digits) for the start value Z(gram(gi))
        ztargetprec = 15 * 3.32192809488736 + 10;

        gram(a, xa, prec);
        prec = z_accuracy(prec, a, ztargetprec);

        acb_set_arb(t, a);
        acb_dirichlet_hardy_z(s, t, NULL, NULL, 1, prec);
        acb_get_real(a, s);

        arb_get_ubound_arf(u, a, prec);
        slong t5a=arf_get_d(u, ARF_RND_NEAR);
        printf("\n");
        printf("Z(gp_%ld) = ", t5a);
        arb_printn(a, 15, ARB_STR_NO_RADIUS);
        printf(" (15 digits accurate) \n\n");
        break;

     case 6:
        //determine the required decimal precision (15 digits) for the start value Z(gram(gi))
        ztargetprec = 15 * 3.32192809488736 + 10;

        prec = z_accuracy(prec, xa, ztargetprec);

        acb_set_arb(t, xa);
        acb_dirichlet_hardy_z(s, t, NULL, NULL, 1, prec);
        acb_get_real(a, s);

        arb_get_ubound_arf(u, xa, prec);
        double t6a=arf_get_d(u, ARF_RND_NEAR);
        printf("\n");
        printf("Z(%f) = ", t6a);
        arb_printn(a, 15, ARB_STR_NO_RADIUS);
        printf(" (15 digits accurate) \n\n");
        break;

        default:
        result = EXIT_FAILURE;
        goto finish;
    }

TIMEIT_ONCE_STOP

finish:
    if (result == EXIT_FAILURE)
    {
        flint_printf("\nRequired inputs:\n");
        flint_printf("%s m mode \n\n", argv[0]);
        flint_printf(
    "This script implements three ways to count zeros of Z(t) < gp_m in the critical strip:\n\n"
    "Mode = 1: Turing's method testing S(gp_m) < 2 using Log(gp_m)+10 samples.\n"
    "Mode = 2: N(gp_m)-formula (unreliable at Bad gram points).\n"
    "Mode = 3: Countour integration (only works for gp_m < 10^9).\n"
    "Mode = 4: Util 1: Return value of gp_m and whether it is a Good or Bad Gram point.\n"
    "Mode = 5: Util 2: Return value of Z(gp_m) at 15 digits accuracy.\n"
    "Mode = 6: Util 3: Return value of Z(m) at 15 digits accuracy.\n\n");
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(xa);
    arb_clear(res);

    acb_clear(s);
    acb_clear(t);

    arf_clear(u);
 
    flint_cleanup();

    return result;
}
