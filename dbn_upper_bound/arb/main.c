/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
#include "acb_poly.h"
 
 
typedef acb_poly_struct *acb_poly_ptr;
typedef const acb_poly_struct *acb_poly_srcptr;
 
void
acb_poly_re_inplace(acb_poly_t p)
{
    slong i;
    for (i = 0; i < acb_poly_length(p); i++)
    {
        acb_ptr z = acb_poly_get_coeff_ptr(p, i);
        arb_zero(acb_imagref(z));
    }
}
 
void
acb_poly_re(acb_poly_t p, const acb_poly_t q)
{
    if (p != q)
        acb_poly_set(p, q);
 
    acb_poly_re_inplace(p);
}
 
void
acb_poly_im_inplace(acb_poly_t p)
{
    slong i;
    for (i = 0; i < acb_poly_length(p); i++)
    {
        acb_ptr z = acb_poly_get_coeff_ptr(p, i);
        arb_set(acb_realref(z), acb_imagref(z));
        arb_zero(acb_imagref(z));
    }
}
 
void
acb_poly_im(acb_poly_t p, const acb_poly_t q)
{
    if (p != q)
        acb_poly_set(p, q);
 
    acb_poly_im_inplace(p);
}

//procedure to calculate the absolute value of an acb_poly function
void
acb_poly_abs_series(acb_poly_t res,
        const acb_poly_t p, slong n, slong prec)
{
    acb_poly_t a, b;
 
    acb_poly_init(a);
    acb_poly_init(b);
 
    acb_poly_re(a, p);
    acb_poly_mullow(a, a, a, n, prec);
 
    acb_poly_im(b, p);
    acb_poly_mullow(b, b, b, n, prec);
 
    acb_poly_add_series(res, a, b, n, prec);
 
    acb_poly_sqrt_series(res, res, n, prec);
    acb_poly_re_inplace(res);
 
    acb_poly_clear(a);
    acb_poly_clear(b);
}


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

//procedure to derive the x-value from the N-value
void
_acb_poly_xNpoly_series(acb_poly_t res,
        const acb_poly_t N, const acb_poly_t t, slong n, slong prec)
{
    acb_poly_t p, q;
    acb_t pi;
 
    acb_init(pi);
    acb_const_pi(pi, prec);
 
    acb_poly_init(p);
    acb_poly_mullow(p, N, N, n, prec);
    acb_poly_scalar_mul(p, p, pi, prec);
    acb_poly_scalar_mul_2exp_si(p, p, 2);
 
    acb_poly_init(q);
    acb_poly_scalar_mul(q, t, pi, prec);
    acb_poly_scalar_mul_2exp_si(q, q, -2);
 
    acb_poly_sub_series(res, p, q, n, prec);
 
    acb_clear(pi);
    acb_poly_clear(p);
    acb_poly_clear(q);
}

//procedure to calculate the modgamma-value
static void
_acb_poly_modgamma_series(acb_poly_t res,
        const acb_poly_t xN, 
        const acb_poly_t y, slong n, slong prec)
{
    acb_poly_t a, b, c;
	acb_t pi,d;
 
    acb_poly_init(a);
    acb_poly_init(b);
    acb_poly_init(c);
	
    acb_init(pi);
    acb_const_pi(pi, prec);
	acb_init(d);
	acb_set_ui(d, 20);
    acb_div_ui(d, d, 1000, prec);

    acb_poly_scalar_mul(a, y, d, prec);
    acb_poly_exp_series(a, a, n, prec);

    acb_poly_scalar_div(b, xN, pi, prec);
    acb_poly_scalar_mul_2exp_si(b, b, -2);

    acb_poly_neg(c, y);
    acb_poly_scalar_mul_2exp_si(c, c, -1);
    acb_poly_pow_series(c,b,c,n,prec);

    acb_poly_mullow(res, a, c, n, prec);
 
    acb_poly_clear(a);
    acb_poly_clear(b);
    acb_poly_clear(c);
    acb_clear(d);
    acb_clear(pi);
}
 
//procedure to calculate the sig exponent-value 
static void
_acb_poly_sig_series(acb_poly_t res,
        const acb_poly_t xN, const acb_poly_t t, const acb_poly_t y,
        slong n, slong prec)
{
    acb_poly_t a,b,c,d,e,f,g,z;

    acb_t pi, three;
    acb_poly_init(a);
    acb_poly_init(b);
    acb_poly_init(c);
    acb_poly_init(d);
    acb_poly_init(e);
    acb_poly_init(f);
    acb_poly_init(g);
    acb_poly_init(z);
    acb_init(three);
    acb_set_ui(three, 3);
    acb_poly_zero(z);
    acb_init(pi);
    acb_const_pi(pi, prec);
 
    acb_poly_add_si(a, y, 1, prec);
    acb_poly_scalar_mul_2exp_si(a, a, -1);

    acb_poly_scalar_div(b, xN, pi, prec);
    acb_poly_scalar_mul_2exp_si(b, b, -2);
    acb_poly_log_series(b, b, n, prec);
    acb_poly_mullow(b, b, t, n, prec);
    acb_poly_scalar_mul_2exp_si(b, b, -2);
    acb_poly_add_series(res, a, b, n, prec);

    acb_poly_mullow(c, xN, xN, n, prec);
    acb_poly_scalar_mul_2exp_si(c, c, 1);
    acb_poly_div_series(c, t, c, n, prec);

    acb_poly_scalar_mul_2exp_si(d, y, 2);
    acb_poly_add_si(e, y, 1, prec);
    acb_poly_mullow(e, e, d, n, prec);
    acb_poly_mullow(f, xN, xN, n, prec);
    acb_poly_div_series(f, e, f, n, prec);

    acb_poly_neg(g, y);
	acb_poly_scalar_mul(g, g, three, prec);
    acb_poly_add_si(g, g, 1, prec);
    acb_poly_add_series(g, g, f, n, prec);

    acb_poly_max_series(g, g, z, n, prec);

    acb_poly_mullow(g, c, g, n, prec);
    acb_poly_sub_series(res, res, g, n, prec);

    acb_poly_clear(a);
    acb_poly_clear(b);
    acb_poly_clear(c);
    acb_poly_clear(d);
    acb_poly_clear(e);
    acb_poly_clear(f);
    acb_poly_clear(g);
    acb_poly_clear(z);
    acb_clear(three);
    acb_clear(pi);

}

//procedure to calculate the bt-value 
static void
_acb_poly_bt_series(acb_poly_t res,
        const acb_t logk, const acb_poly_t t, slong n, slong prec)
{
    acb_poly_scalar_mul(res, t, logk, prec);
    acb_poly_scalar_mul(res, res, logk, prec);
    acb_poly_scalar_mul_2exp_si(res, res, -2);
    acb_poly_exp_series(res, res, n, prec);
}

//procedure to calculate the at-value  
static void
_acb_poly_at_series(acb_poly_t res,
        const acb_t logk, const acb_poly_t t, const acb_poly_t y,
        slong n, slong prec)
{
    acb_poly_t p;
 
    acb_poly_init(p);
    acb_poly_scalar_mul(p, t, logk, prec);
    acb_poly_scalar_mul_2exp_si(p, p, -2);
    acb_poly_add_series(p, p, y, n, prec);
 
    acb_poly_scalar_mul(res, p, logk, prec);
    acb_poly_exp_series(res, res, n, prec);
 
    acb_poly_clear(p);
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

//main procedure to establish the Lemma-bound  
static void
_acb_poly_abbeff_lbound(acb_poly_t res,
        const acb_poly_t t, const acb_poly_t y, const acb_poly_t Npoly,
        acb_poly_srcptr at, acb_poly_srcptr bt, acb_poly_srcptr moll,
        slong nprimes, const ulong *primes, const slong *d,
        slong D, slong N, slong n, slong prec)
{
    slong divs;
    acb_t logk, lognd;
    acb_poly_t ksummand_a, ndsum_a, absndsum_a, ndsummand_a, ndsummand, a, b;
    acb_poly_t ksummand_b, ndsum_b, absndsum_b, ndsummand_b, common, modK, sumnd, sum2;
    acb_poly_t lbound;
    acb_poly_t sig, modgam, xN, kterm, kterm1, ndterm, modmol;
    slong j, k, nd;
 
    divs = 1 << nprimes;
 
    acb_init(logk);
    acb_init(lognd);
 
    acb_poly_init(a);
    acb_poly_init(b);
    acb_poly_init(common);
    acb_poly_init(modK);
    acb_poly_init(sumnd);
    acb_poly_init(sum2);
    acb_poly_init(ksummand_a);
    acb_poly_init(ndsum_a);
    acb_poly_init(absndsum_a);
    acb_poly_init(ndsummand_a);
 
    acb_poly_init(ksummand_b);
    acb_poly_init(ndsum_b);
    acb_poly_init(absndsum_b);
    acb_poly_init(ndsummand_b);
    acb_poly_init(ndsummand);
 
    acb_poly_init(lbound);
    acb_poly_init(sig);
    acb_poly_init(modgam);
    acb_poly_init(xN);
    acb_poly_init(kterm);
    acb_poly_init(kterm1);
    acb_poly_init(ndterm);
    acb_poly_init(modmol);
 
    _acb_poly_xNpoly_series(xN, Npoly, t, n, prec);
    _acb_poly_sig_series(sig, xN, t, y, n, prec);
    _acb_poly_modgamma_series(modgam, xN, y, n, prec);

    //define common
    acb_poly_neg(lbound,modgam);
    acb_poly_add_si(lbound, lbound, 1, prec);
    acb_poly_add_si(b, modgam, 1, prec);
    acb_poly_div_series(common, lbound, b, n, prec);

    //main loop from n=1 to D*N
    for (k = 2; k <= D*N; k++)
    {
        acb_poly_zero(ndsum_a);
        acb_poly_zero(ndsum_b);
        //sub loop from nd=1 to ndivs
        for (nd = 1; nd <= divs; nd++)
        {
            if ((k % d[nd] == 0) && (k <= d[nd]*N))
            {
                j = k / d[nd];
 
                acb_poly_mullow(ndsummand_b, moll + nd, bt + j, n, prec);
                acb_poly_add_series(ndsum_b, ndsum_b, ndsummand_b, n, prec);
 
                acb_poly_mullow(ndsummand_a, moll + nd, at + j, n, prec);
                acb_poly_mullow(ndsummand_a, ndsummand_a, modgam, n, prec);
                acb_poly_add_series(ndsum_a, ndsum_a, ndsummand_a, n, prec);
            }
        }
        if (!acb_poly_is_zero(ndsum_b) || !acb_poly_is_zero(ndsum_a))
        {
            acb_set_si(logk, k);
            acb_log(logk, logk, prec);
            acb_poly_scalar_mul(kterm, sig, logk, prec);
            acb_poly_neg(kterm, kterm);
            acb_poly_exp_series(kterm, kterm, n, prec);

            acb_poly_sub_series(a, ndsum_b, ndsum_a, n, prec);
            acb_poly_abs_series(a, a, n, prec);
            acb_poly_add_series(b, ndsum_b, ndsum_a, n, prec);
            acb_poly_abs_series(b, b, n, prec);
            acb_poly_mullow(b, b, common, n, prec);
            acb_poly_max_series(sumnd, a, b, n, prec);

            acb_poly_mullow(sumnd, kterm, sumnd, n, prec);
            acb_poly_sub_series(lbound, lbound, sumnd, n, prec);
        }
    }

    //define modmol
    acb_poly_zero(modmol);
    for (nd = 1; nd <= divs; nd++)
    {
        acb_set_si(lognd, d[nd]);
        acb_log(lognd, lognd, prec);
        acb_poly_scalar_mul(ndterm, sig, lognd, prec);
        acb_poly_neg(ndterm, ndterm);
        acb_poly_exp_series(ndterm, ndterm, n, prec);
 
        acb_poly_mullow(ndsummand, moll + nd, ndterm, n, prec);
        acb_poly_abs_series(ndsummand, ndsummand, n, prec);
        acb_poly_add_series(modmol, modmol, ndsummand, n, prec);
    }

    //define modK
    acb_poly_add_si(a, xN, -6, prec);
    acb_poly_scalar_mul_2exp_si(a, a, 1);
    acb_poly_div_series(a, y, a, n, prec);
    acb_poly_mullow(modK, t, a, n, prec);

    //define sum2
    acb_poly_zero(sum2);
    for (k = 1; k <= N; k++)
    {
        acb_set_si(logk, k);
        acb_log(logk, logk, prec);
        acb_poly_scalar_mul(kterm, modK, logk, prec);
        acb_poly_exp_series(kterm, kterm, n, prec);
        acb_poly_add_si(kterm, kterm, -1, prec);
        acb_poly_sub_series(b, y, sig, n, prec);
        acb_poly_scalar_mul(kterm1, b, logk, prec);
        acb_poly_exp_series(kterm1, kterm1, n, prec);
        acb_poly_mullow(b, kterm1, kterm, n, prec);
        acb_poly_mullow(b, bt+k, b, n, prec);
        acb_poly_add_series(sum2, sum2, b, n, prec);
    }

    //perform final calculations to complete lbound
    acb_poly_mullow(a, modmol, modgam, n, prec);
    acb_poly_mullow(a, sum2, a, n, prec);
    acb_poly_sub_series(res, lbound, a, n, prec);
    
    acb_clear(logk);
    acb_clear(lognd);
 
    acb_poly_clear(ksummand_a);
    acb_poly_clear(ndsum_a);
    acb_poly_clear(absndsum_a);
    acb_poly_clear(ndsummand_a);
 
    acb_poly_clear(ksummand_b);
    acb_poly_clear(ndsum_b);
    acb_poly_clear(absndsum_b);
    acb_poly_clear(ndsummand_b);
    acb_poly_clear(ndsummand);
 
    acb_poly_clear(lbound);
    acb_poly_clear(xN);
    acb_poly_clear(kterm);
    acb_poly_clear(kterm1);
    acb_poly_clear(sumnd);
    acb_poly_clear(sum2);
    acb_poly_clear(modmol);
    acb_poly_clear(ndterm);
    acb_poly_clear(a);
    acb_poly_clear(b);
    acb_poly_clear(common);
    acb_poly_clear(modK);
}

//initialise and fill arrays for bt, at, moll  
slong run(const acb_t t, const acb_t y,
        slong Na, slong Nb, slong Nstep,
        slong nprimes, slong n, slong digits, slong prec)
{
    slong D, N, Nmax;
    acb_poly_ptr bt, at;
    slong i, j, k;
    slong high_prec;
    acb_t logk, coeff;
    acb_poly_t tpoly, ypoly, rpoly, rdpoly, Npoly;
    acb_poly_ptr moll;
    slong divs;
    ulong *primes;
    slong *d;
    slong completed_steps = 0;
 
    high_prec = digits * 3.32192809488736 + 10;
    Nmax = Nb;
 
    acb_init(logk);
    acb_init(coeff);
    acb_poly_init(tpoly);
    acb_poly_init(ypoly);
    acb_poly_init(rpoly);
    acb_poly_init(rdpoly);
    acb_poly_init(Npoly);
 
    acb_poly_set_acb(tpoly, t);
    acb_poly_set_acb(ypoly, y);
 
    divs = 1 << nprimes;
    primes = flint_malloc(nprimes * sizeof(*primes));
    _first_n_primes(primes, nprimes);
 
    D = 1;
    for (i = 0; i < nprimes; i++)
    {
        D *= primes[i];
    }
 
    bt = flint_malloc((Nmax * D + 1) * sizeof(*bt));
    at = flint_malloc((Nmax * D + 1) * sizeof(*at));
 
    for (k = 1; k <= Nmax * D; k++)
    {
        acb_set_si(logk, k);
        acb_log(logk, logk, prec);
 
        acb_poly_init(bt + k);
        acb_poly_init(at + k);
 
        _acb_poly_bt_series(bt + k, logk, tpoly, n, prec);
        _acb_poly_at_series(at + k, logk, tpoly, ypoly, n, prec);
    }
 
    d = flint_malloc((divs + 1) * sizeof(*d));
    moll = flint_malloc((divs + 1) * sizeof(*moll));
    for (i = 0; i < divs; i++)
    {
        acb_poly_init(moll + i + 1);
        acb_poly_one(moll + i + 1);
        d[i + 1] = 1;
        for (j = 0; j < nprimes; j++)
        {
            if (i & (1 << j))
            {
                d[i + 1] *= primes[j];
                acb_poly_mullow(
                        moll + i + 1, moll + i + 1, bt + primes[j], n, prec);
                acb_poly_neg(moll + i + 1, moll + i + 1);
            }
        }
    }
	
    //main loop from N-start to N-end 
    for (N = Na; N <= Nb; N += Nstep)
    {
        acb_poly_one(Npoly);
        acb_poly_shift_left(Npoly, Npoly, 1);
        acb_poly_add_si(Npoly, Npoly, N, prec);
 
        _acb_poly_abbeff_lbound(rpoly, tpoly, ypoly, Npoly,
                at, bt, moll, nprimes, primes, d,
                D, N, n, prec);
 
        //print output(s)
        acb_poly_set(rdpoly, rpoly);
        for (i = 0; i < n + 1; i++)
        {
            acb_poly_get_coeff_acb(coeff, rdpoly, 0);
            if (acb_rel_accuracy_bits(coeff) < high_prec)
            {
                goto finish;
            }
            acb_poly_derivative(rdpoly, rdpoly, prec);
        }
 
        flint_printf("%ld", N);
        acb_poly_set(rdpoly, rpoly);
        for (i = 0; i < n; i++)
        {
            acb_poly_get_coeff_acb(coeff, rdpoly, 0);
            flint_printf("\t");
            arf_printd(arb_midref(acb_realref(coeff)), digits);
            acb_poly_derivative(rdpoly, rdpoly, prec);
        }
        flint_printf("\n");
 
        completed_steps++;
    }
 
finish:

    //clear all variables and arrays
    acb_clear(logk);
    acb_clear(coeff);
    acb_poly_clear(tpoly);
    acb_poly_clear(ypoly);
    acb_poly_clear(rpoly);
    acb_poly_clear(rdpoly);
    acb_poly_clear(Npoly);
 
    for (k = 1; k <= Nmax * D; k++)
    {
        acb_poly_clear(bt + k);
        acb_poly_clear(at + k);
    }
    flint_free(bt);
    flint_free(at);
 
    for (i = 0; i < divs; i++)
    {
        acb_poly_clear(moll + i + 1);
    }
    flint_free(d);
    flint_free(moll);
    flint_free(primes);
 
    return completed_steps;
}
 
//main initialisation - processing and validating input values 
int main(int argc, char *argv[])
{
    const char *t_str, *y_str;
    const char *Na_str, *Nb_str, *Nstep_str;
    const char *m_str, *n_str, *d_str;
    acb_t t, y;
    slong Na, Nb, Nstep;
    slong m, n, d, nderivs;
    slong total_steps, completed_steps;
    slong prec;
    int result = EXIT_SUCCESS;
    int usage_error = 0;
 
    acb_init(t);
    acb_init(y);
 
    if (argc != 9)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }
 
    t_str = argv[1];
    y_str = argv[2];
    Na_str = argv[3];
    Nb_str = argv[4];
    Nstep_str = argv[5];
    m_str = argv[6];
    n_str = argv[7];
    d_str = argv[8];
//    t_str = "0.4";
//    y_str = "0.2";
//    Na_str = "10";
//    Nb_str = "12";
//    Nstep_str = "1";
//    m_str = "3";
//    n_str = "0";
//    d_str = "30";
 
    Na = atol(Na_str);
    Nb = atol(Nb_str);
    Nstep = atol(Nstep_str);
    m = atol(m_str);
    nderivs = atol(n_str);
    d = atol(d_str);
 
    n = nderivs + 1;
 
    if (Na < 1 || Nb < 1 || Nstep < 1 || Na > Nb ||
        m < 0 || m > 20 || n < 0 || d < 1)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }
 
    prec = 16;
    total_steps = 1 + (Nb - Na) / Nstep;
    completed_steps = 0;
    while (completed_steps < total_steps)
    {
        acb_zero(t);
        acb_zero(y);
        arb_set_str(acb_realref(t), t_str, prec);
        arb_set_str(acb_realref(y), y_str, prec);
        completed_steps += run(
                t, y, Na + Nstep * completed_steps, Nb, Nstep, m, n, d, prec);
        prec *= 2;
    }
 
finish:
 
    if (usage_error && result == EXIT_FAILURE)
    {
        flint_printf("Usage:\n");
        flint_printf("%s t y Na Nb Nstep m n d\n\n", argv[0]);
        flint_printf(
    "This script computes a lower Lemma bound of abs(Aeff + Beff)/abs(Beff0)\n"
    "for N between 'Na' and 'Nb' with steps of size 'Nstep'.\n"
    "It uses an Euler mollification that includes m primes,\n"
    "so that for example m=0 corresponds to no mollification\n"
    "and m=3 corresponds to an Euler mollification that uses\n"
    "the first three primes {2, 3, 5}.\n"
    "Each row of output consists of N followed by the bound\n"
    "and n derivatives with respect to N, so that for example\n"
    "n=2 would print rows consisting of N, the bound,\n"
    "the first derivative of the bound with respect to N,\n"
    "and the second derivative of the bound with respect to N.\n"
    "All output has d significant decimal digits of accuracy.\n");
    }
 
    acb_clear(t);
    acb_clear(y);
 
    flint_cleanup();
    return result;
}
