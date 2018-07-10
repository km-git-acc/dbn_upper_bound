#include <stdio.h>
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

void acb_poly_onei(acb_poly_t res)
{
    acb_poly_one(res);
    acb_onei(acb_poly_get_coeff_ptr(res, 0));
}

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

//procedure to calculate the conjugate of an acb_poly function 
void
acb_poly_conj_series(acb_poly_t res,
        const acb_poly_t p, slong n, slong prec)
{
    acb_poly_t a, b;
    acb_t onei;
 
    acb_poly_init(a);
    acb_poly_init(b);
    acb_init(onei);
    acb_onei(onei);
 
    acb_poly_re(a, p);
    
    acb_poly_im(b, p);
    acb_poly_neg(b, b);
    acb_poly_scalar_mul(b, b, onei, prec);
    acb_poly_add_series(res, a, b, n, prec);

    acb_poly_clear(a);
    acb_poly_clear(b);
    acb_clear(onei);
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
abbeff_series(acb_t res, const arb_t rx, const arb_t ry, const arb_t rt, slong n, slong prec)
{
    slong N, sw;
    sw = 1;
    arb_t ar;
    acb_t one, pi, tmp, logk;
    acb_poly_t x, y, t;
    acb_poly_t s, sc, sprime, scprime;
    acb_poly_t alph1, alph2;
    acb_poly_t hs, hsc;
    acb_poly_t ta1, ta2;
    acb_poly_t A0, B0;
    acb_poly_t a, b, c, d, b1, b2;
    acb_poly_t summand;
    acb_poly_t Asum, Bsum;
    acb_poly_t A, B, C;
    acb_poly_t termC1, termC2;
    acb_poly_t alf, p, polyN, C0p, C0pc, U, Uc;
    acb_poly_t sig, T, Tprime;
    acb_poly_t onei, polyone, two, three, pip;

    acb_poly_init(a); 
    acb_poly_init(b); 
    acb_poly_init(c); 
    acb_poly_init(d); 
    acb_poly_init(x); 
    acb_poly_init(y); 
    acb_poly_init(t); 
    acb_poly_init(sig);
    acb_poly_init(T);
    acb_poly_init(Tprime);
    acb_poly_init(s);
    acb_poly_init(sc);
    acb_poly_init(sprime);
    acb_poly_init(scprime);
    acb_poly_init(alf); 
    acb_poly_init(p); 
    acb_poly_init(polyN); 
    acb_poly_init(C0p); 
    acb_poly_init(C0pc); 
    acb_poly_init(U);
    acb_poly_init(Uc);
    acb_poly_init(termC1);
    acb_poly_init(termC2);


    acb_poly_init(onei);
    acb_poly_onei(onei);

    acb_init(tmp);
    acb_init(one);
    acb_set_si(one, 1);
    acb_init(pi);
    acb_const_pi(pi, prec);

    acb_poly_init(polyone);
    acb_poly_set_si(polyone, 1);
    acb_poly_init(two);
    acb_poly_set_si(two, 2);
    acb_poly_init(three);
    acb_poly_set_si(three, 3);
    acb_poly_init(pip);
    acb_poly_set_acb(pip, pi);

    arb_init(ar);

    //convert arb-inputs to acb_poly variables
    acb_poly_one(x);
    arb_set(acb_realref(acb_poly_get_coeff_ptr(x, 0)), rx);
    acb_poly_one(y);
    arb_set(acb_realref(acb_poly_get_coeff_ptr(y, 0)), ry);
    acb_poly_one(t);
    arb_set(acb_realref(acb_poly_get_coeff_ptr(t, 0)), rt);

    //T
    acb_poly_scalar_mul_2exp_si(T, x, -1);

    //Tprime
    acb_poly_mullow(b, pip, t, n, prec);
    acb_poly_scalar_mul_2exp_si(b, b, -3);
    acb_poly_add_series(Tprime, b, T, n, prec);

    //N
    acb_poly_scalar_mul_2exp_si(a, pip, 1);
    acb_poly_div_series(a, Tprime, a, n, prec);
    acb_poly_sqrt_series(alf, a, n, prec);
    acb_poly_evaluate(tmp, alf, one, prec);
    acb_get_real(ar, tmp);
    arb_floor(ar, ar, prec);
    N = arf_get_d(arb_midref(ar), ARF_RND_DOWN);
    acb_poly_set_si(polyN, N);

    //s & sprime
    acb_poly_neg(a, y);
    acb_poly_add_si(a, a, 1, prec);
    acb_poly_scalar_mul_2exp_si(sig, a, -1);
    acb_poly_mullow(b, onei, T, n, prec);
    acb_poly_add_series(s, sig, b, n, prec);
    acb_poly_mullow(b, onei, Tprime, n, prec);
    acb_poly_add_series(sprime, sig, b, n, prec);

    acb_poly_add_si(sc, s, -1, prec);
    acb_poly_neg(sc, sc);

    acb_poly_add_si(scprime, sprime, -1, prec);
    acb_poly_neg(scprime, scprime);
 
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

    //calculate C-term if required
    acb_poly_init(C);
    acb_poly_zero(C);
    if (sw ==1)
    {
       //p
       acb_poly_sub_series(b, alf, polyN, n, prec);
       acb_poly_scalar_mul_2exp_si(b, b, 1);
       acb_poly_sub_series(p, polyone, b, n, prec);

       //U
       acb_poly_scalar_mul_2exp_si(a, pip, 1);
       acb_poly_div_series(a, Tprime, a, n, prec);
       acb_poly_log_series(a, a, n, prec);
       acb_poly_scalar_mul_2exp_si(b, Tprime, -1);
       acb_poly_mullow(a, a, b, n, prec);
       acb_poly_sub_series(a, a, b, n, prec);
       acb_poly_scalar_mul_2exp_si(b, pip, -3);
       acb_poly_sub_series(a, a, b, n, prec);
       acb_poly_mullow(a, a, onei, n, prec);
       acb_poly_neg(a, a);
       acb_poly_exp_series(U, a, n, prec);
       acb_poly_conj_series(Uc, U, n, prec);

       //C0p
       acb_poly_mullow(a, p, p, n, prec);
       acb_poly_scalar_mul_2exp_si(a, a, -1);
       acb_poly_scalar_mul_2exp_si(b, three, -3);
       acb_poly_add_series(a, a, b, n, prec);
       acb_poly_mullow(a, a, onei, n, prec);
       acb_poly_mullow(a, a, pip, n, prec);
       acb_poly_exp_series(a, a, n, prec);

       acb_poly_mullow(b, p, pip, n, prec);
       acb_poly_scalar_mul_2exp_si(b, b, -1);
       acb_poly_cos_series(b, b, n, prec);
       acb_poly_sqrt_series(c, two, n, prec);
       acb_poly_mullow(b, b, c, n, prec);
       acb_poly_mullow(b, b, onei, n, prec);
       acb_poly_sub_series(a, a, b, n, prec);

       acb_poly_mullow(b, p, pip, n, prec);
       acb_poly_cos_series(b, b, n, prec);
       acb_poly_scalar_mul_2exp_si(b, b, 1);
       acb_poly_div_series(C0p, a, b, n, prec);
       acb_poly_conj_series(C0pc, C0p, n, prec);


       //termC1
       acb_poly_scalar_mul_2exp_si(d, sprime, -1);
       acb_poly_neg(a, d);
       acb_poly_pow_series(a, pip, a, n, prec);

       acb_poly_gamma_series(b, d, n, prec);
       acb_poly_mullow(a, a, b, n, prec);

       acb_poly_neg(b, sig);
       acb_poly_pow_series(b, alf, b, n, prec);
       acb_poly_mullow(a, a, b, n, prec);

       acb_poly_mullow(b, C0p, U, n, prec);
       acb_poly_mullow(termC1, a, b, n, prec);

       //termC2
       acb_poly_scalar_mul_2exp_si(d, scprime, -1);
       acb_poly_neg(a, d);
       acb_poly_pow_series(a, pip, a, n, prec);

       acb_poly_gamma_series(b, d, n, prec);
       acb_poly_mullow(a, a, b, n, prec);

       acb_poly_add_si(b, sig, -1, prec);
       acb_poly_pow_series(b, alf, b, n, prec);
       acb_poly_mullow(a, a, b, n, prec);

       acb_poly_mullow(b, C0pc, Uc, n, prec);
       acb_poly_mullow(termC2, a, b, n, prec);

       //C
       acb_poly_mullow(a, pip, pip, n, prec);
       acb_poly_mullow(a, a, t, n, prec);
       acb_poly_scalar_mul_2exp_si(a, a, -6);
       acb_poly_exp_series(a, a, n, prec);

       acb_poly_add_si(b, sprime, -1, prec);
       acb_poly_mullow(b, b, sprime, n, prec);
       acb_poly_scalar_mul_2exp_si(b, b, -1);
       acb_poly_mullow(a, a, b, n, prec);

       acb_poly_add_series(b, termC1, termC2, n, prec);
       acb_poly_mullow(a, a, b, n, prec);

       acb_poly_set_si(b, -1);
       acb_poly_pow_series(b, b, polyN, n, prec);
       acb_poly_mullow(a, a, b, n, prec);
       acb_poly_scalar_mul_2exp_si(C, a, -3);
    }

    acb_poly_add_series(a, A, B, n, prec);
    acb_poly_sub_series(a, a, C, n, prec);

    acb_poly_div_series(a, a, B0, n, prec);

    acb_poly_evaluate(res, a, one, prec);
 
    arb_clear(ar);

    acb_clear(pi);
    acb_clear(one);
    acb_clear(tmp);
    acb_clear(logk);

    acb_poly_clear(x);
    acb_poly_clear(y);
    acb_poly_clear(t);
    acb_poly_clear(alph1);
    acb_poly_clear(alph2);
    acb_poly_clear(hs);
    acb_poly_clear(hsc);
    acb_poly_clear(ta1);
    acb_poly_clear(ta2);
    acb_poly_clear(A0);
    acb_poly_clear(B0);
    acb_poly_clear(a);
    acb_poly_clear(b);
    acb_poly_clear(c);
    acb_poly_clear(d);
    acb_poly_clear(b1);
    acb_poly_clear(b2);
    acb_poly_clear(summand);
    acb_poly_clear(Asum);
    acb_poly_clear(Bsum);
    acb_poly_clear(A);
    acb_poly_clear(B);
    acb_poly_clear(C);
    acb_poly_clear(onei);
    acb_poly_clear(two);
    acb_poly_clear(three);
    acb_poly_clear(T);
    acb_poly_clear(Tprime);
    acb_poly_clear(s);
    acb_poly_clear(sc);
    acb_poly_clear(sprime);
    acb_poly_clear(scprime);
	acb_poly_clear(U);
    acb_poly_clear(Uc);
    acb_poly_clear(C0p);
    acb_poly_clear(C0pc);
    acb_poly_clear(pip);
    acb_poly_clear(sig);
    acb_poly_clear(alf);
    acb_poly_clear(p);
    acb_poly_clear(termC1);
    acb_poly_clear(termC2);
    acb_poly_clear(polyN);
    acb_poly_clear(polyone);
}

//perform an incremental Muller step
void
acb_Muller_step(acb_t muller, acb_t zk, acb_t zk1, acb_t zk2, acb_t fk, acb_t fk1, acb_t fk2, arb_t accreq, slong prec)

{  
    acb_t a, b, c, d, q, A, B, C, diff1, diff2, diff3;
    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(d);
    acb_init(q);
    acb_init(A);
    acb_init(B);
    acb_init(C);
    acb_init(diff1);
    acb_init(diff2);
    acb_init(diff3);

    arb_t ar, br;
    arb_init(ar);
    arb_init(br);

    acb_zero(muller);

    //avoid drift 
    acb_get_real(ar, zk);
    acb_get_imag(br, zk);
    arb_get_mid_arb(ar, ar);
    arb_get_mid_arb(br, br);
    acb_set_arb_arb(zk, ar, br);

    acb_get_real(ar, zk1);
    acb_get_imag(br, zk1);
    arb_get_mid_arb(ar, ar);
    arb_get_mid_arb(br, br);
    acb_set_arb_arb(zk1, ar, br);

    acb_get_real(ar, zk2);
    acb_get_imag(br, zk2);
    arb_get_mid_arb(ar, ar);
    arb_get_mid_arb(br, br);
    acb_set_arb_arb(zk2, ar, br);

    acb_get_real(ar, fk);
    acb_get_imag(br, fk);
    arb_get_mid_arb(ar, ar);
    arb_get_mid_arb(br, br);
    acb_set_arb_arb(fk, ar, br);

    acb_get_real(ar, fk1);
    acb_get_imag(br, fk1);
    arb_get_mid_arb(ar, ar);
    arb_get_mid_arb(br, br);
    acb_set_arb_arb(fk1, ar, br);

    acb_get_real(ar, fk2);
    acb_get_imag(br, fk2);
    arb_get_mid_arb(ar, ar);
    arb_get_mid_arb(br, br);
    acb_set_arb_arb(fk2, ar, br);

    //if the three function values become equal, break of the process (Muller method can't handle this)
    acb_sub(diff1, fk, fk1, prec);
    acb_sub(diff2, fk1, fk2, prec);
    acb_add(diff3, diff1, diff2, prec);
    acb_abs(ar, diff3, prec);
    if (arb_lt(ar, accreq))
    {
       acb_set_si(muller, -1);
       goto endmulstep;
    }

    acb_sub(a, zk, zk1, prec);
    acb_sub(b, zk1, zk2, prec);
    acb_div(q, a, b, prec);

    acb_mul(A, q, fk, prec);
    acb_add_si(a, q, 1, prec);
    acb_mul(a, a, fk1, prec);
    acb_mul(a, a, q, prec);
    acb_sub(A, A, a, prec);
    acb_mul(a, q, q, prec);
    acb_mul(a, a, fk2, prec);
    acb_add(A, A, a, prec); 

    acb_mul_2exp_si(a, q, 1);
    acb_add_si(a, a, 1, prec);
    acb_mul(B, a, fk, prec);
    acb_add_si(a, q, 1, prec);
    acb_mul(b, a, a, prec);
    acb_mul(b, b, fk1, prec);
    acb_sub(B, B, b, prec);
    acb_mul(a, q, q, prec);
    acb_mul(a, a, fk2, prec);
    acb_add(B, B, a, prec);

    acb_add_si(C, q, 1, prec);
    acb_mul(C, C, fk, prec);

    acb_mul(a, B, B, prec);
    acb_mul_2exp_si(b, A, 2);
    acb_mul(b, b, C, prec);
    acb_sub(a, a, b, prec);

    //avoid drift to allow sqrt to work for negative numbers
    acb_get_real(ar, a);
    acb_get_imag(br, a);
    arb_get_mid_arb(ar, ar);
    arb_get_mid_arb(br, br);
    acb_set_arb_arb(a, ar, br);

    acb_sqrt(d, a, prec);

    //avoid drift
    acb_add(a, B, d, prec);
    acb_abs(ar, a, prec);
    arb_mul(ar, ar, ar, prec);
    acb_sub(b, B, d, prec);
    acb_abs(br, b, prec);
    arb_mul(br, br, br, prec);

    if (arb_gt(ar, br))
    {
        acb_sub(c, zk, zk1, prec);
        acb_mul_2exp_si(c, c, 1);
        acb_mul(c, c, C, prec);
        acb_div(c, c, a, prec);
        acb_sub(muller, zk, c, prec);
    }
    else
    {
        acb_sub(c, zk, zk1, prec);
        acb_mul_2exp_si(c, c, 1);
        acb_mul(c, c, C, prec);
        acb_div(c, c, b, prec);
        acb_sub(muller, zk, c, prec);       
    }

endmulstep:


    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(d);
    acb_clear(q);
    acb_clear(A);
    acb_clear(B);
    acb_clear(C);
    acb_clear(diff1);
    acb_clear(diff2);
    acb_clear(diff3);

    arb_clear(ar);
    arb_clear(br);
}

//execute Muller rootfinding algorithm
void
acb_Muller(acb_t root, acb_t z1, acb_t z2, acb_t z3, arb_t t, arb_t xtestplus10, arb_t accreq, slong nmax, slong prec)

{  
    acb_mat_t xx, f;
    acb_mat_init(xx, nmax+4, 1);
    acb_mat_init(f, nmax+4, 1);

    acb_t a, b, Ht;
    acb_init(a);
    acb_init(b);
    acb_init(Ht);

    arb_t one, abs, realz, imagz, diff;
    arb_init(one);
    arb_init(abs);
    arb_init(realz);
    arb_init(imagz);
    arb_init(diff);

    slong k, rootnotfound;

    //initialise values
    k=3;

    acb_set(acb_mat_entry(xx, 1, 0), z1);
    acb_set(acb_mat_entry(xx, 2, 0), z2);
    acb_set(acb_mat_entry(xx, 3, 0), z3);

    //fill the 3 function values but ensure x and y stay positive
    acb_get_real(realz, z1);
    arb_abs(realz, realz);
    acb_get_imag(imagz, z1);
    arb_abs(imagz, imagz);
    abbeff_series(Ht, realz, imagz, t, 1, prec);
    acb_set(acb_mat_entry(f, 1, 0), Ht);

    acb_get_real(realz, z2);
    arb_abs(realz, realz);
    acb_get_imag(imagz, z2);
    arb_abs(imagz, imagz);
    abbeff_series(Ht, realz, imagz, t, 1, prec);
    acb_set(acb_mat_entry(f, 2, 0), Ht);

    acb_get_real(realz, z3);
    arb_abs(realz, realz);
    acb_get_imag(imagz, z3);
    arb_abs(imagz, imagz);
    abbeff_series(Ht, realz, imagz, t, 1, prec);
    acb_set(acb_mat_entry(f, 3, 0), Ht);


    acb_sub(a, acb_mat_entry(xx, k, 0), acb_mat_entry(xx, k-1, 0), prec);
    acb_abs(diff, a, prec);

    rootnotfound = 0;

    while ((k < nmax+3) && arb_gt(diff, accreq))
    {       
        acb_Muller_step(a, acb_mat_entry(xx, k, 0), acb_mat_entry(xx, k-1, 0), acb_mat_entry(xx, k-2, 0),
                       acb_mat_entry(f, k, 0), acb_mat_entry(f, k-1, 0), acb_mat_entry(f, k-2, 0), accreq, prec);

        acb_abs(abs, a, prec);
        arb_one(one);

        //if the process 'explodes' > test + 10 then ignore this zero and label it as 'notfound' (root=0)
        if (arb_gt(abs, xtestplus10) || arb_eq(abs, one))
        {
            rootnotfound = 1;
            goto exitwhile;
        }

        acb_get_real(realz, a);
        arb_abs(realz, realz);
        acb_get_imag(imagz, a);
        arb_abs(imagz, imagz);
        acb_set_arb_arb(acb_mat_entry(xx, k+1, 0), realz, imagz);
        abbeff_series(Ht, realz, imagz, t, 1, prec);
        acb_set(acb_mat_entry(f, k+1, 0), Ht);

        acb_sub(a, acb_mat_entry(xx, k, 0), acb_mat_entry(xx, k-1, 0), prec);
        acb_abs(diff, a, prec);

        arb_get_mid_arb(diff,diff);

        k = k + 1;
    }

exitwhile:

    if (rootnotfound == 1 || k == nmax + 3)
       acb_zero(root);
    else
       acb_set(root, acb_mat_entry(xx, k - 1, 0));

    acb_mat_clear(xx);
    acb_mat_clear(f);

    acb_clear(a);
    acb_clear(b);
    acb_clear(Ht);

    arb_clear(abs);
    arb_clear(one);
    arb_clear(realz);
    arb_clear(imagz);
    arb_clear(diff);
}

//loop through all t and x values
void
acb_Muller_loop(slong res, arb_t xa, arb_t xb, arb_t ya, arb_t yb, arb_t ta, arb_t tb,
           arb_t tstep, slong digits, slong prec)

{  
    slong nmax;
    nmax = 200;

    arb_t a, b, t, xr, yr, xrprev, yrprev, xdiff, ydiff, accreq, xtest, xtest4, xinc;
    arb_init(a);
    arb_init(b);
    arb_init(t);
    arb_init(xr);
    arb_init(yr);
    arb_init(xrprev);
    arb_init(yrprev);
    arb_init(xdiff);
    arb_init(ydiff);
    arb_init(accreq);
    arb_init(xtest);
    arb_init(xtest4);
    arb_init(xinc);

    acb_t z1, z2, z3, root;
    acb_init(z1);
    acb_init(z2);
    acb_init(z3);
    acb_init(root);
	
    //establish required accuracy
    arb_set_si(a, digits);
    arb_neg(a, a);
    arb_set_si(b, 10);
    arb_pow(accreq, b, a, prec);

    //max difference between two subsequent zeros
    arb_set_si(a, 1);
    arb_set_si(b, 10);
    arb_div(xinc, a, b, prec);

    arb_get_mid_arb(ta, ta);
    arb_set(t, ta);    
    arb_get_mid_arb(tb, tb);
    while(arb_le(t, tb))
    {
        arb_zero(xrprev);
        arb_zero(yrprev);
        arb_add_si(xtest, xa, -1, prec);
        while(arb_le(xtest, xb))
        {
            acb_set_arb(z1, xtest);
            arb_set_si(a, 5);
            arb_set_si(b, 10);
            arb_div(a, a, b, prec);
            acb_add_arb(z2, z1, a, prec);
            acb_add_si(z3, z1, 2, prec);

            //set a limit value (xb+4) above which the Muller process is considered to be derailing.
            arb_add_si(a, xtest, 10, prec);

            acb_Muller(root, z1, z2, z3, t, a, accreq, nmax, prec);
            acb_get_real(xr, root);
            acb_get_imag(yr, root);

            arb_sub(xdiff, xr, xrprev, prec);
            arb_abs(xdiff, xdiff);
            arb_sub(ydiff, yr, yrprev, prec);
            arb_abs(ydiff, ydiff);

            arb_add_si(xtest4, xtest, 4, prec);

            //if the root has not been found or is the same as the previous or is outside the Muller "window" then search+1
            if (arb_is_zero(xr) || (arb_lt(xdiff, accreq) && arb_lt(ydiff, accreq)) || arb_lt(xr, xtest) || arb_gt(xr, xtest4))
            {
                arb_add(xtest, xtest, xinc, prec);
            }
            //else a new root has been found, so set the xtest a tiny bit (xinc) forward (e.g to find a Lehmer pair)
            else
            {
                
                arb_add(xtest, xtest, xinc, prec); 
                //only print the root if doesn't exceed xb
                if (arb_le(xr, xb))
                {
                    //make y properly 0 if it is smaller than the accreq
                    if (arb_lt(yr, accreq))
                    arb_zero(yr);

                    arb_printn(xr, digits, ARB_STR_NO_RADIUS);
                    printf(", ");
                    arb_printn(yr, digits, ARB_STR_NO_RADIUS);
                    printf(", ");
                    double pt = arf_get_d(arb_midref(t), ARF_RND_NEAR);
                    printf("%1.5f", pt);
                    printf("\n");

                    //optional code to verify that indeed a coorect zero has been found
                    //abbeff_series(root, xr, yr, t, 1, prec);
                    //acb_printd(root, 20);
                    //printf("\n");

                    //print the conjugate if the root is complex
                    arb_zero(a);
                    if (arb_gt(yr, a))
                    {
                        arb_printn(xr, digits, ARB_STR_NO_RADIUS);
                        printf(", ");
                        arb_neg(b, yr);
                        arb_printn(b, digits, ARB_STR_NO_RADIUS);
                        printf(", ");
                        printf("%1.5f", pt);
                        printf("\n");
                    }

                    //save the root values to test for double finds
                    arb_set(xrprev, xr);
                    arb_set(yrprev, yr);
                }
           
            }
         //avoid drift in error term
         arb_get_mid_arb(xtest, xtest);     
        }

      arb_add(t, t, tstep, prec);

      //make t an exact decimal for correct comparison with tb
      arb_abs(a, t);
      if (arb_lt(a, accreq))
         arb_zero(t); 
      arb_get_mid_arb(t, t);
    }

    acb_clear(root);
    acb_clear(z1);
    acb_clear(z2);
    acb_clear(z3);

    arb_clear(a);
    arb_clear(b);
    arb_clear(t);
    arb_clear(xr);
    arb_clear(yr);
    arb_clear(xrprev);
    arb_clear(yrprev);
    arb_clear(xdiff);
    arb_clear(ydiff);
    arb_clear(accreq);
    arb_clear(xtest);
    arb_clear(xtest4);
    arb_clear(xinc);
}
 
int main(int argc, char *argv[])
{
    arb_t xa, xb, ya, yb, ta, tb, tstep;
    arb_init(xa);
    arb_init(xb);
    arb_init(ya);
    arb_init(yb);
    arb_init(ta);
    arb_init(tb);
    arb_init(tstep);

    acb_t outcome;
    acb_init(outcome);

    const char *xa_str, *xb_str, *ya_str, *yb_str, *ta_str, *tb_str, *tstep_str, *digits_str;
    slong prec, digits, res;
    int result = EXIT_SUCCESS;
    res=0;

    if (argc != 7)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    xa_str = argv[1];
    xb_str = argv[2];
    ya_str = "0";
    yb_str = "0";
    ta_str = argv[3];
    tb_str = argv[4];
    tstep_str = argv[5];
    digits_str = argv[6];


    //precision 
    digits = atol(digits_str);
    prec = digits * 3.32192809488736 + 30;

    arb_set_str(xa, xa_str, prec);
    arb_set_str(xb, xb_str, prec);
    arb_set_str(ya, ya_str, prec);
    arb_set_str(yb, yb_str, prec);
    arb_set_str(ta, ta_str, prec);
    arb_set_str(tb, tb_str, prec);
    arb_set_str(tstep, tstep_str, prec);

//TIMEIT_ONCE_START
    acb_Muller_loop(res, xa, xb, ya, yb, ta, tb, tstep, digits, prec);
//TIMEIT_ONCE_STOP

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Required inputs:\n");
        flint_printf("%s xa, xb, ta, tb, tstep, digits \n\n", argv[0]);
        flint_printf(
    "This script attempts to find all real and complex roots \n"
    "of Ht_eff(x+iy)/B0 in the domain xa..xb, ta...tb,\n" 
    "Accuracy can be set by digits. \n");
    }

    arb_clear(xa);
    arb_clear(ya);
    arb_clear(ta);
    arb_clear(xb);
    arb_clear(yb);
    arb_clear(tb);
    arb_clear(tstep);

    acb_clear(outcome);
 
    flint_cleanup();

    return result;
}
