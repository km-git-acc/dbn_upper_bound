/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
#include "acb_poly.h"

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

//procedure to calculate the kappa-value
void
_acb_poly_kappa_series(acb_poly_t res, const acb_poly_t x,
        const acb_poly_t y, const acb_poly_t t, slong n, slong prec)
{
    acb_poly_t a, b;
    acb_poly_init(a);
    acb_poly_init(b);

    acb_poly_mullow(a, t, y, n, prec);
    acb_poly_add_si(b, x, -6, prec);
    acb_poly_scalar_mul_2exp_si(b, b, 1);
    acb_poly_div_series(res, a, b, n, prec);
 
    acb_poly_clear(a);
    acb_poly_clear(b);
}

//procedure to calculate the Ress-value
void
_acb_poly_Ress_series(acb_poly_t res, const acb_poly_t x,
        const acb_poly_t y, const acb_poly_t t, slong n, slong prec)
{
    acb_poly_t a, b, c, d, e, f, g, z;
    acb_poly_init(a);
    acb_poly_init(b);
    acb_poly_init(c);
    acb_poly_init(d);
    acb_poly_init(e);
    acb_poly_init(f);
    acb_poly_init(g);
    acb_poly_init(z);
    acb_poly_zero(z);

    acb_t pi, three;
    acb_init(three);
    acb_set_ui(three, 3);
    acb_init(pi);
    acb_const_pi(pi, prec);

    acb_poly_add_si(a, y, 1, prec);
    acb_poly_scalar_mul_2exp_si(a, a, -1);

    acb_poly_scalar_div(b, x, pi, prec);
    acb_poly_scalar_mul_2exp_si(b, b, -2);
    acb_poly_log_series(b, b, n, prec);
    acb_poly_scalar_mul_2exp_si(c, t, -2);
    acb_poly_mullow(c, c, b, n, prec);

    acb_poly_add_series(c, a, c, n, prec);

    acb_poly_scalar_mul_2exp_si(d, y, 2);
    acb_poly_add_si(e, y, 1, prec);
    acb_poly_mullow(e, e, d, n, prec);
    acb_poly_mullow(f, x, x, n, prec);
    acb_poly_div_series(f, e, f, n, prec);

    acb_poly_neg(g, y);
    acb_poly_scalar_mul(g, g, three, prec);
    acb_poly_add_si(g, g, 1, prec);
    acb_poly_add_series(g, g, f, n, prec);

    acb_poly_max_series(g, g, z, n, prec);
    acb_poly_mullow(g, g, t, n, prec);

    acb_poly_mullow(a, x, x, n, prec);
    acb_poly_scalar_mul_2exp_si(a, a, 1);
    acb_poly_div_series(g, g, a, n, prec);

    acb_poly_sub_series(res, c, g, n, prec);
 
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

//procedure to print values of abbeff for the mesh 
static void
print_mesh_series(arb_t res, const arb_t x, const arb_t y, const arb_t t,const arb_t absabb, slong digits)
{

        arf_printd(arb_midref(x), digits);
        printf(", ");
        arf_printd(arb_midref(y), digits);
        printf(", ");
        arf_printd(arb_midref(t), digits);
        printf(", ");
        arf_printd(arb_midref(absabb), digits);
        printf("\n");
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
 
void
ddxbound_series(arb_t res, const arb_t rx, const arb_t ry, const arb_t rt, const slong N, slong n, slong prec)
{
    acb_t one, acbres, logk, pi, onei;
    acb_poly_t a, b, c, d, x ,y, t, kterm, kterm1, Nterm, commsum, pipoly, polyN;
    acb_poly_t bt, ress, modgam, kappa;
    acb_poly_t Asum, Bsum, summandA, summandB, summandAfix;

    acb_init(one);
    acb_init(acbres);
    acb_init(logk);
    acb_init(onei);
    acb_onei(onei);
    acb_init(pi);
    acb_const_pi(pi, prec);

    acb_poly_init(pipoly);
    acb_poly_set_acb(pipoly, pi);
    acb_poly_init(summandA);
    acb_poly_init(summandB);
    acb_poly_init(summandAfix);
    acb_poly_init(Asum);
    acb_poly_init(Bsum);
    acb_poly_init(a);
    acb_poly_init(b);
    acb_poly_init(c);
    acb_poly_init(d);
    acb_poly_init(bt);
    acb_poly_init(ress);
    acb_poly_init(modgam);
    acb_poly_init(kappa);
    acb_poly_init(kterm);
    acb_poly_init(kterm1);
    acb_poly_init(Nterm);
    acb_poly_init(commsum);
    acb_poly_init(polyN);
    acb_poly_zero(polyN);
    acb_poly_add_si(polyN, polyN, N, prec);

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

    _acb_poly_Ress_series(ress, x, y, t, n, prec);
    _acb_poly_modgamma_series(modgam, x, y, n, prec);
    acb_poly_abs_series(modgam, modgam, n, prec);
    _acb_poly_kappa_series(kappa, x, y, t, n, prec);
    acb_poly_abs_series(kappa, kappa, n, prec);

    // common factor t/(4*(x-6))
    acb_poly_add_si(commsum, x, -6, prec);
    acb_poly_scalar_mul_2exp_si(commsum, commsum, 2);
    acb_poly_div_series(commsum, t, commsum, n, prec);

    // fill all k-independent terms for summandA outside loop
    acb_poly_add_si(b, y, 1, prec);
    acb_poly_scalar_mul(c, x, onei, prec);
    acb_poly_add_series(c, c, b, n, prec);
    acb_poly_abs_series(c, c, n, prec);
    acb_poly_div_series(c, c, pipoly, n, prec);
    acb_poly_scalar_mul_2exp_si(c, c, -2);
    acb_poly_log_series(c, c, n, prec);
    acb_poly_add_series(c, c, pipoly, n, prec);
    acb_poly_zero(d);
    acb_poly_add_si(d, d, 3, prec);
    acb_poly_div_series(d, d, x, n, prec);
    acb_poly_add_series(d, d, c, n, prec);
    acb_poly_one(c);
    acb_poly_scalar_mul_2exp_si(c, c, -1);
    acb_poly_add_series(c, c, commsum, n, prec);
    acb_poly_mullow(summandAfix, d, c, n, prec);

    {
        acb_poly_zero(Asum);
        acb_poly_zero(Bsum);
        slong k;
        for (k = 1; k <= N; k++)
        {
            acb_set_si(logk, k);
            acb_log(logk, logk, prec);
            _acb_poly_bt_series(bt, logk, t, n, prec);

            // bt/k^Ress            
            acb_poly_scalar_mul(kterm, ress, logk, prec);
            acb_poly_neg(kterm, kterm);
            acb_poly_exp_series(kterm, kterm, n, prec);
            acb_poly_mullow(kterm, kterm, bt, n, prec);

            //summand B
            acb_poly_one(a);
            acb_poly_scalar_mul(a, a, logk, prec);
            acb_poly_scalar_mul_2exp_si(a, a, -1);

            acb_poly_scalar_mul(b, commsum, logk, prec);
            acb_poly_add_series(b, b, a, n, prec);
            acb_poly_mullow(summandB, kterm, b, n, prec);
            acb_poly_add_series(Bsum, Bsum, summandB, n, prec);

            //summand A
            acb_poly_scalar_mul(kterm1, y, logk, prec);
            acb_poly_exp_series(kterm1, kterm1, n, prec);
            acb_poly_mullow(kterm1, kterm1, kterm, n, prec);

            acb_poly_scalar_mul(a, commsum, logk, prec);
            acb_poly_add_series(summandA, summandAfix, a, n, prec);
            acb_poly_mullow(summandA, kterm1, summandA, n, prec);
            acb_poly_add_series(Asum, Asum, summandA, n, prec);
        }
    }
	
    // Bsum - Asum*modgam*N^kappa 
    acb_poly_pow_series(Nterm, polyN, kappa, n, prec);
    acb_poly_mullow(Nterm, Nterm, Asum, n, prec);
    acb_poly_mullow(Nterm, Nterm, modgam, n, prec);
    acb_poly_add_series(a, Bsum, Nterm, n, prec);

    acb_set_si(one, 1);
    acb_poly_evaluate(acbres, a, one, prec);
    acb_abs(res, acbres, prec);

    acb_clear(one);
    acb_clear(acbres);   
    acb_clear(logk);
    acb_clear(pi);
    acb_clear(onei);
    acb_poly_clear(a);
    acb_poly_clear(b);
    acb_poly_clear(c);
    acb_poly_clear(d);
    acb_poly_clear(x);
    acb_poly_clear(y);
    acb_poly_clear(t);
    acb_poly_clear(bt);
    acb_poly_clear(ress);  
    acb_poly_clear(modgam);
    acb_poly_clear(kappa); 
    acb_poly_clear(summandA);
    acb_poly_clear(summandB);
    acb_poly_clear(summandAfix);
    acb_poly_clear(Asum);
    acb_poly_clear(Bsum);
    acb_poly_clear(kterm);
    acb_poly_clear(kterm1);
    acb_poly_clear(Nterm);
    acb_poly_clear(commsum);
    acb_poly_clear(polyN);
    acb_poly_clear(pipoly);

} 

void
ddtbound_series(arb_t res, const arb_t rx, const arb_t ry, const arb_t rt, const slong N, slong n, slong prec)
{
    acb_t  one, acbres, logk, pi, onei;
    acb_poly_t a, b, c, d, x, y, t, kterm, kterm1, Nterm, commsum, pipoly, kpoly, polyN;
    acb_poly_t bt, ress, modgam, kappa;
    acb_poly_t Asum, Bsum, summandA, summandB, summandAfix;

    acb_init(one);
    acb_init(acbres);
    acb_init(logk);
    acb_init(onei);
    acb_onei(onei);
    acb_init(pi);
    acb_const_pi(pi, prec);

    acb_poly_init(pipoly);
    acb_poly_set_acb(pipoly, pi);
    acb_poly_init(summandA);
    acb_poly_init(summandB);
    acb_poly_init(summandAfix);
    acb_poly_init(Asum);
    acb_poly_init(Bsum);
    acb_poly_init(a);
    acb_poly_init(b);
    acb_poly_init(c);
    acb_poly_init(d);
    acb_poly_init(bt);
    acb_poly_init(ress);
    acb_poly_init(modgam);
    acb_poly_init(kappa);
    acb_poly_init(kterm);
    acb_poly_init(kterm1);
    acb_poly_init(Nterm);
    acb_poly_init(commsum);
    acb_poly_init(kpoly);
    acb_poly_init(polyN);
    acb_poly_zero(polyN);
    acb_poly_add_si(polyN, polyN, N, prec);

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

    _acb_poly_Ress_series(ress, x, y, t, n, prec);
    _acb_poly_modgamma_series(modgam, x, y, n, prec);
    acb_poly_abs_series(modgam, modgam, n, prec);
    _acb_poly_kappa_series(kappa, x, y, t, n, prec);
    acb_poly_abs_series(kappa, kappa, n, prec);

    // common factor 2/(x-6)
	acb_poly_one(a);
    acb_poly_scalar_mul_2exp_si(a, a, 1);
    acb_poly_add_si(commsum, x, -6, prec);
    acb_poly_div_series(commsum, a, commsum, n, prec);

    // fill all k-independent terms for summandA outside loop
    acb_poly_scalar_mul_2exp_si(a, commsum, 2);
    acb_poly_scalar_mul_2exp_si(b, pipoly, -1);
    acb_poly_add_series(b, a, b, n, prec);
    acb_poly_scalar_mul_2exp_si(b, b, -2);

    acb_poly_scalar_mul_2exp_si(c, pipoly, 2);
    acb_poly_div_series(c, x, c, n, prec);
    acb_poly_log_series(c, c, n, prec);
    acb_poly_scalar_mul_2exp_si(d, commsum, 2);
    acb_poly_add_series(d, d, c, n, prec);

    acb_poly_mullow(summandAfix, d, b, n, prec);

    {
        acb_poly_zero(Asum);
        acb_poly_zero(Bsum);
        slong k;
        for (k = 1; k <= N; k++)
        {
            acb_poly_set_si(kpoly, k);
            acb_set_si(logk, k);
            acb_log(logk, logk, prec);
            _acb_poly_bt_series(bt, logk, t, n, prec);

            // bt/k^Ress            
            acb_poly_scalar_mul(kterm, ress, logk, prec);
            acb_poly_neg(kterm, kterm);
            acb_poly_exp_series(kterm, kterm, n, prec);
            acb_poly_mullow(kterm, kterm, bt, n, prec);

            //summand B
            acb_poly_one(a);
            acb_poly_scalar_mul(a, a, logk, prec);
            acb_poly_scalar_mul_2exp_si(a, a, -2);

            acb_poly_mullow(b, pipoly, kpoly, n, prec);
            acb_poly_div_series(c, x, b, n, prec);
            acb_poly_scalar_mul_2exp_si(c, c, -2);
            acb_poly_log_series(c, c, n, prec);
            acb_poly_mullow(c, c, a, n, prec);

            acb_poly_scalar_mul_2exp_si(b, pipoly, -3);
            acb_poly_scalar_mul(b, b, logk, prec);

            acb_poly_scalar_mul(d, commsum, logk, prec);

            acb_poly_add_series(d, d, b, n, prec);
            acb_poly_add_series(d, d, c, n, prec);
            acb_poly_mullow(summandB, kterm, d, n, prec);
            acb_poly_add_series(Bsum, Bsum, summandB, n, prec);

            //summand A
            acb_poly_scalar_mul(kterm1, y, logk, prec);
            acb_poly_exp_series(kterm1, kterm1, n, prec);
            acb_poly_mullow(kterm1, kterm1, kterm, n, prec);

            acb_poly_one(a);
            acb_poly_scalar_mul(a, a, logk, prec);
            acb_poly_scalar_mul_2exp_si(a, a, -2);

            acb_poly_mullow(b, pipoly, kpoly, n, prec);
            acb_poly_div_series(c, x, b, n, prec);
            acb_poly_scalar_mul_2exp_si(c, c, -2);
            acb_poly_log_series(c, c, n, prec);
            acb_poly_mullow(c, c, a, n, prec);

            acb_poly_scalar_mul_2exp_si(b, pipoly, -3);
            acb_poly_scalar_mul(b, b, logk, prec);

            acb_poly_scalar_mul(d, commsum, logk, prec);

            acb_poly_add_series(d, d, b, n, prec);
            acb_poly_add_series(d, d, c, n, prec);
            acb_poly_add_series(summandA, summandAfix, d, n, prec);
            acb_poly_mullow(summandA, kterm1, summandA, n, prec);
            acb_poly_add_series(Asum, Asum, summandA, n, prec);
        }
    }
	
    // Bsum - Asum*modgam*N^kappa
    acb_poly_pow_series(Nterm, polyN, kappa, n, prec);
    acb_poly_mullow(Nterm, Nterm, Asum, n, prec);
    acb_poly_mullow(Nterm, Nterm, modgam, n, prec);
    acb_poly_add_series(a, Bsum, Nterm, n, prec);

    acb_set_si(one, 1);
    acb_poly_evaluate(acbres, a, one, prec);
    acb_abs(res, acbres, prec);

    acb_clear(one);
    acb_clear(acbres);
    acb_clear(logk);
    acb_clear(pi);
    acb_clear(onei);
    acb_poly_clear(a);
    acb_poly_clear(b);
    acb_poly_clear(c);
    acb_poly_clear(d);
    acb_poly_clear(x);
    acb_poly_clear(y);
    acb_poly_clear(t);
    acb_poly_clear(bt);
    acb_poly_clear(ress);  
    acb_poly_clear(modgam);
    acb_poly_clear(kappa); 
    acb_poly_clear(summandA);
    acb_poly_clear(summandB);
    acb_poly_clear(summandAfix);
    acb_poly_clear(Asum);
    acb_poly_clear(Bsum);
    acb_poly_clear(kterm);
    acb_poly_clear(kterm1);
    acb_poly_clear(Nterm);
    acb_poly_clear(commsum);
    acb_poly_clear(kpoly);
    acb_poly_clear(polyN);
    acb_poly_clear(pipoly);

} 

void
mesh_series(arb_t res, const arb_t xa, arb_t xb, const arb_t ya, const arb_t yb, const arb_t ta,
            const arb_t tb, const arb_t c, const slong N, slong n, slong prec, slong digits)
{
    arb_t a, x, y, t, pi, arg, absabb, dtabb, dxabb, windnum; 
    acb_t abb, abbstart, abbprev, argdiv; 

    arb_init(a);
    arb_init(x);
    arb_init(y);
    arb_init(t);
    arb_init(dtabb);
    arb_init(dxabb);
    arb_init(absabb);
    arb_init(arg);
    arb_init(windnum);

    arb_init(pi);
    arb_const_pi(pi, prec);

    acb_init(abb);
    acb_init(abbstart);
    acb_init(abbprev);
    acb_init(argdiv);

    arb_set(t,ta);

    //init the 'walk'
    while (arb_le(t,tb))
    {
        ddxbound_series(dxabb, xa, ya, t, N, n, prec);
        ddtbound_series(dtabb, xa, ya, t, N, n, prec);
        abbeff_series(abbstart, xa, ya, t, N, n, prec);
        acb_set(abbprev, abbstart);
        arb_zero(windnum);

             //walk left of the rectangle
             arb_set(x,xa);
             arb_set(y,ya);
             while (arb_le(y,yb))
             {
                  abbeff_series(abb, x, y, t, N, n, prec);
                  acb_div(argdiv,abb,abbprev, prec);
                  acb_arg(arg, argdiv, prec);
                  arb_add(windnum,windnum,arg,prec);
                  acb_abs(absabb,abb, prec); 
                  arb_sub(a,absabb,c, prec);
                  arb_div(a,a,dxabb, prec);
                  arb_get_mid_arb(a,a);
                  arb_add(y, y, a, prec);
                  acb_set(abbprev, abb);
                  print_mesh_series(a, x, y, t, absabb, digits);
             }

             //walk top of the rectangle
             arb_set(x,xa);
             arb_set(y,yb);
             while (arb_le(x,xb))
             {
                  abbeff_series(abb, x, y, t, N, n, prec);
                  acb_div(argdiv,abb,abbprev, prec);
                  acb_arg(arg, argdiv, prec);
                  arb_add(windnum,windnum,arg,prec);
                  acb_abs(absabb,abb, prec); 
                  arb_sub(a,absabb,c, prec);
                  arb_div(a,a,dxabb, prec);
                  arb_get_mid_arb(a,a);
                  arb_add(x,x,a,prec);
                  acb_set(abbprev, abb);
                  print_mesh_series(a, x, y, t, absabb, digits);
             }

             //walk right of the rectangle
             arb_set(x,xb);
             arb_set(y,yb);
             while (arb_ge(y,ya))
             {
                  abbeff_series(abb, x, y, t, N, n, prec);
                  acb_div(argdiv,abb,abbprev, prec);
                  acb_arg(arg, argdiv, prec);
                  arb_add(windnum,windnum,arg,prec);
                  acb_abs(absabb,abb, prec); 
                  arb_sub(a,absabb,c, prec);
                  arb_div(a,a,dxabb, prec);
                  arb_get_mid_arb(a,a);
                  arb_sub(y,y,a,prec);
                  acb_set(abbprev, abb);
                  print_mesh_series(a, x, y, t, absabb, digits);
             }

             //walk bottom of the rectangle
             arb_set(x,xb);
             arb_set(y,ya);
             while (arb_ge(x,xa))
             {
                  abbeff_series(abb, x, y, t, N, n, prec);
                  acb_div(argdiv,abb,abbprev, prec);
                  acb_arg(arg, argdiv, prec);
                  arb_add(windnum,windnum,arg,prec);
                  acb_abs(absabb,abb, prec); 
                  arb_sub(a,absabb,c, prec);
                  arb_div(a,a,dxabb, prec);
                  arb_get_mid_arb(a,a);
                  arb_sub(x,x,a,prec);
                  acb_set(abbprev, abb);
                  print_mesh_series(a, x, y, t, absabb, digits);
             }

        acb_div(argdiv,abbstart,abbprev, prec);
        acb_arg(arg, argdiv, prec);
        arb_add(windnum,windnum,arg,prec);
        arb_div(windnum, windnum, pi, prec);
        arb_mul_2exp_si(windnum, windnum, -1);

        arb_neg(windnum, windnum);
        arb_set_round(windnum, windnum, prec);

        arb_div(a, c, dtabb, prec);
        arb_add(t, t, a, prec);
    }

    arb_set(res,windnum);

    arb_clear(pi);
    arb_clear(a);
    arb_clear(x);
    arb_clear(y);
    arb_clear(t);
    arb_clear(arg);
    arb_clear(dtabb);
    arb_clear(dxabb);
    arb_clear(absabb);
    arb_clear(windnum);
    acb_clear(abb);
    acb_clear(abbstart);
    acb_clear(abbprev);
    acb_clear(argdiv);

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
    arb_t X, xa, xb, ya, yb, ta, tb, a, c, d, re, im, out, windtot, z;
    acb_t out1;
    const char *X_str, *ya_str, *yb_str, *ta_str, *tb_str, *c_str;
    slong N, prec, digits;
    int result = EXIT_SUCCESS;
    int usage_error = 0;
 
    arb_init(X);
    arb_init(xa);
    arb_init(xb);
    arb_init(ya);
    arb_init(yb);
    arb_init(ta);
    arb_init(tb);
    arb_init(a);
    arb_init(c);
    arb_init(d);
    arb_init(z);
    arb_init(re);
    arb_init(im);
    arb_init(windtot);

    arb_init(out);
    acb_init(out1);
 
    if (argc != 7)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }

    X_str = argv[1];
    ya_str = argv[2];
    yb_str = argv[3];
    ta_str = argv[4];
    tb_str = argv[5];
    c_str = argv[6];

//    X_str = "1000";
//    ya_str = "0.2";
//    yb_str = "1.0";
//    ta_str = "0";
//    tb_str = "0.2";
//    c_str = "0.05";

    //abbeff was calibrated against other software to be accurate at 12 digits, prec=100 for x <= 10^14.
    digits = 12;
    prec = 100;

    arb_set_str(X, X_str, prec);
    arb_set_str(ya, ya_str, prec);
    arb_set_str(yb, yb_str, prec);
    arb_set_str(ta, ta_str, prec);
    arb_set_str(tb, tb_str, prec);
    arb_set_str(c, c_str, prec);
    N = get_N(tb, X, prec);

    arb_set(xa,X);
    arb_add_si(xb, xa, 1, prec);


    arb_zero(z);
    arb_set_d(a, 10000);
    if (arb_lt(X, a) || arb_lt(ta, z) || arb_le(tb, z) || arb_le(ya, z) || arb_le(yb, z) || arb_le(c, z))
        {
            usage_error = 1;
            result = EXIT_FAILURE;
            goto finish;
        }

    if (arb_lt(yb, ya) || arb_lt(tb, ta) )
        {
            usage_error = 1;
            result = EXIT_FAILURE;
            goto finish;
        }


    //situation where yb=ya and tb=ta
    if (arb_equal(yb, ya) && arb_equal(tb, ta))
        {
          flint_printf("\n");
          printf("N = : %ld\n", N);
          flint_printf("\n");
          abbeff_series(out1, X, ya, ta, N, 1, prec);

          flint_printf("(A+B)/B0 eff Re: ");
          arf_printd(arb_midref(acb_realref(out1)), digits);
          flint_printf("\n");
          flint_printf("(A+B)/B0 eff Im: ");
          arf_printd(arb_midref(acb_imagref(out1)), digits);
          flint_printf("\n");
          flint_printf("\n");

          ddxbound_series(out, X, ya, ta, N, 1, prec);
          printf("ddxbound : ");
          arf_printd(arb_midref(out), digits);
          printf("\n");

          ddtbound_series(out, X, ya, ta, N, 1, prec);
          printf("ddtbound : ");
          arf_printd(arb_midref(out), digits);
          printf("\n");
          flint_printf("\n"); 
          goto finish;
        }

    mesh_series(windtot, xa, xb, ya, yb, ta, tb, c, N, 1, prec, digits);
    printf("\n");
    printf("Winding number: \n");
    arf_printd(arb_midref(windtot), digits);
    printf("\n");

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Usage:\n");
        flint_printf("%s X ya yb ta tb c \n\n", argv[0]);
        flint_printf(
    "This script computes the winding number for a 3D 'two plate' Barrier,\n"
    "     placed somewhere at an X > 1000. It computes the function: \n"
    "       H_t(z) (eff) = (Aeff_t(z) + Beff_t(z)) / Beff0_t(z),\n"
    "                  with z=x+iy and ta <= t <= tb. \n"
	"         The script follows a counter clockwise path\n"
    "              along the boundaries of the rectangle: \n"
    "                [X <= x <= X+1] + i[ya <= y <= yb] \n"
    "                starting and ending at (xa + i*ya).\n"
    "                Outputs are: x, y, t and H_t(z).\n"
    "         Note: X, ya, yb, tb, c all must be > 0 and ta >= 0.\n"
    "    If ya=yb and ta=tb, the script will not 'wind' but compute:\n"
    "       N, H_t(z) (eff), ddxbound and ddtbound for X, ya, ta. \n"
    "        Results were calibrated against other software and \n"
	"        should be accurate at 12 digits for all X <= 10^14.\n");
    }
 
    arb_clear(X);
    arb_clear(xa);
    arb_clear(xb);
    arb_clear(ya);
    arb_clear(yb);
    arb_clear(ta);
    arb_clear(tb);
    arb_clear(a);
    arb_clear(c);
    arb_clear(d);
    arb_clear(z);
    arb_clear(re);
    arb_clear(im);
    arb_clear(out);
    arb_clear(windtot);

    acb_clear(out1);
 
    flint_cleanup();

    return result;
}