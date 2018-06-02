/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "acb_poly.h"
#include "acb_mat.h"
#include "flint/profiler.h"


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
 
void
ddxbound_series(arb_t res, const acb_t rx, const acb_t ry, const acb_t rt, slong N, slong n, slong prec)
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
    acb_poly_set_acb(x, rx);
    acb_poly_set_acb(y, ry);
    acb_poly_set_acb(t, rt);

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
ddtbound_series(arb_t res, const acb_t rx, const acb_t ry, const acb_t rt, slong N, slong n, slong prec)
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
    acb_poly_set_acb(x, rx);
    acb_poly_set_acb(y, ry);
    acb_poly_set_acb(t, rt);

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

static void
print_details(slong res, const acb_t ct, acb_t a, acb_t b, acb_t c)

{
    arf_printd(arb_midref(acb_realref(ct)), 20);
    printf(", ");
    arf_printd(arb_midref(acb_realref(a)), 20);
    printf(", ");
    arf_printd(arb_midref(acb_realref(b)), 20);
    printf(", ");
    acb_printn(c, 20, ARB_STR_NO_RADIUS);
    printf("\n");
}

//process a rectangle for each t
void
abbeff_symmetric_rectangle(acb_mat_t ests, const acb_t X, const acb_t y, const acb_t ct, 
	slong num, acb_poly_t finpolyb, acb_poly_t finpolya, acb_t logtn0, acb_t n0acb, slong prt, slong n, slong prec)
{        
    acb_mat_t sarr, thtarr, zarr, afac, bexpo, aexpo, bsums, asums; 
    
    acb_poly_t polya, polyb, polyc, polyd, polys;
    acb_poly_init(polya);
    acb_poly_init(polyb);
    acb_poly_init(polyc);
    acb_poly_init(polyd);
    acb_poly_init(polys);

    acb_t a, b, c, d, one, numacb, vacb;
    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(d);
    acb_init(numacb);
    acb_init(vacb);
    acb_init(one);
    acb_one(one);

    arb_t modabb, minmodabb;
    arb_init(modabb);
    arb_init(minmodabb);

    arb_set_si(minmodabb, 1000);

    slong k, v, result;
    result=0;
    acb_set_ui(numacb, num);
		
		//filling the initial matrices to prep for the walk along the rectangle
		acb_mat_init(thtarr, num, 1);
		acb_mat_init(zarr, num, 1);
		for (v = 0; v < num; v++)
		{
           acb_set_ui(vacb,v); 
		   //thtarr
		   acb_mul_2exp_si(a, one, -1);
		   acb_sub_si(b, numacb, 1, prec);
		   acb_div(b, vacb, b, prec);
		   acb_sub(b, b, a, prec);
		   acb_set(acb_mat_entry(thtarr, v, 0), b);

		   //zarr 
		   acb_neg(a, y); 
		   acb_add_si(a, a, 1, prec);
		   acb_sub_si(b, numacb, 1, prec);
		   acb_div(b, vacb, b, prec);
		   acb_mul_2exp_si(c, a, 1);
		   acb_mul(c, c, b, prec);
		   acb_sub(c, c, a, prec);
		   acb_set(acb_mat_entry(zarr, v, 0), c);
		}

		//run along all four sides of the rectangle

		acb_mat_init(sarr, 4*num-4, 1);
		acb_mat_init(afac, 4*num-4, 1); 
		acb_mat_init(bexpo, 4*num-4, 1);
		acb_mat_init(aexpo, 4*num-4, 1);
		acb_mat_init(bsums, 4*num-4, 1);
		acb_mat_init(asums, 4*num-4, 1);

		//x lower constant
		for (v = 0; v < num-1; v++)
		{ 	
		   //sarr
		   acb_neg(a, y);
		   acb_add_si(a, a, 1, prec);
		   acb_mul_onei(b, X);
		   acb_add(b, b, a, prec);
		   acb_sub(b, b, acb_mat_entry(zarr, v, 0), prec);
		   acb_mul_2exp_si(c, one, -1);
		   acb_mul_onei(c, c);
		   acb_sub(c, b, c, prec);
		   acb_mul_2exp_si(c, c, -1);
		   acb_set(acb_mat_entry(sarr, v, 0), c);

		   //afac
		   acb_mul_2exp_si(a, ct, -2);
		   acb_poly_set_acb(polys, acb_mat_entry(sarr, v, 0));
		   acb_poly_alpha1_series(polya, polys, n, prec);
		   acb_poly_pow_ui(polyb, polya, 2, prec);
		   acb_poly_neg(polys, polys);
		   acb_poly_add_si(polys, polys, 1, prec);
		   acb_poly_alpha1_series(polyc, polys, n, prec);
		   acb_poly_pow_ui(polyd, polyc, 2, prec);
		   acb_poly_sub_series(polyb, polyb, polyd, n, prec);
		   acb_poly_evaluate(b, polyb, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_exp(b, b, prec);
		   acb_poly_set_acb(polys, acb_mat_entry(sarr, v, 0));
		   acb_poly_H01_series(polyb, polys, n, prec);            
		   acb_poly_neg(polys, polys);
		   acb_poly_add_si(polys, polys, 1, prec);
		   acb_poly_H01_series(polyd, polys, n, prec);   
		   acb_poly_div_series(polyb, polyb, polyd, n, prec);
		   acb_poly_evaluate(c, polyb, one, prec);
		   acb_mul(b, b, c, prec);
		   acb_set(acb_mat_entry(afac, v, 0), b);

		   //bexpo
		   acb_mul_2exp_si(a, ct, -1);
		   acb_neg(a, a);
		   acb_poly_evaluate(b, polyc, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_mul_2exp_si(c, acb_mat_entry(zarr,v,0), -1);
		   acb_sub(b, b, c, prec);
		   acb_mul_2exp_si(d, one, -2);
		   acb_mul_onei(d, d);
		   acb_sub(b, b, d, prec);
		   acb_set(acb_mat_entry(bexpo, v, 0), b);

		   //aexpo
		   acb_mul_2exp_si(a, ct, -1);
		   acb_neg(a, a);
		   acb_poly_evaluate(b, polya, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_mul_2exp_si(c, acb_mat_entry(zarr,v,0), -1);
		   acb_add(b, b, c, prec);
		   acb_mul_2exp_si(d, one, -2);
		   acb_mul_onei(d, d);
		   acb_add(b, b, d, prec);
		   acb_set(acb_mat_entry(aexpo, v, 0), b);

		   //bsums
		   acb_mul_2exp_si(d, logtn0, -1);
		   acb_add(a, acb_mat_entry(bexpo, v, 0), d, prec);
		   acb_pow(a, n0acb, a, prec);
		   acb_add(b, acb_mat_entry(bexpo, v, 0), logtn0, prec);
		   acb_poly_evaluate(c, finpolyb, b, prec);
		   acb_mul(b, c, a, prec);
		   acb_set(acb_mat_entry(bsums, v, 0), b);

		   //asums
		   acb_add(a, acb_mat_entry(aexpo, v, 0), d, prec);
		   acb_pow(a, n0acb, a, prec);
		   acb_add(b, acb_mat_entry(aexpo, v, 0), logtn0, prec);
		   acb_poly_evaluate(c, finpolya, b, prec);
		   acb_mul(b, c, a, prec);
		   acb_set(acb_mat_entry(asums, v, 0), b);

		   acb_mul(a, acb_mat_entry(afac, v, 0), acb_mat_entry(asums, v, 0), prec);
		   acb_add(a, a, acb_mat_entry(bsums, v, 0), prec);
		   acb_set(acb_mat_entry(ests, v, 0), a);
           acb_abs(modabb, a, prec);
           if (arb_lt(modabb, minmodabb) )
           {
                 arb_set(minmodabb, modabb);
           } 

		   if (prt==1)
		   {      
              acb_add(a, y, acb_mat_entry(zarr, v, 0), prec);
              acb_mul_2exp_si(d, one, -1);  
              acb_sub(b, X, d, prec);
              acb_set(c, acb_mat_entry(ests, v, 0));
              print_details(result, ct, a, b, c);
		   }

		}

		//y upper constant
		for (v = 0; v < num-1; v++)
		{ 
		   k=num+v-1;

		   //sarr
		   acb_neg(a, y);
		   acb_add_si(a, a, 1, prec);
		   acb_add(b, X, acb_mat_entry(thtarr, v, 0), prec);
		   acb_mul_onei(b, b);
		   acb_add(b, b, a, prec);
		   acb_sub(b, b, a, prec);
		   acb_mul_2exp_si(b, b, -1);
		   acb_set(acb_mat_entry(sarr, k, 0), b);

		   //afac
		   acb_mul_2exp_si(a, ct, -2);
		   acb_poly_set_acb(polys, acb_mat_entry(sarr, k,0));
		   acb_poly_alpha1_series(polya, polys, n, prec);
		   acb_poly_pow_ui(polyb, polya, 2, prec);
		   acb_poly_neg(polys, polys);
		   acb_poly_add_si(polys, polys, 1, prec);
		   acb_poly_alpha1_series(polyc, polys, n, prec);
		   acb_poly_pow_ui(polyd, polyc, 2, prec);
		   acb_poly_sub_series(polyb, polyb, polyd, n, prec);
		   acb_poly_evaluate(b, polyb, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_exp(b, b, prec);
		   acb_poly_set_acb(polys, acb_mat_entry(sarr, k,0));
		   acb_poly_H01_series(polyb, polys, n, prec);            
		   acb_poly_neg(polys, polys);
		   acb_poly_add_si(polys, polys, 1, prec);
		   acb_poly_H01_series(polyd, polys, n, prec);   
		   acb_poly_div_series(polyb, polyb, polyd, n, prec);
		   acb_poly_evaluate(c, polyb, one, prec);
		   acb_mul(b, b, c, prec);
		   acb_set(acb_mat_entry(afac, k, 0), b);

		   //bexpo
		   acb_mul_2exp_si(a, ct, -1);
		   acb_neg(a, a);
		   acb_poly_evaluate(b, polyc, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_mul_2exp_si(c, acb_mat_entry(thtarr,v, 0), -1);
		   acb_mul_onei(c, c);
		   acb_neg(d, y);
		   acb_add_si(d, d, 1, prec);
		   acb_mul_2exp_si(d, d, -1);
		   acb_add(b, b, c, prec);
		   acb_sub(b, b, d, prec);
		   acb_set(acb_mat_entry(bexpo, k, 0), b);

		   //aexpo
		   acb_mul_2exp_si(a, ct, -1);
		   acb_neg(a, a);
		   acb_poly_evaluate(b, polya, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_mul_2exp_si(c, acb_mat_entry(thtarr,v, 0), -1);
		   acb_mul_onei(c, c);
		   acb_neg(d, y);
		   acb_add_si(d, d, 1, prec);
		   acb_mul_2exp_si(d, d, -1);
		   acb_sub(b, b, c, prec);
		   acb_add(b, b, d, prec);
		   acb_set(acb_mat_entry(aexpo, k, 0), b);

		   //bsums
		   acb_mul_2exp_si(d, logtn0, -1);
		   acb_add(a, acb_mat_entry(bexpo, k, 0), d, prec);
		   acb_pow(a, n0acb, a, prec);
		   acb_add(b, acb_mat_entry(bexpo, k, 0), logtn0, prec);
		   acb_poly_evaluate(c, finpolyb, b, prec);
		   acb_mul(b, c, a, prec);
		   acb_set(acb_mat_entry(bsums, k, 0), b);

		   //asums
		   acb_add(a, acb_mat_entry(aexpo, k, 0), d, prec);
		   acb_pow(a, n0acb, a, prec);
		   acb_add(b, acb_mat_entry(aexpo, k, 0), logtn0, prec);
		   acb_poly_evaluate(c, finpolya, b, prec);
		   acb_mul(b, c, a, prec);
		   acb_set(acb_mat_entry(asums, k, 0), b);

		   acb_mul(a, acb_mat_entry(afac, k, 0), acb_mat_entry(asums, k, 0), prec);
		   acb_add(a, a, acb_mat_entry(bsums, k, 0), prec);
		   acb_set(acb_mat_entry(ests, k, 0), a);
           acb_abs(modabb, a, prec);
           if (arb_lt(modabb, minmodabb) )
           {
                 arb_set(minmodabb, modabb);
           } 

		   if (prt==1)
			{      
              acb_add_si(a, y, 1, prec);
              acb_sub(a, a, y, prec);
              acb_add(b, X, acb_mat_entry(thtarr, v, 0), prec);
              acb_set(c, acb_mat_entry(ests, k, 0));
              print_details(result, ct, a, b, c);
			}
		}

		//x upper constant and output to be attached in reverse order
		for (v = num-1; v > 0; v--)
		{ 
		   k=3*num-v-3;

		   //sarr
		   acb_neg(a, y);
		   acb_add_si(a, a, 1, prec);
		   acb_mul_onei(b, X);
		   acb_add(b, b, a, prec);
		   acb_sub(b, b, acb_mat_entry(zarr, v, 0), prec);
		   acb_mul_2exp_si(c, one, -1);
		   acb_mul_onei(c, c);
		   acb_add(c, b, c, prec);
		   acb_mul_2exp_si(c, c, -1);
		   acb_set(acb_mat_entry(sarr, k, 0), c);

		   //afac
		   acb_mul_2exp_si(a, ct, -2);
		   acb_poly_set_acb(polys, acb_mat_entry(sarr, k,0));
		   acb_poly_alpha1_series(polya, polys, n, prec);
		   acb_poly_pow_ui(polyb, polya, 2, prec);
		   acb_poly_neg(polys, polys);
		   acb_poly_add_si(polys, polys, 1, prec);
		   acb_poly_alpha1_series(polyc, polys, n, prec);
		   acb_poly_pow_ui(polyd, polyc, 2, prec);
		   acb_poly_sub_series(polyb, polyb, polyd, n, prec);
		   acb_poly_evaluate(b, polyb, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_exp(b, b, prec);
		   acb_poly_set_acb(polys, acb_mat_entry(sarr, k,0));
		   acb_poly_H01_series(polyb, polys, n, prec);            
		   acb_poly_neg(polys, polys);
		   acb_poly_add_si(polys, polys, 1, prec);
		   acb_poly_H01_series(polyd, polys, n, prec);   
		   acb_poly_div_series(polyb, polyb, polyd, n, prec);
		   acb_poly_evaluate(c, polyb, one, prec);
		   acb_mul(b, b, c, prec);
		   acb_set(acb_mat_entry(afac, k, 0), b);

		   //bexpo
		   acb_mul_2exp_si(a, ct, -1);
		   acb_neg(a, a);
		   acb_poly_evaluate(b, polyc, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_mul_2exp_si(c, acb_mat_entry(zarr,v,0), -1);
		   acb_sub(b, b, c, prec);
		   acb_mul_2exp_si(d, one, -2);
		   acb_mul_onei(d, d);
		   acb_add(b, b, d, prec);
		   acb_set(acb_mat_entry(bexpo, k, 0), b);

		   //aexpo
		   acb_mul_2exp_si(a, ct, -1);
		   acb_neg(a, a);
		   acb_poly_evaluate(b, polya, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_mul_2exp_si(c, acb_mat_entry(zarr, v, 0), -1);
		   acb_add(b, b, c, prec);
		   acb_mul_2exp_si(d, one, -2);
		   acb_mul_onei(d, d);
		   acb_sub(b, b, d, prec);
		   acb_set(acb_mat_entry(aexpo, k, 0), b);

		   //bsums
		   acb_mul_2exp_si(d, logtn0, -1);
		   acb_add(a, acb_mat_entry(bexpo, k, 0), d, prec);
		   acb_pow(a, n0acb, a, prec);
		   acb_add(b, acb_mat_entry(bexpo, k, 0), logtn0, prec);
		   acb_poly_evaluate(c, finpolyb, b, prec);
		   acb_mul(b, c, a, prec);
		   acb_set(acb_mat_entry(bsums, k, 0), b);

		   //asums
		   acb_add(a, acb_mat_entry(aexpo, k, 0), d, prec);
		   acb_pow(a, n0acb, a, prec);
		   acb_add(b, acb_mat_entry(aexpo, k, 0), logtn0, prec);
		   acb_poly_evaluate(c, finpolya, b, prec);
		   acb_mul(b, c, a, prec);
		   acb_set(acb_mat_entry(asums, k, 0), b);

		   acb_mul(a, acb_mat_entry(afac, k, 0), acb_mat_entry(asums, k, 0), prec);
		   acb_add(a, a, acb_mat_entry(bsums, k, 0), prec);
		   acb_set(acb_mat_entry(ests, k, 0), a);
           acb_abs(modabb, a, prec);
           if (arb_lt(modabb, minmodabb) )
           {
                 arb_set(minmodabb, modabb);
           } 

		   if (prt==1)
			{      
              acb_add(a, y, acb_mat_entry(zarr, v, 0), prec);
              acb_mul_2exp_si(d, one, -1);  
              acb_add(b, X, d, prec);
              acb_set(c, acb_mat_entry(ests, k, 0));
              print_details(result, ct, a, b, c);
			}
		}

		//y lower constant and output to be attached in reverse order
		for (v = num-1; v > 0; v--)
		{ 
		   k=4*num-v-4;

		   //sarr
		   acb_neg(a, y);
		   acb_add_si(a, a, 1, prec);
		   acb_add(b, X, acb_mat_entry(thtarr, v, 0), prec);
		   acb_mul_onei(b, b);
		   acb_add(b, b, a, prec);
		   acb_add(b, b, a, prec);
		   acb_mul_2exp_si(b, b, -1);
		   acb_set(acb_mat_entry(sarr, k, 0), b);

		   //afac
		   acb_mul_2exp_si(a, ct, -2);
		   acb_poly_set_acb(polys, acb_mat_entry(sarr, k,0));
		   acb_poly_alpha1_series(polya, polys, n, prec);
		   acb_poly_pow_ui(polyb, polya, 2, prec);
		   acb_poly_neg(polys, polys);
		   acb_poly_add_si(polys, polys, 1, prec);
		   acb_poly_alpha1_series(polyc, polys, n, prec);
		   acb_poly_pow_ui(polyd, polyc, 2, prec);
		   acb_poly_sub_series(polyb, polyb, polyd, n, prec);
		   acb_poly_evaluate(b, polyb, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_exp(b, b, prec);
		   acb_poly_set_acb(polys, acb_mat_entry(sarr, k,0));
		   acb_poly_H01_series(polyb, polys, n, prec);            
		   acb_poly_neg(polys, polys);
		   acb_poly_add_si(polys, polys, 1, prec);
		   acb_poly_H01_series(polyd, polys, n, prec);   
		   acb_poly_div_series(polyb, polyb, polyd, n, prec);
		   acb_poly_evaluate(c, polyb, one, prec);
		   acb_mul(b, b, c, prec);
		   acb_set(acb_mat_entry(afac, k, 0), b);

		   //bexpo
		   acb_mul_2exp_si(a, ct, -1);
		   acb_neg(a, a);
		   acb_poly_evaluate(b, polyc, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_mul_2exp_si(c, acb_mat_entry(thtarr,v, 0), -1);
		   acb_mul_onei(c, c);
		   acb_neg(d, y);
		   acb_add_si(d, d, 1, prec);
		   acb_mul_2exp_si(d, d, -1);
		   acb_add(b, b, c, prec);
		   acb_add(b, b, d, prec);
		   acb_set(acb_mat_entry(bexpo, k, 0), b);

		   //aexpo
		   acb_mul_2exp_si(a, ct, -1);
		   acb_neg(a, a);
		   acb_poly_evaluate(b, polya, one, prec);
		   acb_mul(b, b, a, prec);
		   acb_mul_2exp_si(c, acb_mat_entry(thtarr,v, 0), -1);
		   acb_mul_onei(c, c);
		   acb_neg(d, y);
		   acb_add_si(d, d, 1, prec);
		   acb_mul_2exp_si(d, d, -1);
		   acb_sub(b, b, c, prec);
		   acb_sub(b, b, d, prec);
		   acb_set(acb_mat_entry(aexpo, k, 0), b);

		   //bsums
		   acb_mul_2exp_si(d, logtn0, -1);
		   acb_add(a, acb_mat_entry(bexpo, k, 0), d, prec);
		   acb_pow(a, n0acb, a, prec);
		   acb_add(b, acb_mat_entry(bexpo, k, 0), logtn0, prec);
		   acb_poly_evaluate(c, finpolyb, b, prec);
		   acb_mul(b, c, a, prec);
		   acb_set(acb_mat_entry(bsums, k, 0), b);

		   //asums
		   acb_add(a, acb_mat_entry(aexpo, k, 0), d, prec);
		   acb_pow(a, n0acb, a, prec);
		   acb_add(b, acb_mat_entry(aexpo, k, 0), logtn0, prec);
		   acb_poly_evaluate(c, finpolya, b, prec);
		   acb_mul(b, c, a, prec);
		   acb_set(acb_mat_entry(asums, k, 0), b);

		   acb_mul(a, acb_mat_entry(afac, k, 0), acb_mat_entry(asums, k, 0), prec);
		   acb_add(a, a, acb_mat_entry(bsums, k, 0), prec);
		   acb_set(acb_mat_entry(ests, k, 0), a);
           acb_abs(modabb, a, prec);
           if (arb_lt(modabb, minmodabb) )
           {
                 arb_set(minmodabb, modabb);
           } 

		   if (prt==1)
			{      
              acb_add_si(a, y,-1, prec);
              acb_add(a, a, y, prec);
              acb_add(b, X, acb_mat_entry(thtarr, v, 0), prec);
              acb_set(c, acb_mat_entry(ests, k, 0));
              print_details(result, ct, a, b, c);
			}
		}

    acb_set_arb(a, minmodabb);
    acb_set(acb_mat_entry(ests, 4*num-4, 0), a);

    acb_mat_clear(thtarr);
    acb_mat_clear(sarr);
    acb_mat_clear(zarr);
    acb_mat_clear(afac);
    acb_mat_clear(bexpo);
    acb_mat_clear(aexpo);
    acb_mat_clear(bsums);
    acb_mat_clear(asums);

    acb_poly_clear(polya);
    acb_poly_clear(polyb);
    acb_poly_clear(polyc);
    acb_poly_clear(polyd);
    acb_poly_clear(polys);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(d);
    acb_clear(one);
    acb_clear(numacb);
    acb_clear(vacb);

    arb_clear(modabb);
    arb_clear(minmodabb);

}

//procedure to transform the matrix into a poly
static void
mattopoly(acb_poly_t polyres, const arb_t tval, const slong expterms, const slong taylorterms, 
           acb_mat_t mattmp, slong prec)
{

    arb_t ab;
    acb_t a, b;

    arb_init(ab);
    acb_init(a);
    acb_init(b);

    slong e, t;

    acb_poly_zero(polyres);

    for (e = 0; e < expterms; e++)
    {
        for (t = 0; t < taylorterms; t++)
        {
            arb_pow_ui(ab, tval, t, prec);
            acb_set_arb(a, ab);
            acb_mul(b, a, acb_mat_entry(mattmp, e, t), prec);
            acb_poly_get_coeff_acb(a, polyres, e);
            acb_add(b, b, a, prec);
            acb_poly_set_coeff_acb(polyres, e, b);
        }
    }

    arb_clear(ab);
    acb_clear(a);
    acb_clear(b);
}

void 
abbeff_t_loop(slong res, const acb_t cX, const acb_t cy, const arb_t ts, const arb_t te, 
           const slong H, const slong taylorterms, const slong expterms, acb_mat_t finalmatb, acb_mat_t finalmata,
           slong prt, slong n, slong prec)

{
    acb_mat_t ests;

    arb_t pi, ab, harb, Harb, n0, t, modabb, minmodabb, windnum, windtot, dxabb, dtabb;
    arb_init(ab);
    arb_init(harb);
    arb_init(Harb);
    arb_init(n0);
    arb_init(modabb);
    arb_init(minmodabb);
    arb_init(windnum);
    arb_init(windtot);
    arb_init(dxabb);
    arb_init(dtabb);
    arb_init(t);
    arb_init(pi);
    arb_const_pi(pi, prec);

    acb_t a, b, X, y, ct, one, hacb, logtn0, n0acb, argdiv;
    acb_init(a);
    acb_init(b);
    acb_init(X);
    acb_init(y);
    acb_init(ct);
    acb_init(argdiv);
    acb_init(hacb);
    acb_init(logtn0);
    acb_init(n0acb);
    acb_init(one);
    acb_one(one);

    acb_poly_t finpolyb, finpolya;
    acb_poly_init(finpolyb);
    acb_poly_init(finpolya);
   
    arb_set_si(Harb, H);

    printf("\n");
    printf("Processing the barrier for X= ");
    arf_printd(arb_midref(acb_realref(cX)), 20);
    printf("...");
    acb_add_si(a, cX, 1, prec);
    arf_printd(arb_midref(acb_realref(a)), 20);
    printf(" (N = %ld)", H);
    printf(", y0 = ");
    arf_printd(arb_midref(acb_realref(cy)), 10);
    printf("...1 ");
    printf(", t = ");
    arf_printd(arb_midref(ts), 10);
    printf("...");
    arf_printd(arb_midref(te), 10);
    printf("\n");
    printf("\n");

    slong idx, count, rectmesh, sidemesh, num;

    //change X and y to midpoints
    acb_mul_2exp_si(a, one, -1);  
    acb_add(X, cX, a, prec); 
    acb_add_si(a, cy, 1, prec);
    acb_mul_2exp_si(y, a, -1);

    //filling all the initial matrices (once off)

    //n0
    arb_mul_2exp_si(n0, Harb, -1);
    acb_set_arb(n0acb, n0);

    //perform the t-loop over all the X..X+1, y0..1 rectangles
    arb_set(t, ts);
    arb_zero(windtot);
    count = 0;
    acb_set_si(ct,0);

    while(arb_le(t, te))
    {
        count=count+1;
        acb_set_arb(ct, t);

        //establish the ddx and ddt bounds for each t
        ddxbound_series(dxabb, cX, cy, ct, H, 1, prec);
        ddtbound_series(dtabb, cX, cy, ct, H, 1, prec);

        //logtn0
        acb_mul_2exp_si(a, ct, -1);
        acb_log(b, n0acb, prec);
        acb_mul(logtn0, b, a, prec);

        mattopoly(finpolyb, t, expterms, taylorterms, finalmatb, prec);
        mattopoly(finpolya, t, expterms, taylorterms, finalmata, prec);

        arb_ceil(ab, dxabb, prec);
        sidemesh = arf_get_d(arb_midref(ab), ARF_RND_NEAR);
        rectmesh = 4 * sidemesh - 4;
        num = sidemesh;

        acb_mat_init(ests, 4*num-3, 1);

        //evaluate the rectangle
        abbeff_symmetric_rectangle(ests, X, y, ct, num, finpolyb, finpolya, logtn0, n0acb, prt, n, prec);
 
        //calculate and print the winding number for this x,y rectangle
        arb_zero(windnum);
        for (idx = 0; idx < rectmesh-1; idx++)
        { 
            acb_div(argdiv, acb_mat_entry(ests, idx, 0), acb_mat_entry(ests, idx+1, 0), prec);
            acb_arg(ab, argdiv, prec);
            arb_add(windnum, windnum, ab, prec);
        }   

        acb_div(argdiv, acb_mat_entry(ests, rectmesh-1, 0), acb_mat_entry(ests, 0, 0), prec);
        acb_arg(ab, argdiv, prec);
        arb_add(windnum, windnum, ab, prec);
        arb_div(windnum, windnum, pi, prec);
        arb_mul_2exp_si(windnum, windnum, -1);

        arb_add(windtot, windtot, windnum, prec);

        acb_get_real(minmodabb, acb_mat_entry(ests, 4*num-4, 0));

        //print rectangle summary
        printf("Rectangle(%ld) : ", count);
        arf_printd(arb_midref(t), 20);
        printf(", ");
        arf_printd(arb_midref(dtabb), 20);
        printf(", ");
        arf_printd(arb_midref(dxabb), 20);
        printf(", ");
        arf_printd(arb_midref(windnum), 10);
        printf(", ");
        arf_printd(arb_midref(minmodabb), 20);
        printf(", ");
        printf("%ld\n", rectmesh);

        //establish the next t
        arb_div(ab, minmodabb, dtabb, prec);
        arb_mul_2exp_si(ab, ab, -1);
        arb_add(t, t, ab, prec);

        //prevent drift in the error term of t
        arb_get_mid_arb(t,t);

        acb_mat_clear(ests);
    }

    flint_printf("\n");
    printf("Overall winding number: %f \n", arf_get_d(arb_midref(windtot), ARF_RND_NEAR));
    printf("\n");

    arb_clear(ab);
    arb_clear(harb);
    arb_clear(Harb);
    arb_clear(n0);
    arb_clear(modabb);
    arb_clear(minmodabb);
    arb_clear(windnum);
    arb_clear(windtot);
    arb_clear(dxabb);
    arb_clear(dtabb);
    arb_clear(t);
    arb_clear(pi);

    acb_clear(a);
    acb_clear(b);
    acb_clear(X);
    acb_clear(y);
    acb_clear(ct);
    acb_clear(argdiv);
    acb_clear(hacb);
    acb_clear(logtn0);
    acb_clear(n0acb);
    acb_clear(one);

    acb_poly_clear(finpolya);
    acb_poly_clear(finpolyb);
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

const char* getfield(char* line, int num)
{
    const char* tok;
    for (tok = strtok(line, ",");
            tok && *tok;
            tok = strtok(NULL, ",\n"))
    {
        if (!--num)
            return tok;
    }
    return NULL;
}
 
int main(int argc, char *argv[])
{
    arb_t X, y0, ts, te, re, im;
    acb_t a, cX, cy0;

    const char *ts_str, *te_str, *expterms_str, *taylorterms_str, *prt_str;
    slong e, t, N, prec, expterms, taylorterms, prt, res, linesize;
    int result = EXIT_SUCCESS;
    res=0;

    arb_init(X);
    arb_init(y0);
    arb_init(ts);
    arb_init(te);
    arb_init(re);
    arb_init(im);

    acb_init(a);
    acb_init(cX);
    acb_init(cy0);

    //precision
    prec = 128;

    linesize=50000;
    char line[linesize];

    if (argc != 7)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    ts_str = argv[1];
    te_str = argv[2];
    expterms_str = argv[3];
    taylorterms_str = argv[4];
    prt_str = argv[5];

    arb_set_str(ts, ts_str, prec);
    arb_set_str(te, te_str, prec);

    expterms = atol(expterms_str);
    taylorterms = atol(taylorterms_str);
    prt = atol(prt_str);

    acb_mat_t finalmatb, finalmata;
    acb_mat_init(finalmatb, expterms, taylorterms);
    acb_mat_init(finalmata, expterms, taylorterms);



    FILE *f = fopen(argv[6], "r");

    //recover X, y0
    fgets(line, linesize, f);
    char* str = strdup(line);
    arb_set_str(X, getfield(str, 1), prec);
    acb_set_arb(cX, X);
    str = strdup(line);
    arb_set_str(y0, getfield(str, 2), prec);
    acb_set_arb(cy0, y0);
    free(str);

    //fill finalmatb
    printf("Filling finalmatb\n");
    e=1;
    while (e <= expterms)
    {
        fgets(line, linesize, f);
        t=1;
        while (t <= 2*taylorterms)
        {
            char* str = strdup(line);
            arb_set_str(re, getfield(str, t), prec);
            str = strdup(line);
            arb_set_str(im, getfield(str, t+1), prec);
            acb_set_arb_arb(a, re, im);
            acb_set(acb_mat_entry(finalmatb, e-1, (t-1)/2), a);
            free(str);
            t=t+2;
        }
        e=e+1;
    }

    //scroll to ***
    while(!(strncmp (line,"***",3) == 0))
    {
        fgets(line, linesize, f);
        e=e+1;
    }
    //finalmatb
	printf("Filling finalmata\n");
    e=1;
    while (e <= expterms)
    {
        fgets(line, linesize, f);
        t=1;
        while (t <= 2*taylorterms)
        {
            char* str = strdup(line);
            arb_set_str(re, getfield(str, t), prec);
            str = strdup(line);
            arb_set_str(im, getfield(str, t+1), prec);
            acb_set_arb_arb(a, re, im);
            acb_set(acb_mat_entry(finalmata, e-1, (t-1)/2), a);
            free(str);
            t=t+2;
        }
        e=e+1;
    }

    fclose(f);

    N = get_N(ts, X, prec);

TIMEIT_ONCE_START

    abbeff_t_loop(res, cX, cy0, ts, te, N, taylorterms, expterms, finalmatb, finalmata, prt, 1, prec);

TIMEIT_ONCE_STOP

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Required inputs:\n");
        flint_printf("%s ts te expterms taylorterms Prt \n\n", argv[0]);
        flint_printf(
    "This script computes the winding number for a '3D-Barrier',\n"
    "that runs along rectangle: [X <= x <= X+1] + i[y0 <= y <= 1],\n"
    "and along: [ts <= t <= te]. It takes X, y0  from the input file.\n"
    "The two terms are the number of Taylor terms used to approximate\n"
    "the Exp-terms. With Prt the output can be controlled:\n"
    " 0 = print rectangle summary only, 1 = print full details.\n"
    "The ACB-precision has been fixed at 128 after outcomes were\n" 
	"calibrated against other software. Outputs should be accurate \n"
	"to at least 12 digits for all X <= 10^15.\n");
    }

    arb_clear(X);
    arb_clear(y0);
    arb_clear(ts);
    arb_clear(te);
    arb_clear(re);
    arb_clear(im);

    acb_clear(a);
    acb_clear(cX);
    acb_clear(cy0);

    acb_mat_clear(finalmatb);
    acb_mat_clear(finalmata);
 
    flint_cleanup();

    return result;
}
