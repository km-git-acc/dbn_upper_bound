/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
#include "acb_poly.h"
#include "acb_mat.h"

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
abbeff_multieval_xseries(slong res, const acb_t X, const acb_t y, const acb_t t, 
           slong N, slong num, slong H, slong expterms, slong npre, slong digits, slong n, slong prec)
{

    acb_mat_t sarr, afac, n0hmat, n0sumsbmat, n0sumsamat, ybexpo, yaexpo, ybexpolarge, yaexpolarge;
    acb_mat_t n0thtbmat, n0thtamat, bsumarr1, bsumarr2, bsumarr3, bsums, asumarr1, asumarr2, asumarr3, asums, ests;
    acb_mat_t n0, thtarr, n0htmat, n0hlogmat;

    acb_t a, b, c, d, sumA, sumB, sumA1, sumB1, one;
    acb_t expoacb, npreacb, numacb, hacb, Hacb, idxacb, nacb, vacb;

    acb_poly_t polya, polyb, polyc, polyd, polys;

    acb_poly_init(polya);
    acb_poly_init(polyb);
    acb_poly_init(polyc);
    acb_poly_init(polyd);
    acb_poly_init(polys);

    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(d);
    acb_init(sumA);
    acb_init(sumB);
    acb_init(sumA1);
    acb_init(sumB1);
    acb_init(expoacb);
    acb_init(npreacb);
    acb_init(numacb);
    acb_init(hacb);
    acb_init(Hacb);
    acb_init(idxacb);
    acb_init(nacb);
    acb_init(vacb);
    acb_init(one);
    acb_one(one);
   
    acb_set_ui(npreacb,npre);
    acb_set_ui(numacb,num);
    acb_set_ui(Hacb,H);

    slong h, idx, n1, v, w, numn0, expo;
    
    numn0 = round((N-npre)/H);
	
    acb_mat_init(n0, numn0, 1);
    for (v = 0; v < numn0; v++)
    {
       //n0
       acb_set_ui(vacb,v); 
       acb_mul(a, vacb, Hacb, prec);
       acb_add(a, a, npreacb, prec);
       acb_set(acb_mat_entry(n0, v, 0), a);
    }

    acb_mat_init(thtarr, num, 1);
    acb_mat_init(sarr, num, 1);
    acb_mat_init(afac, num, 1);
    acb_mat_init(ybexpo, num, 1);
    acb_mat_init(yaexpo, num, 1);
    acb_mat_init(ybexpolarge, num, 1);
    acb_mat_init(yaexpolarge, num, 1);
    for (v = 0; v < num; v++)
    {
      //thtarr
       acb_set_ui(vacb,v); 
       acb_neg(a, y);
       acb_add_si(a, a, 1, prec);
       acb_mul(a, a, vacb, prec);
       acb_sub_si(b, numacb, 1, prec);
       acb_div(a, a, b, prec);
       acb_set(acb_mat_entry(thtarr, v, 0), a);

       //sarr
       acb_set_ui(vacb,v); 
       acb_neg(a, y);
       acb_add_si(a, a, 1, prec);
       acb_mul_onei(b, X);
       acb_add(b, b, a, prec);
       acb_sub(b, b, acb_mat_entry(thtarr, v, 0), prec);
       acb_mul_2exp_si(b, b, -1);
       acb_set(acb_mat_entry(sarr, v, 0), b);

       //afac
       acb_mul_2exp_si(a, t, -2);
       acb_poly_set_acb(polys, acb_mat_entry(sarr,v,0));
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
       acb_poly_set_acb(polys, acb_mat_entry(sarr,v,0));
       acb_poly_H01_series(polyb, polys, n, prec);            
       acb_poly_neg(polys, polys);
       acb_poly_add_si(polys, polys, 1, prec);
       acb_poly_H01_series(polyd, polys, n, prec);   
       acb_poly_div_series(polyb, polyb, polyd, n, prec);
       acb_poly_evaluate(c, polyb, one, prec);
       acb_mul(b, b, c, prec);
       acb_set(acb_mat_entry(afac, v, 0), b);

       //ybexpo
       acb_add_si(a, y, 1, prec);
       acb_mul_2exp_si(a, a, -1);
       acb_neg(a, a);
       acb_poly_evaluate(b, polyc, one, prec);
       acb_mul_2exp_si(b, b, -1);
       acb_mul(b, b, t, prec);
       acb_mul_2exp_si(c, acb_mat_entry(thtarr,v,0), -1);
       acb_sub(a, a, b, prec);
       acb_sub(a, a, c, prec);
       acb_set(acb_mat_entry(ybexpo, v, 0), a);

       //ybexpolarge
       acb_mul_onei(d, X);
       acb_mul_2exp_si(d, d, -1);
       acb_add(a, a, d, prec);
       acb_set(acb_mat_entry(ybexpolarge, v, 0), a);

       //yaexpo
       acb_add_si(a, y, -1, prec);
       acb_mul_2exp_si(a, a, -1);
       acb_poly_evaluate(b, polya, one, prec);
       acb_mul_2exp_si(b, b, -1);
       acb_mul(b, b, t, prec);
       acb_sub(a, a, b, prec);
       acb_add(a, a, c, prec);
       acb_set(acb_mat_entry(yaexpo, v, 0), a);

       //yaexpolarge
       acb_sub(a, a, d, prec);
       acb_set(acb_mat_entry(yaexpolarge, v, 0), a);

    }

    acb_mat_init(n0hmat, numn0, H);
    acb_mat_init(n0htmat, numn0, H);
    acb_mat_init(n0hlogmat, numn0, H);
    for (v = 0; v < numn0; v++)
    {
        acb_set_ui(vacb,v);
        for (h = 0; h < H; h++)
        {
            //n0hmat
            acb_set_ui(hacb,h);
            acb_add_si(hacb, hacb, 1, prec);
            acb_add(a, acb_mat_entry(n0, v, 0), hacb, prec);
            acb_mul_onei(b, X);
            acb_mul_2exp_si(b, b, -1);
            acb_pow(c, a, b, prec);
            acb_set(acb_mat_entry(n0hmat, v, h), c);

            //n0htmat
            acb_div(a, hacb, acb_mat_entry(n0, v, 0), prec);
            acb_add_si(a, a, 1, prec);
            acb_log(a, a, prec);
            acb_pow_ui(a, a, 2, prec);
            acb_mul_2exp_si(b, t, -2);
            acb_mul(b, a, b, prec);
            acb_exp(b, b, prec);
            acb_set(acb_mat_entry(n0htmat, v, h), b);

            //n0hlogmat
            acb_div(a, hacb, acb_mat_entry(n0, v, 0), prec);
            acb_add_si(a, a, 1, prec);
            acb_log(a, a, prec);
            acb_set(acb_mat_entry(n0hlogmat, v, h), a);
        }
    }

    acb_mat_init(n0sumsbmat, numn0, expterms);
    acb_mat_init(n0sumsamat, numn0, expterms);
    for (v = 0; v < numn0; v++)
    {
        acb_set_ui(vacb,v);
        for (expo = 0; expo < expterms; expo++)
        {
            acb_set_ui(expoacb,expo);
            acb_zero(sumA);
            acb_zero(sumB);
            for (h = 0; h < H; h++)
            {
                acb_pow(a, acb_mat_entry(n0hlogmat, v, h), expoacb, prec);
                acb_mul(a, a, acb_mat_entry(n0htmat, v, h), prec);   
                acb_mul(b, a, acb_mat_entry(n0hmat, v, h), prec);  
                acb_div(c, a, acb_mat_entry(n0hmat, v, h), prec);
                acb_add(sumB, sumB, b, prec);
                acb_add(sumA, sumA, c, prec);
            }
            //n0sumsbmat
            acb_set(acb_mat_entry(n0sumsbmat, v, expo), sumB);

            //n0sumsamat
            acb_set(acb_mat_entry(n0sumsamat, v, expo), sumA);
        }
    }

    acb_mat_init(n0thtbmat, numn0, num);
    acb_mat_init(n0thtamat, numn0, num);
    for (w = 0; w < numn0; w++)
    {
        acb_set_ui(vacb,v);
        for (v = 0; v < num; v++)
        {
            //n0thtbmat
            acb_log(a, acb_mat_entry(n0, w, 0), prec);
            acb_mul_2exp_si(b, t, -1);
            acb_mul(b, b, a, prec);
            acb_add(c, b, acb_mat_entry(ybexpo, v, 0), prec);
            acb_set(acb_mat_entry(n0thtbmat, w, v), c);

            //n0thtamat
            acb_add(c, b, acb_mat_entry(yaexpo, v, 0), prec);
            acb_set(acb_mat_entry(n0thtamat, w, v), c);
        }
    }

    acb_mat_init(bsumarr1, num, 1);
    acb_mat_init(asumarr1, num, 1);
    acb_mat_init(bsumarr2, num, 1);
    acb_mat_init(asumarr2, num, 1);
    acb_mat_init(bsumarr3, num, 1);
    acb_mat_init(asumarr3, num, 1);
    for (v = 0; v < num; v++)
    {
        acb_zero(sumA);
        acb_zero(sumB);
        for (n1 = 1; n1 <= npre; n1++)
        {
            acb_set_si(nacb, n1);
            acb_log(b, nacb, prec);
            acb_mul_2exp_si(b, b, -2);
            acb_mul(b, b, t, prec);
            acb_add(c, b, acb_mat_entry(ybexpolarge, v, 0), prec);
            acb_pow(a, nacb, c, prec);
            acb_add(sumB, sumB, a, prec);
            acb_add(c, b, acb_mat_entry(yaexpolarge, v, 0), prec);
            acb_pow(a, nacb, c, prec);
            acb_add(sumA, sumA, a, prec);
        }
        //bsumarr1
        acb_set(acb_mat_entry(bsumarr1, v, 0), sumB);

        //asumarr1
        acb_set(acb_mat_entry(asumarr1, v, 0), sumA);

        acb_zero(sumA);
        acb_zero(sumB);
        for (n1 = 0; n1 < numn0 ; n1++)
        {
            acb_zero(sumA1);
            acb_zero(sumB1);
            for (idx = 0; idx < expterms ; idx++)
            {
                acb_set_si(idxacb, idx);
                acb_pow(a, acb_mat_entry(n0thtbmat, n1, v), idxacb, prec);
                acb_add_si(d, idxacb, 1, prec); 
                acb_gamma(b, d, prec);
                acb_div(a, a, b, prec);
                acb_mul(a, a, acb_mat_entry(n0sumsbmat, n1, idx), prec);
                acb_add(sumB1, sumB1, a, prec);
                acb_pow(a, acb_mat_entry(n0thtamat, n1, v), idxacb, prec); 
                acb_div(a, a, b, prec);
                acb_mul(a, a, acb_mat_entry(n0sumsamat, n1, idx), prec);
                acb_add(sumA1, sumA1, a, prec);
            }
            acb_set(d, acb_mat_entry(n0, n1, 0));
            acb_log(a, d, prec);
            acb_mul_2exp_si(b, t, -2);
            acb_mul(b, a, b, prec);
            acb_sub(c, acb_mat_entry(n0thtbmat, n1, v), b, prec);
            acb_pow(c, d, c, prec);
            acb_mul(c, c, sumB1, prec);
            acb_add(sumB, sumB, c, prec);
            acb_sub(c, acb_mat_entry(n0thtamat, n1, v), b, prec);
            acb_pow(c, d, c, prec);
            acb_mul(c, c, sumA1, prec);
            acb_add(sumA, sumA, c, prec);
        }
        //bsumarr2
        acb_set(acb_mat_entry(bsumarr2, v, 0), sumB);

        //asumarr2
        acb_set(acb_mat_entry(asumarr2, v, 0), sumA);

        acb_zero(sumA);
        acb_zero(sumB);
        for (n1 = npre+(numn0-1)*H+H+1; n1 <= N ; n1++)
        {
            acb_set_si(nacb, n1);
            acb_log(b, nacb, prec);
            acb_mul_2exp_si(b, b, -2);
            acb_mul(b, b, t, prec);
            acb_add(c, b, acb_mat_entry(ybexpolarge, v, 0), prec);
            acb_pow(a, nacb, c, prec);
            acb_add(sumB, sumB, a, prec);
            acb_add(c, b, acb_mat_entry(yaexpolarge, v, 0), prec);
            acb_pow(a, nacb, c, prec);
            acb_add(sumA, sumA, a, prec);
        }
        //bsumarr3
        acb_set(acb_mat_entry(bsumarr3, v, 0), sumB);

        //asumarr3
        acb_set(acb_mat_entry(asumarr3, v, 0), sumA);
    }

    acb_mat_init(bsums, num, 1);
    acb_mat_init(asums, num, 1);
    //bsums
    acb_mat_add(bsums, bsumarr1, bsumarr2, prec);
    acb_mat_add(bsums, bsums, bsumarr3, prec);
    //asums
    acb_mat_add(asums, asumarr1, asumarr2, prec);
    acb_mat_add(asums, asums, asumarr3, prec);

    acb_mat_init(ests, num, 1);
    //ests
    acb_mat_mul_entrywise(ests, afac, asums, prec);
    acb_mat_add(ests, ests, bsums, prec);

    for (v = 0; v < num; v++)
    { 

       acb_printn(acb_mat_entry(ests, v, 0), digits, ARB_STR_NO_RADIUS);
       flint_printf("\n");
    }
    res = 0;

    acb_mat_clear(n0);
    acb_mat_clear(thtarr);
    acb_mat_clear(n0htmat);
    acb_mat_clear(n0hlogmat);
    acb_mat_clear(sarr);
    acb_mat_clear(afac);
    acb_mat_clear(n0hmat);
    acb_mat_clear(n0sumsbmat);
    acb_mat_clear(n0sumsamat);
    acb_mat_clear(ybexpo);
    acb_mat_clear(yaexpo);
    acb_mat_clear(ybexpolarge);
    acb_mat_clear(yaexpolarge);
    acb_mat_clear(n0thtbmat);
    acb_mat_clear(n0thtamat);
    acb_mat_clear(bsumarr1);
    acb_mat_clear(bsumarr2);
    acb_mat_clear(bsumarr3);
    acb_mat_clear(bsums);
    acb_mat_clear(asumarr1);
    acb_mat_clear(asumarr2);
    acb_mat_clear(asumarr3);
    acb_mat_clear(asums);
    acb_mat_clear(ests);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(d);
    acb_clear(sumA);
    acb_clear(sumB);
    acb_clear(sumA1);
    acb_clear(sumB1);
    acb_clear(one);
    acb_clear(expoacb);
    acb_clear(npreacb);
    acb_clear(numacb);
    acb_clear(hacb);
    acb_clear(Hacb);
    acb_clear(idxacb);
    acb_clear(nacb);
    acb_clear(vacb);

    acb_poly_clear(polya);
    acb_poly_clear(polyb);
    acb_poly_clear(polyc);
    acb_poly_clear(polyd);
    acb_poly_clear(polys);
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
    arb_t x, y, t, a,  d, re, im, out;
    acb_t cx, cy, ct, out1;
    const char *x_str, *y_str, *t_str, *num_str, *H_str, *expterms_str, *npre_str, *d_str;

    slong N, prec, digits, res, num, H, expterms, npre;
    int result = EXIT_SUCCESS;
    int usage_error = 0;
    res = 1;
    arb_init(x);
    arb_init(y);
    arb_init(t);
    arb_init(a);
    arb_init(d);
    arb_init(re);
    arb_init(im);
    arb_init(out);

    acb_init(cx);
    acb_init(cy);
    acb_init(ct);

    acb_init(out1);
 
    if (argc != 9)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }

    x_str = argv[1];
    y_str = argv[2];
    t_str = argv[3];
    num_str = argv[4];
    H_str = argv[5];
    expterms_str = argv[6];
    npre_str = argv[7];
    d_str = argv[8];

    num = atol(num_str);
    H = atol(H_str);
    expterms = atol(expterms_str);
    npre = atol(npre_str);
    digits = atol(d_str);
    prec = 128;

    arb_set_str(x, x_str, prec);
    arb_set_str(y, y_str, prec);
    arb_set_str(t, t_str, prec);
    acb_set_arb(cx, x);
    acb_set_arb(cy, y);
    acb_set_arb(ct, t);

    N = get_N(t, x, prec);
    printf("\n");
    printf("N = : %ld\n",N);
    printf("\n");

    abbeff_multieval_xseries(res, cx, cy, ct, N, num, H, expterms, npre, digits, 1, prec);

//    flint_printf("\n");
//    printf("N = : %ld\n", N);
//    flint_printf("\n");
//    abbeff_series(out1, x, y, t, N, 1, prec);
//    acb_abs(out,out1,prec);
//    flint_printf("(A+B)/B0 eff Re: ");
//    arf_printd(arb_midref(acb_realref(out1)), digits);
//    flint_printf("\n");
//    flint_printf("(A+B)/B0 eff Im: ");
//    arf_printd(arb_midref(acb_imagref(out1)), digits);
//    flint_printf("\n");
//
//    flint_printf("|(A+B)/B0|     : ");
//    arf_printd(arb_midref(out), digits);
//    flint_printf("\n");

finish:
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Usage:\n");
        flint_printf("%s x y t num H expterms npre digits \n\n", argv[0]);
        flint_printf(
    "  This script computes: abbeff_multieval_general_xconst\n"
    "           which is a fast approaximation for:\n"
    "    H_t(z) (eff) = (Aeff_t(z) + Beff_t(z))/Beff0_t(z) \n"
    "             (more informtation will follow) \n"
    "                         \n");
    }
 
    arb_clear(x);
    arb_clear(y);
    arb_clear(t);
    arb_clear(a);
    arb_clear(d);
    arb_clear(re);
    arb_clear(im);
    arb_clear(out);
 
   acb_clear(out1);
    acb_clear(cx);
    acb_clear(cy);
    acb_clear(ct);
 
    flint_cleanup();

    return result;
}
