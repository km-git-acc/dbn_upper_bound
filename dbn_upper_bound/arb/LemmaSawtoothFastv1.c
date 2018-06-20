/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
#include "arb_poly.h"
#include "acb_poly.h"
#include "flint/profiler.h"
 
 
typedef arb_poly_struct *arb_poly_ptr;
typedef const arb_poly_struct *arb_poly_srcptr;
 

//procedure to calculate the absolute value of an arb_poly function
void
arb_poly_abs_series(arb_poly_t res, const arb_poly_t p, slong n, slong prec)
{
    arb_t a, one;
    arb_init(a);
    arb_init(one);

    arb_set_si(one, 1);

    arb_poly_evaluate(a, p, one, prec);

    arb_abs(a, a);

    arb_poly_set_arb(res, a);

    arb_clear(a);
    arb_clear(one);
}

//procedure to calculate the absolute value of an acb_poly function
void
acb_poly_abs_series(arb_poly_t res, const acb_poly_t p, slong n, slong prec)
{
    acb_t a, one; 
    acb_init(a);
    acb_init(one);

    arb_t abs;
    arb_init(abs);
 
    acb_set_si(one, 1);
    acb_poly_evaluate(a, p, one, prec);
    acb_abs(abs, a, prec);

    arb_poly_set_arb(res, abs);
 
    arb_clear(abs);
    acb_clear(a);
    acb_clear(one);
}


//procedure to estabish the max value of the real parts of two arb_poly functions
void
arb_poly_max_series(arb_poly_t res, const arb_poly_t p, const arb_poly_t q, slong n, slong prec)
{
    arb_t a, b, one;
 
    arb_init(a);
    arb_init(b);
    arb_init(one);
    arb_set_si(one, 1);

    arb_poly_evaluate(a, p, one, prec);

    arb_poly_evaluate(b, q, one, prec);

    arb_max(a, a, b, prec);

    arb_poly_set_arb(res, a);

    arb_clear(a);
    arb_clear(b);
    arb_clear(one);
}

//procedure to derive the x-value from the N-value
void
_arb_poly_xNpoly_series(arb_poly_t res, const arb_poly_t N, const arb_poly_t t, slong n, slong prec)
{
    arb_poly_t p, q;
    arb_t pi;
 
    arb_init(pi);
    arb_const_pi(pi, prec);
 
    arb_poly_init(p);
    arb_poly_mullow(p, N, N, n, prec);
    arb_poly_scalar_mul(p, p, pi, prec);
    arb_poly_scalar_mul_2exp_si(p, p, 2);
 
    arb_poly_init(q);
    arb_poly_scalar_mul(q, t, pi, prec);
    arb_poly_scalar_mul_2exp_si(q, q, -2);
 
    arb_poly_sub_series(res, p, q, n, prec);
 
    arb_clear(pi);
    arb_poly_clear(p);
    arb_poly_clear(q);
}

//procedure to calculate the modgamma-value
static void
_arb_poly_modgamend_series(arb_poly_t res, const arb_poly_t Nend, const arb_poly_t t, const arb_poly_t y, slong n, slong prec)
{
    arb_poly_t a, b, c, pip;
	arb_t d, pi;
 
    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(pip);

    arb_init(d);
    arb_init(pi);
    arb_const_pi(pi, prec);

    arb_set_ui(d, 20);
    arb_div_ui(d, d, 1000, prec);

    arb_poly_set_arb(pip, pi);

    arb_poly_scalar_mul(a, y, d, prec);
    arb_poly_exp_series(a, a, n, prec);

    arb_poly_pow_ui(b, Nend, 2, prec);
    arb_poly_scalar_mul_2exp_si(c, pip, 2);
    arb_poly_mullow(b, b, c, n, prec);

    arb_poly_scalar_mul_2exp_si(c, pip, -2);
    arb_poly_mullow(c, t, c, n, prec);

    arb_poly_sub_series(c, b, c, n, prec);
    arb_poly_div_series(c, c, pip, n, prec);
    arb_poly_scalar_mul_2exp_si(c, c, -2);

    arb_poly_neg(b, y);
    arb_poly_scalar_mul_2exp_si(b, b, -1);
    arb_poly_pow_series(c, c, b, n, prec);

    arb_poly_mullow(res, a, c, n, prec);
 
    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(pip);

    arb_clear(d);
    arb_clear(pi);
}


//procedure to calculate the modgamma-value
static void
_arb_poly_modgamma_series(arb_poly_t res, const arb_poly_t xN, const arb_poly_t y, slong n, slong prec)
{
    arb_poly_t a, b, c;
    arb_t d, pi;
 
    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);

    arb_init(d);
    arb_init(pi);
    arb_const_pi(pi, prec);

    arb_set_ui(d, 20);
    arb_div_ui(d, d, 1000, prec);

    arb_poly_scalar_mul(a, y, d, prec);
    arb_poly_exp_series(a, a, n, prec);

    arb_poly_scalar_div(b, xN, pi, prec);
    arb_poly_scalar_mul_2exp_si(b, b, -2);

    arb_poly_neg(c, y);
    arb_poly_scalar_mul_2exp_si(c, c, -1);
    arb_poly_pow_series(c,b,c,n,prec);

    arb_poly_mullow(res, a, c, n, prec);
 
    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_clear(d);
    arb_clear(pi);
}
 
//procedure to calculate the sig exponent-value 
static void
_arb_poly_sig_series(arb_poly_t res, const arb_poly_t xN, const arb_poly_t t, const arb_poly_t y, slong n, slong prec)
{
    arb_poly_t a, b, c, d, e, f, g, z;

    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(d);
    arb_poly_init(e);
    arb_poly_init(f);
    arb_poly_init(g);
    arb_poly_init(z);

    arb_t pi, three;
    arb_init(three);
    arb_set_ui(three, 3);
    arb_poly_zero(z);
    arb_init(pi);
    arb_const_pi(pi, prec);
 
    arb_poly_add_si(a, y, 1, prec);
    arb_poly_scalar_mul_2exp_si(a, a, -1);

    arb_poly_scalar_div(b, xN, pi, prec);
    arb_poly_scalar_mul_2exp_si(b, b, -2);
    arb_poly_log_series(b, b, n, prec);
    arb_poly_mullow(b, b, t, n, prec);
    arb_poly_scalar_mul_2exp_si(b, b, -2);
    arb_poly_add_series(res, a, b, n, prec);

    arb_poly_mullow(c, xN, xN, n, prec);
    arb_poly_scalar_mul_2exp_si(c, c, 1);
    arb_poly_div_series(c, t, c, n, prec);

    arb_poly_scalar_mul_2exp_si(d, y, 2);
    arb_poly_add_si(e, y, 1, prec);
    arb_poly_mullow(e, e, d, n, prec);
    arb_poly_mullow(f, xN, xN, n, prec);
    arb_poly_div_series(f, e, f, n, prec);

    arb_poly_neg(g, y);
	arb_poly_scalar_mul(g, g, three, prec);
    arb_poly_add_si(g, g, 1, prec);
    arb_poly_add_series(g, g, f, n, prec);

    arb_poly_max_series(g, g, z, n, prec);

    arb_poly_mullow(g, c, g, n, prec);
    arb_poly_sub_series(res, res, g, n, prec);

    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(d);
    arb_poly_clear(e);
    arb_poly_clear(f);
    arb_poly_clear(g);
    arb_poly_clear(z);
    arb_clear(three);
    arb_clear(pi);

}

//procedure to calculate the bt-value 
static void
_arb_poly_bt_series(arb_poly_t res, const arb_t logk, const arb_poly_t t, slong n, slong prec)
{
    arb_poly_scalar_mul(res, t, logk, prec);
    arb_poly_scalar_mul(res, res, logk, prec);
    arb_poly_scalar_mul_2exp_si(res, res, -2);
    arb_poly_exp_series(res, res, n, prec);
}

//procedure to calculate the at-value  
static void
_arb_poly_at_series(arb_poly_t res, const arb_t logk, const arb_poly_t t, const arb_poly_t y, slong n, slong prec)
{
    arb_poly_t p;
    arb_poly_init(p);

    arb_poly_scalar_mul(p, t, logk, prec);
    arb_poly_scalar_mul_2exp_si(p, p, -2);
    arb_poly_add_series(p, p, y, n, prec);
 
    arb_poly_scalar_mul(res, p, logk, prec);
    arb_poly_exp_series(res, res, n, prec);
 
    arb_poly_clear(p);
}

//procedure to establish some max values  
static void
lbound_numerator(arb_poly_t summandnd, const arb_poly_t ndsum_b, const arb_poly_t ndsum_a, 
                 arb_poly_t modgam, arb_poly_t modgamend, slong n, slong prec)
{
    arb_poly_t a, b, c, M1, M2, M3;
    arb_poly_init(a); 
    arb_poly_init(b); 
    arb_poly_init(c);
    arb_poly_init(M1);
    arb_poly_init(M2);
    arb_poly_init(M3);

    arb_poly_one(a);
    arb_poly_sub_series(b, a, modgamend, n, prec);
    arb_poly_add_series(c, a, modgamend, n, prec);
    arb_poly_div_series(a, b, c, n, prec);
    arb_poly_add_series(b, ndsum_b, ndsum_a, n, prec);
    arb_poly_abs_series(b, b, n, prec);
    arb_poly_mullow(M1, b, a, n, prec);

    arb_poly_sub_series(a, ndsum_b, ndsum_a, n, prec);
    arb_poly_abs_series(M2, a, n, prec);

    arb_poly_div_series(a, modgamend, modgam, n, prec);
    arb_poly_mullow(b, a, ndsum_a, n, prec);
    arb_poly_sub_series(a, ndsum_b, b, n, prec);
    arb_poly_abs_series(M3, a, n, prec);

    arb_poly_max_series(a, M1, M2, n, prec);
    arb_poly_max_series(summandnd, a, M3, n, prec);

    arb_poly_clear(a); 
    arb_poly_clear(b); 
    arb_poly_clear(c); 
    arb_poly_clear(M1);
    arb_poly_clear(M2);
    arb_poly_clear(M3);
}

//procedure to calculate the eA upper error bound 
static void
_arb_poly_eA_series(arb_poly_t eout, const arb_poly_t xN, const arb_poly_t t, const arb_poly_t y, slong N, slong n, slong prec)
{
    arb_poly_t a, b, c, kpoly, kterm, Npoly, zerod313, threed33, summand;
    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(kpoly);
    arb_poly_init(kterm);
    arb_poly_init(Npoly);
    arb_poly_init(zerod313);
    arb_poly_init(threed33);
    arb_poly_init(summand);

    arb_t logk, pi;
    arb_init(logk);
    arb_init(pi);
    arb_const_pi(pi, prec);

    slong k;

    arb_poly_set_si(Npoly, N);

    arb_poly_set_si(a, 313);
    arb_poly_set_si(b, 1000);
    arb_poly_div_series(zerod313, a, b, n, prec);

    arb_poly_set_si(a, 333);
    arb_poly_set_si(b, 100);
    arb_poly_div_series(threed33, a, b, n, prec);

    arb_poly_zero(eout);
    for (k = 1; k <= N; k++)
    {
        arb_poly_set_si(kpoly, k);
        arb_poly_pow_series(a, kpoly, y, n, prec);

        arb_set_si(logk, k);
        arb_log(logk, logk, prec);

        _arb_poly_bt_series(b, logk, t, n, prec);

        arb_poly_mullow(a, a, b, n, prec);

        _arb_poly_sig_series(c, xN, t, y, n, prec);

        arb_poly_scalar_mul(kterm, c, logk, prec);
        arb_poly_neg(kterm, kterm);
        arb_poly_exp_series(kterm, kterm, n, prec);
        arb_poly_mullow(a, a, kterm, n, prec);

        arb_poly_pow_ui(b, t, 2, prec);
        arb_poly_scalar_mul_2exp_si(b, b, -5);

        arb_poly_pow_ui(c, kpoly, 2, prec);  
        arb_poly_scalar_mul_2exp_si(c, c, 2);
        arb_poly_scalar_mul(c, c, pi, prec);   
        arb_poly_div_series(c, xN, c, n, prec);    
        arb_poly_log_series(c, c, n, prec); 
        arb_poly_pow_ui(c, c, 2, prec);  

        arb_poly_mullow(c, c, b, n, prec);
        arb_poly_add_series(c, c, zerod313, n, prec);  

        arb_poly_scalar_mul_2exp_si(b, xN, -1);
        arb_poly_sub_series(b, b, threed33, n, prec);  

        arb_poly_div_series(c, c, b, n, prec); 
        arb_poly_exp_series(c, c, n, prec); 
        arb_poly_add_si(c, c, -1, prec);

        arb_poly_mullow(summand, c, a, n, prec);

        arb_poly_add_series(eout, eout, summand, n, prec);
    }

    //define modK
    arb_poly_add_si(a, xN, -6, prec);
    arb_poly_scalar_mul_2exp_si(a, a, 1);
    arb_poly_div_series(a, y, a, n, prec);
    arb_poly_mullow(a, t, a, n, prec);
    arb_poly_abs_series(a, a, n, prec);
    arb_poly_pow_series(a, Npoly, a, n, prec);

    _arb_poly_modgamma_series(b, xN, y, n, prec);
    arb_poly_abs_series(b, b, n, prec);

    arb_poly_mullow(b, b, a, n, prec);
    arb_poly_mullow(eout, eout, b, n, prec);

    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(kpoly);
    arb_poly_clear(kterm);
    arb_poly_clear(Npoly);
    arb_poly_clear(zerod313);
    arb_poly_clear(threed33);
    arb_poly_clear(summand);

    arb_clear(pi);
    arb_clear(logk);
}

//procedure to calculate the eB upper error bound 
static void
_arb_poly_eB_series(arb_poly_t eout, const arb_poly_t xN, const arb_poly_t t, const arb_poly_t y, slong N, slong n, slong prec)
{
    arb_poly_t a, b, c, kpoly, kterm, Npoly, zerod313, threed33, summand;
    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(kpoly);
    arb_poly_init(kterm);
    arb_poly_init(Npoly);
    arb_poly_init(zerod313);
    arb_poly_init(threed33);
    arb_poly_init(summand);

    arb_t logk, pi;
    arb_init(logk);
    arb_init(pi);
    arb_const_pi(pi, prec);

    slong k;

    arb_poly_set_si(Npoly, N);

    arb_poly_set_si(a, 313);
    arb_poly_set_si(b, 1000);
    arb_poly_div_series(zerod313, a, b, n, prec);

    arb_poly_set_si(a, 333);
    arb_poly_set_si(b, 100);
    arb_poly_div_series(threed33, a, b, n, prec);

    arb_poly_zero(eout);
    for (k = 1; k <= N; k++)
    {
        arb_poly_set_si(kpoly, k);
        arb_set_si(logk, k);
        arb_log(logk, logk, prec);

        _arb_poly_bt_series(a, logk, t, n, prec);

        _arb_poly_sig_series(c, xN, t, y, n, prec);
        arb_poly_scalar_mul(kterm, c, logk, prec);
        arb_poly_neg(kterm, kterm);
        arb_poly_exp_series(kterm, kterm, n, prec);

        arb_poly_mullow(a, a, kterm, n, prec);

        arb_poly_pow_ui(b, t, 2, prec);
        arb_poly_scalar_mul_2exp_si(b, b, -5);

        arb_poly_pow_ui(c, kpoly, 2, prec);  
        arb_poly_scalar_mul_2exp_si(c, c, 2);
        arb_poly_scalar_mul(c, c, pi, prec);   
        arb_poly_div_series(c, xN, c, n, prec);    
        arb_poly_log_series(c, c, n, prec); 
        arb_poly_pow_ui(c, c, 2, prec);  

        arb_poly_mullow(c, c, b, n, prec);
        arb_poly_add_series(c, c, zerod313, n, prec);  

        arb_poly_scalar_mul_2exp_si(b, xN, -1);
        arb_poly_sub_series(b, b, threed33, n, prec);  

        arb_poly_div_series(c, c, b, n, prec); 
        arb_poly_exp_series(c, c, n, prec); 
        arb_poly_add_si(c, c, -1, prec);

        arb_poly_mullow(summand, c, a, n, prec);

        arb_poly_add_series(eout, eout, summand, n, prec);
    }

    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(kpoly);
    arb_poly_clear(kterm);
    arb_poly_clear(Npoly);
    arb_poly_clear(zerod313);
    arb_poly_clear(threed33);
    arb_poly_clear(summand);
    arb_clear(pi);
    arb_clear(logk);
}

//procedure to calculate the eBs (step) upper error bound 
static void
_arb_poly_eBs_series(arb_poly_t eout, const arb_poly_t xN0, const arb_poly_t t, const arb_poly_t y, slong N0, slong Nstep, slong n, slong prec)
{
    arb_poly_t a, b, c, pip, kpoly, kterm, Npoly, zerod313, threed33, summand;
    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(pip);
    arb_poly_init(kpoly);
    arb_poly_init(kterm);
    arb_poly_init(Npoly);
    arb_poly_init(zerod313);
    arb_poly_init(threed33);
    arb_poly_init(summand);

    arb_t logk, pi;
    arb_init(logk);
    arb_init(pi);
    arb_const_pi(pi, prec);
    arb_poly_set_arb(pip, pi);

    slong k;

    arb_poly_set_si(Npoly, N0);

    arb_poly_set_si(a, 313);
    arb_poly_set_si(b, 1000);
    arb_poly_div_series(zerod313, a, b, n, prec);

    arb_poly_set_si(a, 333);
    arb_poly_set_si(b, 100);
    arb_poly_div_series(threed33, a, b, n, prec);

    arb_poly_zero(eout);
    for (k = 1; k <= (N0 + Nstep); k++)
    {
        arb_poly_set_si(kpoly, k);
        arb_set_si(logk, k);
        arb_log(logk, logk, prec);

        _arb_poly_bt_series(a, logk, t, n, prec);

        _arb_poly_sig_series(c, xN0, t, y, n, prec);
        arb_poly_scalar_mul(kterm, c, logk, prec);
        arb_poly_neg(kterm, kterm);
        arb_poly_exp_series(kterm, kterm, n, prec);

        arb_poly_mullow(a, a, kterm, n, prec);

        arb_poly_pow_ui(b, t, 2, prec);
        arb_poly_scalar_mul_2exp_si(b, b, -5);

        arb_poly_scalar_mul_2exp_si(c, pip, 2);  
        arb_poly_div_series(c, xN0, c, n, prec);    
        arb_poly_log_series(c, c, n, prec); 
        arb_poly_pow_ui(c, c, 2, prec);  

        arb_poly_mullow(c, c, b, n, prec);
        arb_poly_add_series(c, c, zerod313, n, prec);  

        arb_poly_scalar_mul_2exp_si(b, xN0, -1);
        arb_poly_sub_series(b, b, threed33, n, prec);  

        arb_poly_div_series(c, c, b, n, prec); 
        arb_poly_exp_series(c, c, n, prec); 
        arb_poly_add_si(c, c, -1, prec);

        arb_poly_mullow(summand, c, a, n, prec);

        arb_poly_add_series(eout, eout, summand, n, prec);
    }

    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(pip);
    arb_poly_clear(kpoly);
    arb_poly_clear(kterm);
    arb_poly_clear(Npoly);
    arb_poly_clear(zerod313);
    arb_poly_clear(threed33);
    arb_poly_clear(summand);
    arb_clear(pi);
    arb_clear(logk);
}


//procedure to calculate the eC0 upper error bound 
static void
_arb_poly_eC0_series(arb_poly_t eout, const arb_poly_t xN, const arb_poly_t t, const arb_poly_t y, slong N, slong n, slong prec)
{
    arb_poly_t a, b, c, d, pip, Npoly, logxdiv4pi, ecterm1, ecterm2, three, zerod125, oned24, threed58, sixd66, sixd92, eightd52;
    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(d);
    arb_poly_init(pip);
    arb_poly_init(Npoly);
    arb_poly_init(logxdiv4pi);
    arb_poly_init(ecterm1);
    arb_poly_init(ecterm2);
    arb_poly_init(three);
    arb_poly_init(zerod125);
    arb_poly_init(oned24);
    arb_poly_init(threed58);
    arb_poly_init(sixd66);
    arb_poly_init(sixd92);
    arb_poly_init(eightd52);

    arb_t pi;
    arb_init(pi);
    arb_const_pi(pi, prec);
    arb_poly_set_arb(pip, pi);

    acb_poly_t tmp, tmp1;
    acb_poly_init(tmp);
    acb_poly_init(tmp1);

    arb_poly_set_si(Npoly, N);

    arb_poly_one(a);   
    arb_poly_scalar_mul_2exp_si(zerod125, a, -3);

    arb_poly_set_si(three, 3);

    arb_poly_set_si(a, 124);
    arb_poly_set_si(b, 100);
    arb_poly_div_series(oned24, a, b, n, prec);

    arb_poly_set_si(a, 358);
    arb_poly_set_si(b, 100);
    arb_poly_div_series(threed58, a, b, n, prec);

    arb_poly_set_si(a, 666);
    arb_poly_set_si(b, 100);
    arb_poly_div_series(sixd66, a, b, n, prec);

    arb_poly_set_si(a, 692);
    arb_poly_set_si(b, 100);
    arb_poly_div_series(sixd92, a, b, n, prec);

    arb_poly_set_si(a, 852);
    arb_poly_set_si(b, 100);
    arb_poly_div_series(eightd52, a, b, n, prec);
  
    arb_poly_div_series(a, xN, pip, n, prec);  
    arb_poly_scalar_mul_2exp_si(a, a, -2);  
    arb_poly_log_series(logxdiv4pi, a, n, prec); 
    arb_poly_pow_ui(a, logxdiv4pi, 2, prec);  

    arb_poly_neg(b, t);
    arb_poly_scalar_mul_2exp_si(b, b, -4);  
    arb_poly_mullow(b, b, a, n, prec);

    arb_poly_zero(a); 
    arb_poly_scalar_mul_2exp_si(d, pip, -1); 
    acb_poly_set2_arb_poly(tmp, a, d); 
    acb_poly_set_arb_poly(tmp1, logxdiv4pi);
    acb_poly_add_series(tmp, tmp, tmp1, n, prec);
    acb_poly_abs_series(c, tmp, n, prec);
    arb_poly_mullow(c, three, c, n, prec);
    arb_poly_add_series(c, c, threed58, n, prec);

    arb_poly_sub_series(a, xN, eightd52, n, prec);
    arb_poly_div_series(a, c, a, n, prec);
    arb_poly_add_series(a, a, b, n, prec);
    arb_poly_exp_series(ecterm1, a, n, prec);

    arb_poly_pow_series(a, three, y, n, prec);
    arb_poly_neg(b, y);
    arb_poly_pow_series(b, three, b, n, prec);
    arb_poly_add_series(a, a, b, n, prec);
    arb_poly_mullow(a, a, oned24, n, prec);

    arb_poly_sub_series(b, Npoly, zerod125, n, prec);
    arb_poly_div_series(a, a, b, n, prec);

    arb_poly_sub_series(b, xN, sixd66, n, prec);
    arb_poly_div_series(b, sixd92, b, n, prec);

    arb_poly_add_series(ecterm2, b, a, n, prec); 
    arb_poly_add_si(ecterm2, ecterm2, 1, prec); 

    arb_poly_add_si(a, y, 1, prec);   
    arb_poly_scalar_mul_2exp_si(a, a, -2);  
    arb_poly_neg(a, a); 

    arb_poly_div_series(b, xN, pip, n, prec);  
    arb_poly_scalar_mul_2exp_si(b, b, -2);

    arb_poly_pow_series(a, b, a, n, prec);

    arb_poly_mullow(eout, a, ecterm1, n, prec);
    arb_poly_mullow(eout, eout, ecterm2, n, prec);

    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(d);
    arb_poly_clear(pip);
    arb_poly_clear(Npoly);
    arb_poly_clear(logxdiv4pi);
    arb_poly_clear(ecterm1);
    arb_poly_clear(ecterm2);
    arb_poly_clear(three);
    arb_poly_clear(zerod125);
    arb_poly_clear(oned24);
    arb_poly_clear(threed58);
    arb_poly_clear(sixd66);
    arb_poly_clear(sixd92);
    arb_poly_clear(eightd52);

    acb_poly_clear(tmp);
    acb_poly_clear(tmp1);

    arb_clear(pi);
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

//main procedure to establish incremental 'sawtooth' Lemma-bound  
void
arb_poly_abbeff_sawtooth_lbound(arb_poly_t lbound, const arb_poly_t t, const arb_poly_t y, 
                                arb_poly_srcptr moll,
                                slong nprimes, const ulong *primes, const slong *d,
                                slong D, slong nextN, arb_poly_t modgamend, arb_poly_t prevbound, slong n, slong prec)
{
    slong divs;

    arb_poly_t ndsum_a, ndsummand_a, ndsummand, a, b, btpoly, atpoly;
    arb_poly_t ndsum_b, ndsummand_b, modK, sumnd, summandnd, sum2;
    arb_poly_t nextNpoly;
    arb_poly_t sig, modgam, xN, kterm, kterm1, ndterm, modmol;

    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(btpoly);
    arb_poly_init(atpoly);
    arb_poly_init(modK);
    arb_poly_init(sumnd);
    arb_poly_init(summandnd);
    arb_poly_init(sum2);
    arb_poly_init(ndsum_a);
    arb_poly_init(ndsummand_a);
 
    arb_poly_init(ndsum_b);
    arb_poly_init(ndsummand_b);
    arb_poly_init(ndsummand);
 
    arb_poly_init(nextNpoly);

    arb_poly_init(sig);
    arb_poly_init(modgam);
    arb_poly_init(xN);
    arb_poly_init(kterm);
    arb_poly_init(kterm1);
    arb_poly_init(ndterm);
    arb_poly_init(modmol);

    slong j, k, v, nd;
 
    divs = 1 << nprimes;
	
    arb_t logk, logj, lognd; 
    arb_init(logk);
    arb_init(logj);
    arb_init(lognd);

        //define common
        arb_poly_set_si(nextNpoly, nextN);

        _arb_poly_xNpoly_series(xN, nextNpoly, t, n, prec);
        _arb_poly_sig_series(sig, xN, t, y, n, prec);
        _arb_poly_modgamma_series(modgam, xN, y, n, prec);

        //offset = D*(nextN-1);

        arb_poly_zero(sumnd);
        for (v = 1; v <= divs; v++)
        {

           for (k = 1; k <= d[v]; k++)
           {

              //sub loop from nd=1 to ndivs
              arb_poly_zero(ndsum_a);
              arb_poly_zero(ndsum_b);
              for (nd = 1; nd <= divs; nd++)
              {
                  if (((k+d[v]*(nextN-1)) % d[nd] == 0) && ((k+d[v]*(nextN-1)) <= d[nd]*nextN))
                  {
                      j = (k+d[v]*(nextN-1)) / d[nd];
 
                      arb_set_si(logj, j);
                      arb_log(logj, logj, prec);
 
                      _arb_poly_bt_series(btpoly, logj, t, n, prec);
                      arb_poly_mullow(ndsummand_b, moll + nd, btpoly, n, prec);
                      arb_poly_add_series(ndsum_b, ndsum_b, ndsummand_b, n, prec);

                       _arb_poly_at_series(atpoly, logj, t, y, n, prec);
                      arb_poly_mullow(ndsummand_a, moll + nd, atpoly, n, prec);
                      arb_poly_mullow(ndsummand_a, ndsummand_a, modgam, n, prec);
                      arb_poly_add_series(ndsum_a, ndsum_a, ndsummand_a, n, prec);
                  }
              }

              if (!arb_poly_is_zero(ndsum_b) || !arb_poly_is_zero(ndsum_a))
              {
                  arb_set_si(logk, k+d[v]*(nextN-1));
                  arb_log(logk, logk, prec);
                  arb_poly_scalar_mul(kterm, sig, logk, prec);
                  arb_poly_neg(kterm, kterm);
                  arb_poly_exp_series(kterm, kterm, n, prec);

                  lbound_numerator(summandnd, ndsum_b, ndsum_a, modgam, modgamend, n, prec);

                  arb_poly_mullow(summandnd, kterm, summandnd, n, prec);
                  arb_poly_add_series(sumnd, sumnd, summandnd, n, prec);
               }
            }
        }

       //define modmol
       arb_poly_zero(modmol);
       for (nd = 1; nd <= divs; nd++)
       {
           arb_set_si(lognd, d[nd]);
           arb_log(lognd, lognd, prec);
           arb_poly_scalar_mul(ndterm, sig, lognd, prec);
           arb_poly_neg(ndterm, ndterm);
           arb_poly_exp_series(ndterm, ndterm, n, prec);
 
           arb_poly_mullow(ndsummand, moll + nd, ndterm, n, prec);
           arb_poly_abs_series(ndsummand, ndsummand, n, prec);
           arb_poly_add_series(modmol, modmol, ndsummand, n, prec);
        }

       //define modK
       arb_poly_add_si(a, xN, -6, prec);
       arb_poly_scalar_mul_2exp_si(a, a, 1);
       arb_poly_div_series(a, y, a, n, prec);
       arb_poly_mullow(modK, t, a, n, prec);

       //define sum2
       arb_set_si(logk, nextN);
       arb_log(logk, logk, prec);
       arb_poly_scalar_mul(kterm, modK, logk, prec);
       arb_poly_exp_series(kterm, kterm, n, prec);
       arb_poly_add_si(kterm, kterm, -1, prec);
       arb_poly_sub_series(b, y, sig, n, prec);
       arb_poly_scalar_mul(kterm1, b, logk, prec);
       arb_poly_exp_series(kterm1, kterm1, n, prec);
       arb_poly_mullow(b, kterm1, kterm, n, prec);
       _arb_poly_bt_series(btpoly, logk, t, n, prec);
       arb_poly_mullow(sum2, btpoly, b, n, prec);

       //perform final calculations to complete lbound
       arb_poly_div_series(a, sumnd, modmol, n, prec);
       arb_poly_mullow(b, modgam, sum2, n, prec);
       arb_poly_add_series(a, a, b, n, prec);
       arb_poly_sub_series(lbound, prevbound, a, n, prec);

    arb_clear(logk);
    arb_clear(logj);
    arb_clear(lognd);
 
    arb_poly_clear(ndsum_a);
    arb_poly_clear(ndsummand_a);
 
    arb_poly_clear(ndsum_b);
    arb_poly_clear(ndsummand_b);
    arb_poly_clear(ndsummand);
 
    arb_poly_clear(nextNpoly);

    arb_poly_clear(xN);
    arb_poly_clear(kterm);
    arb_poly_clear(kterm1);
    arb_poly_clear(sumnd);
    arb_poly_clear(summandnd);
    arb_poly_clear(sum2);
    arb_poly_clear(modmol);
    arb_poly_clear(modgam);
    arb_poly_clear(ndterm);
    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(btpoly);
    arb_poly_clear(atpoly);
    arb_poly_clear(modK);
    arb_poly_clear(sig);
}

//main procedure to establish the 'base' Lemma-bound  
void
arb_poly_abbeff_base_lbound(arb_poly_t lbound,
        const arb_poly_t t, const arb_poly_t y,
        arb_poly_srcptr moll,
        slong nprimes, const ulong *primes, const slong *d,
        slong D, slong Nmin, arb_poly_t modgamend, slong n, slong prec)
{
    slong divs;
    arb_poly_t ndsum_a, ndsummand_a, ndsummand, a, b, atpoly, btpoly;
    arb_poly_t ndsum_b, ndsummand_b, modK, sumnd, summandnd, sum2;
    arb_poly_t Nminpoly, Nmaxpoly;
    arb_poly_t sig, modgam, xN, kterm, kterm1, ndterm, modmol;

    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(btpoly);
    arb_poly_init(atpoly);
    arb_poly_init(modK);
    arb_poly_init(sumnd);
    arb_poly_init(summandnd);
    arb_poly_init(sum2);
    arb_poly_init(ndsum_a);
    arb_poly_init(ndsummand_a);
 
    arb_poly_init(ndsum_b);
    arb_poly_init(ndsummand_b);
    arb_poly_init(ndsummand);
 
    arb_poly_init(Nminpoly);
    arb_poly_init(Nmaxpoly);

    arb_poly_init(sig);
    arb_poly_init(modgam);
    arb_poly_init(xN);
    arb_poly_init(kterm);
    arb_poly_init(kterm1);
    arb_poly_init(ndterm);
    arb_poly_init(modmol);

    slong j, k, nd;
 
    divs = 1 << nprimes;
 
    arb_t logk, logj, lognd;
    arb_init(logk);
    arb_init(logj);
    arb_init(lognd);

    //define common
    arb_poly_set_si(Nminpoly, Nmin);

    _arb_poly_xNpoly_series(xN, Nminpoly, t, n, prec);
    _arb_poly_sig_series(sig, xN, t, y, n, prec);
    _arb_poly_modgamma_series(modgam, xN, y, n, prec);

    //main loop from n=1 to D*N
    arb_poly_zero(sumnd);
    for (k = 1; k <= D*Nmin-1; k++)
    {
        arb_poly_zero(ndsum_a);
        arb_poly_zero(ndsum_b);
        //sub loop from nd=1 to ndivs
        for (nd = 1; nd <= divs; nd++)
        {
            if (((k+1) % d[nd] == 0) && ((k+1) <= d[nd]*Nmin))
            {
                j = (k+1) / d[nd];
 
                arb_set_si(logj, j);
                arb_log(logj, logj, prec);
 
                _arb_poly_bt_series(btpoly, logj, t, n, prec);
                arb_poly_mullow(ndsummand_b, moll + nd, btpoly, n, prec);
                arb_poly_add_series(ndsum_b, ndsum_b, ndsummand_b, n, prec);

                _arb_poly_at_series(atpoly, logj, t, y, n, prec);
                arb_poly_mullow(ndsummand_a, moll + nd, atpoly, n, prec);
                arb_poly_mullow(ndsummand_a, ndsummand_a, modgam, n, prec);
                arb_poly_add_series(ndsum_a, ndsum_a, ndsummand_a, n, prec);
            }
        }

        if (!arb_poly_is_zero(ndsum_b) || !arb_poly_is_zero(ndsum_a))
        {
            arb_set_si(logk, k+1);
            arb_log(logk, logk, prec);
            arb_poly_scalar_mul(kterm, sig, logk, prec);
            arb_poly_neg(kterm, kterm);
            arb_poly_exp_series(kterm, kterm, n, prec);

            lbound_numerator(summandnd, ndsum_b, ndsum_a, modgam, modgamend, n, prec);

            arb_poly_mullow(summandnd, kterm, summandnd, n, prec);
            arb_poly_add_series(sumnd, sumnd, summandnd, n, prec);
        }
    }

    //define modmol
    arb_poly_zero(modmol);
    for (nd = 1; nd <= divs; nd++)
    {
        arb_set_si(lognd, d[nd]);
        arb_log(lognd, lognd, prec);
        arb_poly_scalar_mul(ndterm, sig, lognd, prec);
        arb_poly_neg(ndterm, ndterm);
        arb_poly_exp_series(ndterm, ndterm, n, prec);
 
        arb_poly_mullow(ndsummand, moll + nd, ndterm, n, prec);
        arb_poly_abs_series(ndsummand, ndsummand, n, prec);
        arb_poly_add_series(modmol, modmol, ndsummand, n, prec);
    }

    //define modK
    arb_poly_add_si(a, xN, -6, prec);
    arb_poly_scalar_mul_2exp_si(a, a, 1);
    arb_poly_div_series(a, y, a, n, prec);
    arb_poly_mullow(modK, t, a, n, prec);

    //define sum2
    arb_poly_zero(sum2);
    for (k = 1; k <= Nmin; k++)
    {
        arb_set_si(logk, k);
        arb_log(logk, logk, prec);
        arb_poly_scalar_mul(kterm, modK, logk, prec);
        arb_poly_exp_series(kterm, kterm, n, prec);
        arb_poly_add_si(kterm, kterm, -1, prec);
        arb_poly_sub_series(b, y, sig, n, prec);
        arb_poly_scalar_mul(kterm1, b, logk, prec);
        arb_poly_exp_series(kterm1, kterm1, n, prec);
        arb_poly_mullow(b, kterm1, kterm, n, prec);
       _arb_poly_bt_series(btpoly, logk, t, n, prec);
        arb_poly_mullow(b, btpoly, b, n, prec);
        arb_poly_add_series(sum2, sum2, b, n, prec);
    }

    //perform final calculations to complete lbound
    arb_poly_one(a);
    arb_poly_sub_series(a, a, modgam, n, prec);
    arb_poly_sub_series(a, a, sumnd, n, prec);
    arb_poly_div_series(a, a, modmol, n, prec);
    arb_poly_mullow(b, modgam, sum2, n, prec);
    arb_poly_sub_series(lbound, a, b, n, prec);

    arb_clear(logk);
    arb_clear(logj);
    arb_clear(lognd);
 
    arb_poly_clear(ndsum_a);
    arb_poly_clear(ndsummand_a);
 
    arb_poly_clear(ndsum_b);
    arb_poly_clear(ndsummand_b);
    arb_poly_clear(ndsummand);
 
    arb_poly_clear(Nminpoly);
    arb_poly_clear(Nmaxpoly);

    arb_poly_clear(xN);
    arb_poly_clear(kterm);
    arb_poly_clear(kterm1);
    arb_poly_clear(sumnd);
    arb_poly_clear(summandnd);
    arb_poly_clear(sum2);
    arb_poly_clear(modmol);
    arb_poly_clear(modgam);
    arb_poly_clear(ndterm);
    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(btpoly);
    arb_poly_clear(atpoly);
    arb_poly_clear(modK);
    arb_poly_clear(sig);
}

//initialise and fill arrays for bt, at, moll. Process main loop  
void static
main_loop(slong res, const arb_t t, const arb_t y, slong Na, slong Nb, slong Nprint,
          slong nprimes, arb_t threshold, slong mode, slong digits, slong n, slong prec)
{

    slong D, Nmax, nextN;
    //arb_poly_ptr bt, at;
    slong i, j;
    slong high_prec;

    arb_poly_ptr moll;
    slong divs;
    ulong *primes;
    slong *d;

    divs=0; D=0; 
 
    high_prec = digits * 3.32192809488736 + 10;

    Nmax = Nb;
 
    arb_t ar, a, logk, logp, coeff;
    arb_init(ar);
    arb_init(a);
    arb_init(logk);
    arb_init(logp);
    arb_init(coeff);

    arb_poly_t tpoly, ypoly, btpoly, xN, xN0, eA, eB, eC0, eTot, lboundfull, lboundinc, modgamend, Napoly, Nbpoly, nextNpoly, prevbound;
    arb_poly_init(tpoly);
    arb_poly_init(ypoly);
    arb_poly_init(btpoly);
    arb_poly_init(xN0);
    arb_poly_init(xN);
    arb_poly_init(eA);
    arb_poly_init(eB);
    arb_poly_init(eC0);
    arb_poly_init(eTot);
    arb_poly_init(lboundfull);
    arb_poly_init(lboundinc);
    arb_poly_init(modgamend);
    arb_poly_init(Napoly);
    arb_poly_init(Nbpoly);
    arb_poly_init(nextNpoly);
    arb_poly_init(prevbound);
 
    arb_poly_set_arb(tpoly, t);
    arb_poly_set_arb(ypoly, y);
    arb_poly_set_si(Napoly, Na);
    arb_poly_set_si(Nbpoly, Nb);
    _arb_poly_xNpoly_series(xN0, Napoly, tpoly, n, prec);

    //only print the error terms
    if (mode == 2)
    {
       nextN = Na;
       while (nextN <= Nb)
       {
           arb_poly_set_si(Nbpoly, nextN);
           _arb_poly_xNpoly_series(xN, Nbpoly, tpoly, n, prec);

           _arb_poly_eA_series(eA, xN, tpoly, ypoly, nextN, n, prec);
           _arb_poly_eB_series(eB, xN, tpoly, ypoly, nextN, n, prec);
           _arb_poly_eC0_series(eC0, xN, tpoly, ypoly, nextN, n, prec);

           arb_poly_add_series(eTot, eA, eB, n, prec); 
           arb_poly_add_series(eTot, eTot, eC0, n, prec); 

           if (nextN % Nprint == 0)
           {
               flint_printf("%ld, ", nextN);
               arb_poly_get_coeff_arb(coeff, eA, 0);
               arf_printd(arb_midref(coeff), digits);
               flint_printf(", ");
               arb_poly_get_coeff_arb(coeff, eB, 0);
               arf_printd(arb_midref(coeff), digits);
               flint_printf(", ");
               arb_poly_get_coeff_arb(coeff, eC0, 0);
               arf_printd(arb_midref(coeff), digits);
               flint_printf(", ");
               arb_poly_get_coeff_arb(coeff, eTot, 0);
               arf_printd(arb_midref(coeff), digits);
               flint_printf("\n");
           }
           nextN = nextN + Nprint; 
       }
       goto end;
    } 
 
    divs = 1 << nprimes;
    primes = flint_malloc(nprimes * sizeof(*primes));
    _first_n_primes(primes, nprimes);
 
    D = 1;
    for (i = 0; i < nprimes; i++)
    {
        D *= primes[i];
    }
 
    //bt = flint_malloc((Nmax * D + 1) * sizeof(*bt));
    //at = flint_malloc((Nmax * D + 1) * sizeof(*at));
 
//    for (k = 1; k <= Nmax * D; k++)
//    {
//        arb_set_si(logk, k);
//        arb_log(logk, logk, prec);
// 
//        arb_poly_init(bt + k);
//        arb_poly_init(at + k);
// 
//        _arb_poly_bt_series(bt + k, logk, tpoly, n, prec);
//        _arb_poly_at_series(at + k, logk, tpoly, ypoly, n, prec);
//    }
 
    d = flint_malloc((divs + 1) * sizeof(*d));
    moll = flint_malloc((divs + 1) * sizeof(*moll));
    for (i = 0; i < divs; i++)
    {
        arb_poly_init(moll + i + 1);
        arb_poly_one(moll + i + 1);
        d[i + 1] = 1;
        for (j = 0; j < nprimes; j++)
        {
            if (i & (1 << j))
            {
                d[i + 1] *= primes[j];
                arb_set_si(logp, primes[j]);
                arb_log(logp, logp, prec);
 
                _arb_poly_bt_series(btpoly, logp, tpoly, n, prec);
                arb_poly_mullow(
                        moll + i + 1, moll + i + 1, btpoly, n, prec);
                arb_poly_neg(moll + i + 1, moll + i + 1);
            }
        }
    }

    //main loop from N-start to N-end 
    nextN = Na;
    _arb_poly_modgamend_series(modgamend, Nbpoly, tpoly, ypoly, n, prec);

    arb_poly_abbeff_base_lbound(lboundfull, tpoly, ypoly, moll, nprimes, primes, d, D, Na, modgamend, 1, prec);
    arb_poly_set(prevbound, lboundfull);
    flint_printf("%ld, ", nextN);
    arb_poly_get_coeff_arb(coeff, lboundfull, 0);
    arf_printd(arb_midref(coeff), digits);
    if (mode == 1)
    {
       flint_printf(", ");
       _arb_poly_xNpoly_series(xN0, Napoly, tpoly, n, prec);
       _arb_poly_eBs_series(eB, xN0, tpoly, ypoly, Na, Nprint, n, prec);
       _arb_poly_eC0_series(eC0, xN0, tpoly, ypoly, Na, n, prec); 
       arb_poly_add_series(eTot, eC0, eB, n, prec); 
       arb_poly_add_series(eTot, eTot, eB, n, prec);  
       arb_poly_add_series(eTot, eTot, eB, n, prec);   
       arb_poly_get_coeff_arb(coeff, eTot, 0);
       arf_printd(arb_midref(coeff), digits);
    }
    flint_printf(", (init)\n");


    //main loop from N-start to N-end 
    while (nextN < Nb)
    {
        nextN = nextN + 1;   

        arb_poly_abbeff_sawtooth_lbound(lboundinc, tpoly, ypoly, moll, nprimes, primes, d, D, nextN, modgamend, prevbound, 1, prec);

        arb_one(a);
        arb_poly_evaluate(a, lboundinc, a, prec);

        if (arb_lt(a, threshold))
        {
            arb_poly_abbeff_base_lbound(lboundfull, tpoly, ypoly, moll, nprimes, primes, d, D, nextN, modgamend, 1, prec);
            arb_poly_set(prevbound, lboundfull);
            flint_printf("%ld, ", nextN);
            arb_poly_get_coeff_arb(coeff, lboundfull, 0);
            arf_printd(arb_midref((coeff)), digits);
            if (mode == 1)
            {
               flint_printf(", ");
               arb_poly_set_si(nextNpoly, nextN);
               _arb_poly_xNpoly_series(xN0, nextNpoly, tpoly, n, prec);
               _arb_poly_eB_series(eB, xN0, tpoly, ypoly, nextN, n, prec);
               _arb_poly_eC0_series(eC0, xN0, tpoly, ypoly, nextN, n, prec); 
               arb_poly_add_series(eTot, eC0, eB, n, prec); 
               arb_poly_add_series(eTot, eTot, eB, n, prec);  
               arb_poly_add_series(eTot, eTot, eB, n, prec);   
               arb_poly_get_coeff_arb(coeff, eTot, 0);
               arf_printd(arb_midref(coeff), digits);
            }
            flint_printf(", (recalc)\n");
        }
        else
        {
            arb_poly_set(prevbound, lboundinc);
            if (nextN % Nprint == 0)
            {
                flint_printf("%ld, ", nextN);
                arb_poly_get_coeff_arb(coeff, lboundinc, 0);
                arf_printd(arb_midref(coeff), digits);
                if (mode == 1)
                {
                  flint_printf(", ");
                  arb_poly_set_si(nextNpoly, nextN);
                  _arb_poly_xNpoly_series(xN0, nextNpoly, tpoly, n, prec);
                  _arb_poly_eBs_series(eB, xN0, tpoly, ypoly, nextN, Nprint, n, prec);
                  arb_poly_add_series(eTot, eC0, eB, n, prec); 
                  arb_poly_add_series(eTot, eTot, eB, n, prec);  
                  arb_poly_add_series(eTot, eTot, eB, n, prec);   
                  arb_poly_get_coeff_arb(coeff, eTot, 0);
                  arf_printd(arb_midref(coeff), digits);
                }
                flint_printf("\n");
            }
        }
    }

    //clear all arrays
//    for (k = 1; k <= Nmax * D; k++)
//    {
//        arb_poly_clear(bt + k);
//        arb_poly_clear(at + k);
//    }
//    flint_free(bt);
//    flint_free(at);
 
    for (i = 0; i < divs; i++)
    {
        arb_poly_clear(moll + i + 1);
    }
    flint_free(d);
    flint_free(moll);
    flint_free(primes);

end:
    //clear all variables
    arb_clear(ar);
    arb_clear(a);
    arb_clear(logk);
    arb_clear(logp);
    arb_clear(coeff);

    arb_poly_clear(tpoly);
    arb_poly_clear(ypoly);
    arb_poly_clear(btpoly);
    arb_poly_clear(xN);
    arb_poly_clear(xN0);
    arb_poly_clear(eA);
    arb_poly_clear(eB);
    arb_poly_clear(eC0);
    arb_poly_clear(eTot);
    arb_poly_clear(lboundfull);
    arb_poly_clear(lboundinc);
    arb_poly_clear(modgamend);
    arb_poly_clear(Napoly);
    arb_poly_clear(Nbpoly);
    arb_poly_clear(nextNpoly);
    arb_poly_clear(prevbound);
}
 
//main initialisation - processing and validating input values 
int main(int argc, char *argv[])
{
    const char *t_str, *y_str;
    const char *Na_str, *Nb_str, *Nprint_str;
    const char *m_str, *thres_str, *mode_str, *d_str;

    slong Na, Nb, Nprint;
    slong m, mode, d;
    slong prec, res;
    int result = EXIT_SUCCESS;
    int usage_error = 0;
    res = 0;
 
    arb_t t, y, threshold;
    arb_init(t);
    arb_init(y);
    arb_init(threshold);
 
    if (argc != 10)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }
 
    t_str = argv[1];
    y_str = argv[2];
    Na_str = argv[3];
    Nb_str = argv[4];
    Nprint_str = argv[5];
    m_str = argv[6];
    thres_str = argv[7];
    mode_str = argv[8];
    d_str = argv[9];
 
    Na = atol(Na_str);
    Nb = atol(Nb_str);
    Nprint = atol(Nprint_str);
    m = atol(m_str);
    mode = atol(mode_str);
    d = atol(d_str);
 
    if (Na < 1 || Nb < 1 || Nprint < 1 || Na > Nb ||
        m < 0 || m > 20 || mode < 0 ||d < 1)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }
 
    prec = d * 3.32192809488736 + 30;

    arb_set_str(t, t_str, prec);
    arb_set_str(y, y_str, prec);
    arb_set_str(threshold, thres_str, prec);

TIMEIT_ONCE_START

    main_loop(res, t, y, Na, Nb, Nprint, m, threshold, mode, d, 1, prec);

TIMEIT_ONCE_STOP 
 
finish:
 
    if (usage_error && result == EXIT_FAILURE)
    {
        flint_printf("Usage:\n");
        flint_printf("%s t y Na Nb Nprint m threshold mode d\n\n", argv[0]);
        flint_printf(
    "This script computes a lower Lemma bound of abs(Aeff + Beff)/abs(Beff0)\n"
    "for N between 'Na' and 'Nb'.\n"
    "It uses an Euler mollification that includes m primes,\n"
    "so that for example m=0 corresponds to no mollification\n"
    "and m=3 corresponds to an Euler mollification that uses\n"
    "the first three primes {2, 3, 5}.\n"
    "When the incremental bounds reach a certain 'threshold' then\n"
    "the fully accurate bound is recalculated (labelled 'recalc'\n"
    "Each row of output consists of N followed by the bound\n"
    "and output is only printed when 'N mod Nprint == 0'.\n"
    "The script can be run in 3 'modes': .\n"
    "Mode 0 = no upper error bounds are calculated.'\n"
    "Mode 1 = conservative error bounds are calculated and shown at 'Nprint'.\n"
    "Mode 2 = accurate eA, eB, eCO and total error bounds are printed 'Nprint'.\n"
    "Number of decimal digits can be set by 6.\n");
    }
 
    arb_clear(t);
    arb_clear(y);
    arb_clear(threshold);
 
    flint_cleanup();
    return result;
}
