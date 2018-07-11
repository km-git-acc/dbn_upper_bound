/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath
 
    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
 
#include "arb_poly.h"
#include "acb_poly.h"
#include "arb_hypgeom.h"
#include "flint/profiler.h"
 
typedef arb_poly_struct *arb_poly_ptr;
typedef const arb_poly_struct *arb_poly_srcptr;


//evaluate the closed form of the tailintegral
void
arb_poly_tailint(arb_poly_t result, arb_poly_t aterm, arb_poly_t bterm, arb_poly_t ilim, arb_poly_t ulim, 
                 slong n, slong prec)
{  
    arb_poly_t a, b, c, d, e, pipoly;
 
    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(d);
    arb_poly_init(e);
    arb_poly_init(pipoly);

    arb_t pi;
    arb_init(pi);
    arb_const_pi(pi, prec);
    arb_poly_set_arb(pipoly, pi);

    arb_poly_sqrt_series(a, pipoly, n, prec);
    arb_poly_sqrt_series(b, bterm, n, prec);
    arb_poly_div_series(a, a, b, n, prec);
    arb_poly_scalar_mul_2exp_si(a, a, -1);


    arb_poly_add_si(b, aterm, -1, prec);
    arb_poly_mullow(b, b, b, n, prec);
    arb_poly_div_series(b, b, bterm, n, prec);
    arb_poly_scalar_mul_2exp_si(b, b, -2);
    arb_poly_neg(b, b);
    arb_poly_exp_series(b, b, n, prec);
    arb_poly_mullow(a, a, b, n, prec);

    arb_poly_scalar_mul_2exp_si(b, bterm, 1);  
    arb_poly_log_series(c, ulim, n, prec); 
    arb_poly_mullow(b, b, c, n, prec);
    arb_poly_sub_series(b, b, aterm, n, prec);
    arb_poly_add_si(b, b, 1, prec);
    arb_poly_sqrt_series(c, bterm, n, prec);
    arb_poly_div_series(b, b, c, n, prec);
    arb_poly_scalar_mul_2exp_si(b, b, -1);
    arb_hypgeom_erfi_series(d, b, 1, prec);

    arb_poly_scalar_mul_2exp_si(b, bterm, 1);  
    arb_poly_log_series(c, ilim, n, prec); 
    arb_poly_mullow(b, b, c, n, prec);
    arb_poly_sub_series(b, b,aterm, n, prec);
    arb_poly_add_si(b, b, 1, prec);
    arb_poly_sqrt_series(c, bterm, n, prec);
    arb_poly_div_series(b, b, c, n, prec);
    arb_poly_scalar_mul_2exp_si(b, b, -1);
    arb_hypgeom_erfi_series(e, b, 1, prec);

    arb_poly_sub_series(b, d, e, n, prec);
    arb_poly_mullow(result, a, b, n, prec);

    arb_clear(pi);

    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(d);
    arb_poly_clear(e);
    arb_poly_clear(pipoly);
}

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
                 arb_poly_t modgam, slong n, slong prec)
{
    arb_poly_t a, b, c, M1, M2;
    arb_poly_init(a); 
    arb_poly_init(b); 
    arb_poly_init(c);
    arb_poly_init(M1);
    arb_poly_init(M2);

    arb_poly_one(a);
    arb_poly_sub_series(b, a, modgam, n, prec);
    arb_poly_add_series(c, a, modgam, n, prec);
    arb_poly_div_series(a, b, c, n, prec);
    arb_poly_add_series(b, ndsum_b, ndsum_a, n, prec);
    arb_poly_abs_series(b, b, n, prec);
    arb_poly_mullow(M1, b, a, n, prec);

    arb_poly_sub_series(a, ndsum_b, ndsum_a, n, prec);
    arb_poly_abs_series(M2, a, n, prec);

    arb_poly_max_series(summandnd, M1, M2, n, prec);

    arb_poly_clear(a); 
    arb_poly_clear(b); 
    arb_poly_clear(c); 
    arb_poly_clear(M1);
    arb_poly_clear(M2);

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

//main procedure to establish the approximate Triangle-bound moll 2  
void
arb_poly_abbeff_approx_tbound_moll2(arb_poly_t tbound, const arb_poly_t t, const arb_poly_t y,
		                            slong N, slong N0, slong n, slong prec)
{
    slong k;

    arb_poly_t a, b, c, d, two;
    arb_poly_t aterm, bterm, ilim, ulim;
    arb_poly_t s1, s1main, s1summand, s1tail, s1tail1, s1tail2, boundint, boundint1, boundint2, polyN, polyk, modb2;

    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(d);
    arb_poly_init(two);
    arb_poly_init(aterm);
    arb_poly_init(bterm);
    arb_poly_init(ilim);
    arb_poly_init(ulim);
    arb_poly_init(s1);
    arb_poly_init(s1main);
    arb_poly_init(s1summand);
    arb_poly_init(s1tail);
    arb_poly_init(s1tail1);
    arb_poly_init(s1tail2);
    arb_poly_init(boundint);
    arb_poly_init(boundint1);
    arb_poly_init(boundint2);
    arb_poly_init(polyN);
    arb_poly_init(polyk);
    arb_poly_init(modb2);

    //init
    arb_poly_set_si(two, 2);
    arb_poly_set_si(polyN, N);

    //s1main
    arb_poly_zero(s1main);
    for (k = 1; k <= N0; k++)
    {
       arb_poly_set_si(polyk, k);
       arb_poly_div_series(a, polyk, polyN, n, prec);
       arb_poly_pow_series(a, a, y, n, prec);
       arb_poly_add_si(a, a, 1, prec);

       arb_poly_add_si(b, y, 1, prec);
       arb_poly_neg(b, b);
       arb_poly_scalar_mul_2exp_si(b, b, -1);

       arb_poly_mullow(c, polyN, polyN, n, prec);
       arb_poly_div_series(c, c, polyk, n, prec);
       arb_poly_log_series(c, c, n, prec);
       arb_poly_mullow(c, c, t, n, prec);
       arb_poly_scalar_mul_2exp_si(c, c, -2);
       arb_poly_sub_series(b, b, c, n, prec);
       arb_poly_pow_series(b, polyk, b, n, prec);

       arb_poly_mullow(s1summand, a, b, n, prec);
       arb_poly_add_series(s1main, s1main, s1summand, n, prec);
    }

    //s1tail1
    arb_poly_add_si(a, y, 1, prec);
    arb_poly_scalar_mul_2exp_si(a, a, -1);
    arb_poly_log_series(b, polyN, n, prec);
    arb_poly_scalar_mul_2exp_si(c, t, -1);
    arb_poly_mullow(b, b, c, n, prec);
    arb_poly_add_series(aterm, a, b, n, prec);
    arb_poly_scalar_mul_2exp_si(bterm, t, -2);
    arb_poly_set_si(ilim, N0);
    arb_poly_set(ulim, polyN);
    arb_poly_tailint(s1tail1, aterm, bterm, ilim, ulim, n, prec);

    //s1tail2
    arb_poly_neg(a, y);
    arb_poly_add_si(a, a, 1, prec);
    arb_poly_scalar_mul_2exp_si(a, a, -1);
    arb_poly_log_series(b, polyN, n, prec);
    arb_poly_scalar_mul_2exp_si(c, t, -1);
    arb_poly_mullow(b, b, c, n, prec);
    arb_poly_add_series(aterm, a, b, n, prec);
    arb_poly_scalar_mul_2exp_si(bterm, t, -2);
    arb_poly_set_si(ilim, N0);
    arb_poly_set_si(ulim, N);
    arb_poly_tailint(s1tail2, aterm, bterm, ilim, ulim, n, prec);

    //s1tail
    arb_poly_neg(a, y);
    arb_poly_pow_series(a, polyN, a, n, prec);
    arb_poly_mullow(a, a, s1tail2, n, prec);
    arb_poly_add_series(s1tail, s1tail1, a, n, prec);

    //s1
    arb_poly_add_series(s1, s1main, s1tail, n, prec);

    //modb2
    arb_poly_add_si(b, y, 1, prec);
    arb_poly_neg(b, b);
    arb_poly_scalar_mul_2exp_si(b, b, -1);
    arb_poly_mullow(c, polyN, polyN, n, prec);
    arb_poly_div_series(c, c, two, n, prec);
    arb_poly_log_series(c, c, n, prec);
    arb_poly_mullow(c, c, t, n, prec);
    arb_poly_scalar_mul_2exp_si(c, c, -2);
    arb_poly_sub_series(b, b, c, n, prec);
    arb_poly_pow_series(modb2, two, b, n, prec);

    //boundint1
    arb_poly_add_si(a, y, 1, prec);
    arb_poly_scalar_mul_2exp_si(a, a, -1);
    arb_poly_log_series(b, polyN, n, prec);
    arb_poly_scalar_mul_2exp_si(c, t, -1);
    arb_poly_mullow(b, b, c, n, prec);
    arb_poly_add_series(aterm, a, b, n, prec);
    arb_poly_scalar_mul_2exp_si(bterm, t, -2);
    arb_poly_scalar_mul_2exp_si(ilim, polyN, -1);
    arb_poly_set(ulim, polyN);
    arb_poly_tailint(boundint1, aterm, bterm, ilim, ulim, n, prec);

    //bounint2
    arb_poly_neg(a, y);
    arb_poly_add_si(a, a, 1, prec);
    arb_poly_scalar_mul_2exp_si(a, a, -1);
    arb_poly_log_series(b, polyN, n, prec);
    arb_poly_scalar_mul_2exp_si(c, t, -1);
    arb_poly_mullow(b, b, c, n, prec);
    arb_poly_add_series(aterm, a, b, n, prec);
    arb_poly_scalar_mul_2exp_si(bterm, t, -2);
    arb_poly_scalar_mul_2exp_si(ilim, polyN, -1);
    arb_poly_set_si(ulim, N);
    arb_poly_tailint(boundint2, aterm, bterm, ilim, ulim, n, prec);

    //boundint
    arb_poly_neg(a, y);
    arb_poly_pow_series(a, polyN, a, n, prec);
    arb_poly_mullow(a, a, boundint2, n, prec);
    arb_poly_add_series(boundint, boundint1, a, n, prec);

    //complete the tbound
    arb_poly_mullow(a, modb2, boundint, n, prec);
    arb_poly_scalar_mul_2exp_si(a, a, 1);

    arb_poly_neg(b, modb2);
    arb_poly_add_si(b, b, 1, prec);
    arb_poly_mullow(b, b, s1, n, prec);

    arb_poly_add_series(c, a, b, n, prec);
    arb_poly_sub_series(tbound, two, c, n, prec);
 
    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(d);
    arb_poly_clear(two);
    arb_poly_clear(aterm);
    arb_poly_clear(bterm);
    arb_poly_clear(ilim);
    arb_poly_clear(ulim);
    arb_poly_clear(s1);
    arb_poly_clear(s1main);
    arb_poly_clear(s1summand);
    arb_poly_clear(s1tail);
    arb_poly_clear(s1tail1);
    arb_poly_clear(s1tail2);
    arb_poly_clear(boundint);
    arb_poly_clear(boundint1);
    arb_poly_clear(boundint2);
    arb_poly_clear(polyN);
    arb_poly_clear(polyk);
    arb_poly_clear(modb2);
}

//main procedure to establish the approximate Triangle-bound moll0
void
arb_poly_abbeff_approx_tbound_moll0(arb_poly_t tbound, const arb_poly_t t, const arb_poly_t y,
        arb_poly_srcptr moll, slong nprimes, const ulong *primes, const slong *d,
        slong D, slong N, slong N0, slong n, slong prec)
{
    slong divs, j, k, nd;

    arb_poly_t a, b, c, e, atpoly, btpoly, dndpoly, polylogN, polylogdnd, Npoly, N0poly;
    arb_poly_t ndsum_a, ndsummand_a, tailsuma, tailsummanda, sumndb, summandndb, ndsummand;
    arb_poly_t ndsum_b, ndsummand_b, tailsumb, tailsummandb, sumnda, summandnda,sumnd, summandnd, sum2;
    arb_poly_t tdiv2, tdiv4, tdiv4logN, oneminsigb, oneminsiga, tdiv4logNdivnd, mollabsbt, finalterm;
    arb_poly_t sig, modgam, modK, xN, kterm, kterm1, ndterm, modmol;

    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(e);
    arb_poly_init(btpoly);
    arb_poly_init(atpoly);
    arb_poly_init(dndpoly);
    arb_poly_init(Npoly);
    arb_poly_init(N0poly);
    arb_poly_init(polylogN);
    arb_poly_init(polylogdnd);
    arb_poly_init(modK);
    arb_poly_init(sumnd);
    arb_poly_init(ndsummand);
    arb_poly_init(summandnd);
    arb_poly_init(sum2);
	
    arb_poly_init(sumnda);
    arb_poly_init(summandnda);
    arb_poly_init(ndsum_a);
    arb_poly_init(ndsummand_a);
    arb_poly_init(sumndb);
    arb_poly_init(summandndb);
    arb_poly_init(ndsum_b);
    arb_poly_init(ndsummand_b);

    arb_poly_init(tailsuma);
    arb_poly_init(tailsummanda);
    arb_poly_init(tailsumb);
    arb_poly_init(tailsummandb);

    arb_poly_init(tdiv2);
    arb_poly_init(tdiv4);
    arb_poly_init(tdiv4logN);
    arb_poly_init(oneminsigb);
    arb_poly_init(oneminsiga);
    arb_poly_init(tdiv4logNdivnd);
    arb_poly_init(mollabsbt);
    arb_poly_init(finalterm);

    arb_poly_init(sig);
    arb_poly_init(modgam);
    arb_poly_init(xN);
    arb_poly_init(kterm);
    arb_poly_init(kterm1);
    arb_poly_init(ndterm);
    arb_poly_init(modmol);
 
    divs = 1 << nprimes;
 
    arb_t dnd, logk, logj, logN, lognd, logdnd;
    arb_init(dnd);
    arb_init(logk);
    arb_init(logj);
    arb_init(logN);
    arb_init(lognd);
    arb_init(logdnd);

    //define common
    arb_poly_set_si(Npoly, N);
    arb_poly_set_si(N0poly, N0);

    _arb_poly_xNpoly_series(xN, Npoly, t, n, prec);
    _arb_poly_sig_series(sig, xN, t, y, n, prec);
    _arb_poly_modgamma_series(modgam, xN, y, n, prec);

    //main loop from n=1 to N0
    arb_poly_zero(sumnda);
    arb_poly_zero(sumndb);
    for (k = 2; k <= N0; k++)
    {
        arb_poly_zero(ndsum_a);
        arb_poly_zero(ndsum_b);
        //sub loop from nd=1 to ndivs
        for (nd = 1; nd <= divs; nd++)
        {
            if (((k) % d[nd] == 0) && ((k) <= d[nd]*N))
            {
                j = (k) / d[nd];
 
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
            arb_set_si(logk, k);
            arb_log(logk, logk, prec);
            arb_poly_scalar_mul(kterm, sig, logk, prec);
            arb_poly_neg(kterm, kterm);
            arb_poly_exp_series(kterm, kterm, n, prec);

            arb_poly_mullow(summandndb, kterm, ndsum_b, n, prec);
            arb_poly_abs_series(summandndb, summandndb, n, prec);
            arb_poly_add_series(sumndb, sumndb, summandndb, n, prec);

            arb_poly_mullow(summandnda, kterm, ndsum_a, n, prec);
            arb_poly_abs_series(summandnda, summandnda, n, prec);
            arb_poly_add_series(sumnda, sumnda, summandnda, n, prec);
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

    arb_set_si(logN, N);
    arb_log(logN, logN, prec);

    arb_poly_scalar_mul(kterm, modK, logN, prec);
    arb_poly_exp_series(kterm, kterm, n, prec);
    arb_poly_add_si(kterm, kterm, -1, prec);

    _arb_poly_bt_series(btpoly, logN, t, n, prec);
    arb_poly_mullow(b, btpoly, kterm, n, prec);
    arb_poly_mullow(sum2, b, Npoly, n, prec);   

    //calculate tailsums b and a
    arb_poly_zero(tailsuma);
    arb_poly_zero(tailsumb);

       arb_poly_set_arb(polylogN, logN);
       arb_poly_scalar_mul_2exp_si(tdiv2, t, -1);
       arb_poly_scalar_mul_2exp_si(tdiv4, t, -2);

       arb_poly_log_series(a, N0poly, n, prec);
       arb_poly_scalar_mul_2exp_si(b, polylogN, 1);
       arb_poly_sub_series(b, b, a, n, prec);
       arb_poly_mullow(tdiv4logN, tdiv4, b, n, prec);

       arb_poly_add_si(b, y, 1, prec);
       arb_poly_scalar_mul_2exp_si(b, b, -1);
       arb_poly_neg(b, b);
       arb_poly_add_si(oneminsigb, b, 1, prec);

       arb_poly_neg(a, y);
       arb_poly_add_si(b, a, 1, prec);
       arb_poly_scalar_mul_2exp_si(b, b, -1);
       arb_poly_neg(b, b);
       arb_poly_add_si(oneminsiga, b, 1, prec);

       arb_poly_zero(tailsuma);
       arb_poly_zero(tailsumb);
       for (nd = 1; nd <= divs; nd++)
       {
           //common terms
           arb_set_si(dnd, d[nd]);
           arb_log(logdnd, dnd, prec);
           arb_poly_set_arb(polylogdnd, logdnd);
           arb_poly_set_arb(dndpoly, dnd);
           arb_poly_sub_series(a, polylogN, polylogdnd, n, prec);
           arb_poly_mullow(tdiv4logNdivnd, tdiv4, a, n, prec);

           arb_poly_abs_series(a, moll + nd, n, prec);
           _arb_poly_bt_series(btpoly, logdnd, t, n, prec);
           arb_poly_mullow(mollabsbt, a, btpoly, n, prec);

           arb_poly_div_series(b, Npoly, N0poly, n, prec);
           arb_poly_log_series(b, b, 1, prec);

           arb_poly_add_series(b, b, polylogdnd, n, prec);
           arb_poly_scalar_div(b, b, dnd, prec);
           arb_poly_scalar_mul_2exp_si(finalterm, b, -1);

           //tail b-sum
           arb_poly_scalar_mul(c, tdiv2, logdnd, prec);
           arb_poly_pow_series(c, N0poly, c, n, prec);
           arb_poly_div_series(a, mollabsbt, c, n, prec);

           arb_poly_sub_series(b, oneminsigb, tdiv4logN, n, prec);
           arb_poly_pow_series(b, N0poly, b, n, prec);
           arb_poly_scalar_mul(c, Npoly, dnd, prec);
           arb_poly_sub_series(e, oneminsigb, tdiv4logNdivnd, n, prec);
           arb_poly_pow_series(c, c, e, n, prec);
           arb_poly_add_series(b, b, c, n, prec);

           arb_poly_mullow(b, b, a, n, prec);
           arb_poly_mullow(tailsummandb, b, finalterm, n, prec);
           arb_poly_add_series(tailsumb, tailsumb, tailsummandb, n, prec);

           //tail a-sum
           arb_poly_pow_series(a, dndpoly, y, n, prec);
           arb_poly_scalar_mul(c, tdiv2, logdnd, prec);
           arb_poly_pow_series(c, N0poly, c, n, prec);
           arb_poly_mullow(c, c, a, n, prec);
           arb_poly_div_series(a, mollabsbt, c, n, prec);

           arb_poly_sub_series(b, oneminsiga, tdiv4logN, n, prec);
           arb_poly_pow_series(b, N0poly, b, n, prec);
           arb_poly_scalar_mul(c, Npoly, dnd, prec);
           arb_poly_sub_series(e, oneminsiga, tdiv4logNdivnd, n, prec);
           arb_poly_pow_series(c, c, e, n, prec);
           arb_poly_add_series(b, b, c, n, prec);

           arb_poly_mullow(b, b, a, n, prec);
           arb_poly_mullow(tailsummanda, b, finalterm, n, prec);
           arb_poly_add_series(tailsuma, tailsuma, tailsummanda, n, prec);
       }
       arb_poly_mullow(tailsuma, tailsuma, modgam, n, prec);

    arb_poly_add_series(sumndb, sumndb, tailsumb, n, prec);
    arb_poly_add_series(sumnda, sumnda, tailsuma, n, prec);

    //perform final calculations to complete lbound
    arb_poly_one(a);
    arb_poly_sub_series(a, a, modgam, n, prec);
    arb_poly_add_series(b, sumnda, sumndb, n, prec);
    arb_poly_sub_series(a, a, b, n, prec);
    arb_poly_div_series(tbound, a, modmol, n, prec);
    arb_poly_mullow(b, modgam, sum2, n, prec);
    arb_poly_sub_series(tbound, tbound, b, n, prec);

    arb_clear(dnd);
    arb_clear(logk);
    arb_clear(logj);
    arb_clear(logN);
    arb_clear(lognd);
    arb_clear(logdnd);
 
    arb_poly_clear(sumnda);
    arb_poly_clear(summandnda);
    arb_poly_clear(ndsum_a);
    arb_poly_clear(ndsummand_a);
    arb_poly_clear(sumndb);
    arb_poly_clear(summandndb);
    arb_poly_clear(ndsum_b);
    arb_poly_clear(ndsummand_b);

    arb_poly_clear(tailsuma);
    arb_poly_clear(tailsummanda);
    arb_poly_clear(tailsumb);
    arb_poly_clear(tailsummandb);

    arb_poly_clear(tdiv2);
    arb_poly_clear(tdiv4);
    arb_poly_clear(tdiv4logN);
    arb_poly_clear(oneminsigb);
    arb_poly_clear(oneminsiga);
    arb_poly_clear(tdiv4logNdivnd);
    arb_poly_clear(mollabsbt);
    arb_poly_clear(finalterm);

    arb_poly_clear(xN);
    arb_poly_clear(kterm);
    arb_poly_clear(kterm1);
    arb_poly_clear(sumnd);
    arb_poly_clear(ndsummand);
    arb_poly_clear(summandnd);
    arb_poly_clear(sum2);
    arb_poly_clear(modmol);
    arb_poly_clear(modgam);
    arb_poly_clear(ndterm);
    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(e);
    arb_poly_clear(btpoly);
    arb_poly_clear(atpoly);
    arb_poly_clear(dndpoly);
    arb_poly_clear(polylogN);
    arb_poly_clear(polylogdnd);
    arb_poly_clear(Npoly);
    arb_poly_clear(N0poly);
    arb_poly_clear(modK);
    arb_poly_clear(sig);
}

//main procedure to establish the Lemma-bound  
void
arb_poly_abbeff_lbound(arb_poly_t lbound, const arb_poly_t t, const arb_poly_t y,
        arb_poly_srcptr moll, slong nprimes, const ulong *primes, const slong *d,
        slong D, slong N, slong tail, slong n, slong prec)
{
    slong divs, Nend;

    arb_poly_t a, b, c, e, atpoly, btpoly, dndpoly, polylogN, polylogdnd, Npoly;
    arb_poly_t ndsum_a, ndsummand_a, tailsuma, tailsummanda, ndsummand;
    arb_poly_t ndsum_b, ndsummand_b, tailsumb, tailsummandb, sumnd, summandnd, sum2;
    arb_poly_t tdiv2, tdiv4, tdiv4logN, oneminsigb, oneminsiga, tdiv4logNdivnd, mollabsbt, finalterm;
    arb_poly_t sig, modgam, modK, xN, kterm, kterm1, ndterm, modmol;

    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(e);
    arb_poly_init(btpoly);
    arb_poly_init(atpoly);
    arb_poly_init(dndpoly);
    arb_poly_init(Npoly);
    arb_poly_init(polylogN);
    arb_poly_init(polylogdnd);
    arb_poly_init(modK);
    arb_poly_init(sumnd);
    arb_poly_init(ndsummand);
    arb_poly_init(summandnd);
    arb_poly_init(sum2);
	
    arb_poly_init(ndsum_a);
    arb_poly_init(ndsummand_a);
    arb_poly_init(ndsum_b);
    arb_poly_init(ndsummand_b);

    arb_poly_init(tailsuma);
    arb_poly_init(tailsummanda);
    arb_poly_init(tailsumb);
    arb_poly_init(tailsummandb);

    arb_poly_init(tdiv2);
    arb_poly_init(tdiv4);
    arb_poly_init(tdiv4logN);
    arb_poly_init(oneminsigb);
    arb_poly_init(oneminsiga);
    arb_poly_init(tdiv4logNdivnd);
    arb_poly_init(mollabsbt);
    arb_poly_init(finalterm);

    arb_poly_init(sig);
    arb_poly_init(modgam);
    arb_poly_init(xN);
    arb_poly_init(kterm);
    arb_poly_init(kterm1);
    arb_poly_init(ndterm);
    arb_poly_init(modmol);

    slong j, k, nd;
 
    divs = 1 << nprimes;
 
    arb_t dnd, logk, logj, logN, lognd, logdnd;
    arb_init(dnd);
    arb_init(logk);
    arb_init(logj);
    arb_init(logN);
    arb_init(lognd);
    arb_init(logdnd);

    //define common
    arb_poly_set_si(Npoly, N);

    _arb_poly_xNpoly_series(xN, Npoly, t, n, prec);
    _arb_poly_sig_series(sig, xN, t, y, n, prec);
    _arb_poly_modgamma_series(modgam, xN, y, n, prec);

    if (tail == 1)
    {
       Nend = N-1;
    }
    else
    {
       Nend = D*N-1; 
    }

    //main loop from n=1 to Nend
    arb_poly_zero(sumnd);
    for (k = 1; k <= Nend; k++)
    {
        arb_poly_zero(ndsum_a);
        arb_poly_zero(ndsum_b);
        //sub loop from nd=1 to ndivs
        for (nd = 1; nd <= divs; nd++)
        {
            if (((k+1) % d[nd] == 0) && ((k+1) <= d[nd]*N))
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

            lbound_numerator(summandnd, ndsum_b, ndsum_a, modgam, n, prec);

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
    for (k = 1; k <= N; k++)
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

    arb_poly_zero(tailsuma);
    arb_poly_zero(tailsumb);
    //calculate tailsums b and a
    if (tail == 1)
    {
       arb_set_si(logN, N);
       arb_log(logN, logN, prec);
       arb_poly_set_arb(polylogN, logN);
       arb_poly_scalar_mul_2exp_si(tdiv2, t, -1);
       arb_poly_scalar_mul_2exp_si(tdiv4, t, -2);
       arb_poly_mullow(tdiv4logN, tdiv4, polylogN, n, prec);

       arb_poly_add_si(b, y, 1, prec);
       arb_poly_scalar_mul_2exp_si(b, b, -1);
       arb_poly_neg(b, b);
       arb_poly_add_si(oneminsigb, b, 1, prec);

       arb_poly_neg(a, y);
       arb_poly_add_si(b, a, 1, prec);
       arb_poly_scalar_mul_2exp_si(b, b, -1);
       arb_poly_neg(b, b);
       arb_poly_add_si(oneminsiga, b, 1, prec);

       arb_poly_zero(tailsuma);
       arb_poly_zero(tailsumb);
       for (nd = 2; nd <= divs; nd++)
       {
           //common terms
           arb_set_si(dnd, d[nd]);
           arb_log(logdnd, dnd, prec);
           arb_poly_set_arb(polylogdnd, logdnd);
           arb_poly_set_arb(dndpoly, dnd);
           arb_poly_sub_series(a, polylogN, polylogdnd, n, prec);
           arb_poly_mullow(tdiv4logNdivnd, tdiv4, a, n, prec);

           arb_poly_abs_series(a, moll + nd, n, prec);
           _arb_poly_bt_series(btpoly, logdnd, t, n, prec);
           arb_poly_mullow(mollabsbt, a, btpoly, n, prec);

           arb_poly_scalar_div(a, polylogdnd, dnd, prec);
           arb_poly_scalar_mul_2exp_si(finalterm, a, -1);

           //tail b-sum
           arb_poly_scalar_mul(c, tdiv2, logdnd, prec);
           arb_poly_pow_series(c, Npoly, c, n, prec);
           arb_poly_div_series(a, mollabsbt, c, n, prec);

           arb_poly_sub_series(b, oneminsigb, tdiv4logN, n, prec);
           arb_poly_pow_series(b, Npoly, b, n, prec);
           arb_poly_scalar_mul(c, Npoly, dnd, prec);
           arb_poly_sub_series(e, oneminsigb, tdiv4logNdivnd, n, prec);
           arb_poly_pow_series(c, c, e, n, prec);
           arb_poly_add_series(b, b, c, n, prec);

           arb_poly_mullow(b, b, a, n, prec);
           arb_poly_mullow(tailsummandb, b, finalterm, n, prec);
           arb_poly_add_series(tailsumb, tailsumb, tailsummandb, n, prec);

           //tail a-sum
           arb_poly_pow_series(a, dndpoly, y, n, prec);
           arb_poly_scalar_mul(c, tdiv2, logdnd, prec);
           arb_poly_pow_series(c, Npoly, c, n, prec);
           arb_poly_mullow(c, c, a, n, prec);
           arb_poly_div_series(a, mollabsbt, c, n, prec);

           arb_poly_sub_series(b, oneminsiga, tdiv4logN, n, prec);
           arb_poly_pow_series(b, Npoly, b, n, prec);
           arb_poly_scalar_mul(c, Npoly, dnd, prec);
           arb_poly_sub_series(e, oneminsiga, tdiv4logNdivnd, n, prec);
           arb_poly_pow_series(c, c, e, n, prec);
           arb_poly_add_series(b, b, c, n, prec);

           arb_poly_mullow(b, b, a, n, prec);
           arb_poly_mullow(tailsummanda, b, finalterm, n, prec);
           arb_poly_add_series(tailsuma, tailsuma, tailsummanda, n, prec);
       }
       arb_poly_mullow(tailsuma, tailsuma, modgam, n, prec);
    }

    //perform final calculations to complete lbound
    arb_poly_one(a);
    arb_poly_sub_series(a, a, modgam, n, prec);
    arb_poly_sub_series(a, a, sumnd, n, prec);
    arb_poly_sub_series(a, a, tailsumb, n, prec);
    arb_poly_sub_series(a, a, tailsuma, n, prec);
    arb_poly_div_series(a, a, modmol, n, prec);
    arb_poly_mullow(b, modgam, sum2, n, prec);
    arb_poly_sub_series(lbound, a, b, n, prec);
 
    arb_clear(dnd);
    arb_clear(logk);
    arb_clear(logj);
    arb_clear(logN);
    arb_clear(lognd);
    arb_clear(logdnd);
 
    arb_poly_clear(ndsum_a);
    arb_poly_clear(ndsummand_a);
    arb_poly_clear(ndsum_b);
    arb_poly_clear(ndsummand_b);

    arb_poly_clear(tailsuma);
    arb_poly_clear(tailsummanda);
    arb_poly_clear(tailsumb);
    arb_poly_clear(tailsummandb);

    arb_poly_clear(tdiv2);
    arb_poly_clear(tdiv4);
    arb_poly_clear(tdiv4logN);
    arb_poly_clear(oneminsigb);
    arb_poly_clear(oneminsiga);
    arb_poly_clear(tdiv4logNdivnd);
    arb_poly_clear(mollabsbt);
    arb_poly_clear(finalterm);

    arb_poly_clear(xN);
    arb_poly_clear(kterm);
    arb_poly_clear(kterm1);
    arb_poly_clear(sumnd);
    arb_poly_clear(ndsummand);
    arb_poly_clear(summandnd);
    arb_poly_clear(sum2);
    arb_poly_clear(modmol);
    arb_poly_clear(modgam);
    arb_poly_clear(ndterm);
    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(e);
    arb_poly_clear(btpoly);
    arb_poly_clear(atpoly);
    arb_poly_clear(dndpoly);
    arb_poly_clear(polylogN);
    arb_poly_clear(polylogdnd);
    arb_poly_clear(Npoly);
    arb_poly_clear(modK);
    arb_poly_clear(sig);
}

//main procedure to establish the Triangle-bound  
void
arb_poly_abbeff_tbound(arb_poly_t tbound, const arb_poly_t t, const arb_poly_t y,
        arb_poly_srcptr moll, slong nprimes, const ulong *primes, const slong *d,
        slong D, slong N, slong tail, slong n, slong prec)
{
    slong divs, Nend, j, k, nd;

    arb_poly_t a, b, c, e, atpoly, btpoly, dndpoly, polylogN, polylogdnd, Npoly;
    arb_poly_t ndsum_a, ndsummand_a, tailsuma, tailsummanda, sumndb, summandndb, ndsummand;
    arb_poly_t ndsum_b, ndsummand_b, tailsumb, tailsummandb, sumnda, summandnda,sumnd, summandnd, sum2;
    arb_poly_t tdiv2, tdiv4, tdiv4logN, oneminsigb, oneminsiga, tdiv4logNdivnd, mollabsbt, finalterm;
    arb_poly_t sig, modgam, modK, xN, kterm, kterm1, ndterm, modmol;

    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    arb_poly_init(e);
    arb_poly_init(btpoly);
    arb_poly_init(atpoly);
    arb_poly_init(dndpoly);
    arb_poly_init(Npoly);
    arb_poly_init(polylogN);
    arb_poly_init(polylogdnd);
    arb_poly_init(modK);
    arb_poly_init(sumnd);
    arb_poly_init(ndsummand);
    arb_poly_init(summandnd);
    arb_poly_init(sum2);
	
    arb_poly_init(sumnda);
    arb_poly_init(summandnda);
    arb_poly_init(ndsum_a);
    arb_poly_init(ndsummand_a);
    arb_poly_init(sumndb);
    arb_poly_init(summandndb);
    arb_poly_init(ndsum_b);
    arb_poly_init(ndsummand_b);

    arb_poly_init(tailsuma);
    arb_poly_init(tailsummanda);
    arb_poly_init(tailsumb);
    arb_poly_init(tailsummandb);

    arb_poly_init(tdiv2);
    arb_poly_init(tdiv4);
    arb_poly_init(tdiv4logN);
    arb_poly_init(oneminsigb);
    arb_poly_init(oneminsiga);
    arb_poly_init(tdiv4logNdivnd);
    arb_poly_init(mollabsbt);
    arb_poly_init(finalterm);

    arb_poly_init(sig);
    arb_poly_init(modgam);
    arb_poly_init(xN);
    arb_poly_init(kterm);
    arb_poly_init(kterm1);
    arb_poly_init(ndterm);
    arb_poly_init(modmol);
 
    divs = 1 << nprimes;
 
    arb_t dnd, logk, logj, logN, lognd, logdnd;
    arb_init(dnd);
    arb_init(logk);
    arb_init(logj);
    arb_init(logN);
    arb_init(lognd);
    arb_init(logdnd);

    //define common
    arb_poly_set_si(Npoly, N);

    _arb_poly_xNpoly_series(xN, Npoly, t, n, prec);
    _arb_poly_sig_series(sig, xN, t, y, n, prec);
    _arb_poly_modgamma_series(modgam, xN, y, n, prec);

    if (tail == 1)
    {
       Nend = N;
    }
    else
    {
       Nend = D*N; 
    }

    //main loop from n=1 to D*N
    arb_poly_zero(sumnda);
    arb_poly_zero(sumndb);
    for (k = 2; k <= Nend; k++)
    {
        arb_poly_zero(ndsum_a);
        arb_poly_zero(ndsum_b);
        //sub loop from nd=1 to ndivs
        for (nd = 1; nd <= divs; nd++)
        {
            if (((k) % d[nd] == 0) && ((k) <= d[nd]*N))
            {
                j = (k) / d[nd];
 
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
            arb_set_si(logk, k);
            arb_log(logk, logk, prec);
            arb_poly_scalar_mul(kterm, sig, logk, prec);
            arb_poly_neg(kterm, kterm);
            arb_poly_exp_series(kterm, kterm, n, prec);

            arb_poly_mullow(summandndb, kterm, ndsum_b, n, prec);
            arb_poly_abs_series(summandndb, summandndb, n, prec);
            arb_poly_add_series(sumndb, sumndb, summandndb, n, prec);

            arb_poly_mullow(summandnda, kterm, ndsum_a, n, prec);
            arb_poly_abs_series(summandnda, summandnda, n, prec);
            arb_poly_add_series(sumnda, sumnda, summandnda, n, prec);
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
    for (k = 1; k <= N; k++)
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

    arb_poly_zero(tailsuma);
    arb_poly_zero(tailsumb);
    //calculate tailsums b and a
    if (tail == 1)
    {
       arb_set_si(logN, N);
       arb_log(logN, logN, prec);
       arb_poly_set_arb(polylogN, logN);
       arb_poly_scalar_mul_2exp_si(tdiv2, t, -1);
       arb_poly_scalar_mul_2exp_si(tdiv4, t, -2);
       arb_poly_mullow(tdiv4logN, tdiv4, polylogN, n, prec);

       arb_poly_add_si(b, y, 1, prec);
       arb_poly_scalar_mul_2exp_si(b, b, -1);
       arb_poly_neg(b, b);
       arb_poly_add_si(oneminsigb, b, 1, prec);

       arb_poly_neg(a, y);
       arb_poly_add_si(b, a, 1, prec);
       arb_poly_scalar_mul_2exp_si(b, b, -1);
       arb_poly_neg(b, b);
       arb_poly_add_si(oneminsiga, b, 1, prec);

       arb_poly_zero(tailsuma);
       arb_poly_zero(tailsumb);
       for (nd = 2; nd <= divs; nd++)
       {
           //common terms
           arb_set_si(dnd, d[nd]);
           arb_log(logdnd, dnd, prec);
           arb_poly_set_arb(polylogdnd, logdnd);
           arb_poly_set_arb(dndpoly, dnd);
           arb_poly_sub_series(a, polylogN, polylogdnd, n, prec);
           arb_poly_mullow(tdiv4logNdivnd, tdiv4, a, n, prec);

           arb_poly_abs_series(a, moll + nd, n, prec);
           _arb_poly_bt_series(btpoly, logdnd, t, n, prec);
           arb_poly_mullow(mollabsbt, a, btpoly, n, prec);

           arb_poly_scalar_div(a, polylogdnd, dnd, prec);
           arb_poly_scalar_mul_2exp_si(finalterm, a, -1);

           //tail b-sum
           arb_poly_scalar_mul(c, tdiv2, logdnd, prec);
           arb_poly_pow_series(c, Npoly, c, n, prec);
           arb_poly_div_series(a, mollabsbt, c, n, prec);

           arb_poly_sub_series(b, oneminsigb, tdiv4logN, n, prec);
           arb_poly_pow_series(b, Npoly, b, n, prec);
           arb_poly_scalar_mul(c, Npoly, dnd, prec);
           arb_poly_sub_series(e, oneminsigb, tdiv4logNdivnd, n, prec);
           arb_poly_pow_series(c, c, e, n, prec);
           arb_poly_add_series(b, b, c, n, prec);

           arb_poly_mullow(b, b, a, n, prec);
           arb_poly_mullow(tailsummandb, b, finalterm, n, prec);
           arb_poly_add_series(tailsumb, tailsumb, tailsummandb, n, prec);

           //tail a-sum
           arb_poly_pow_series(a, dndpoly, y, n, prec);
           arb_poly_scalar_mul(c, tdiv2, logdnd, prec);
           arb_poly_pow_series(c, Npoly, c, n, prec);
           arb_poly_mullow(c, c, a, n, prec);
           arb_poly_div_series(a, mollabsbt, c, n, prec);

           arb_poly_sub_series(b, oneminsiga, tdiv4logN, n, prec);
           arb_poly_pow_series(b, Npoly, b, n, prec);
           arb_poly_scalar_mul(c, Npoly, dnd, prec);
           arb_poly_sub_series(e, oneminsiga, tdiv4logNdivnd, n, prec);
           arb_poly_pow_series(c, c, e, n, prec);
           arb_poly_add_series(b, b, c, n, prec);

           arb_poly_mullow(b, b, a, n, prec);
           arb_poly_mullow(tailsummanda, b, finalterm, n, prec);
           arb_poly_add_series(tailsuma, tailsuma, tailsummanda, n, prec);
       }
       arb_poly_mullow(tailsuma, tailsuma, modgam, n, prec);
    }

    arb_poly_add_series(sumndb, sumndb, tailsumb, n, prec);
    arb_poly_add_series(sumnda, sumnda, tailsuma, n, prec);

    //perform final calculations to complete lbound
    arb_poly_one(a);
    arb_poly_sub_series(a, a, modgam, n, prec);
    arb_poly_add_series(b, sumnda, sumndb, n, prec);
    arb_poly_sub_series(a, a, b, n, prec);
    arb_poly_div_series(tbound, a, modmol, n, prec);
    arb_poly_mullow(b, modgam, sum2, n, prec);
    arb_poly_sub_series(tbound, tbound, b, n, prec);

    arb_clear(dnd);
    arb_clear(logk);
    arb_clear(logj);
    arb_clear(logN);
    arb_clear(lognd);
    arb_clear(logdnd);
 
    arb_poly_clear(sumnda);
    arb_poly_clear(summandnda);
    arb_poly_clear(ndsum_a);
    arb_poly_clear(ndsummand_a);
    arb_poly_clear(sumndb);
    arb_poly_clear(summandndb);
    arb_poly_clear(ndsum_b);
    arb_poly_clear(ndsummand_b);

    arb_poly_clear(tailsuma);
    arb_poly_clear(tailsummanda);
    arb_poly_clear(tailsumb);
    arb_poly_clear(tailsummandb);

    arb_poly_clear(tdiv2);
    arb_poly_clear(tdiv4);
    arb_poly_clear(tdiv4logN);
    arb_poly_clear(oneminsigb);
    arb_poly_clear(oneminsiga);
    arb_poly_clear(tdiv4logNdivnd);
    arb_poly_clear(mollabsbt);
    arb_poly_clear(finalterm);

    arb_poly_clear(xN);
    arb_poly_clear(kterm);
    arb_poly_clear(kterm1);
    arb_poly_clear(sumnd);
    arb_poly_clear(ndsummand);
    arb_poly_clear(summandnd);
    arb_poly_clear(sum2);
    arb_poly_clear(modmol);
    arb_poly_clear(modgam);
    arb_poly_clear(ndterm);
    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    arb_poly_clear(e);
    arb_poly_clear(btpoly);
    arb_poly_clear(atpoly);
    arb_poly_clear(dndpoly);
    arb_poly_clear(polylogN);
    arb_poly_clear(polylogdnd);
    arb_poly_clear(Npoly);
    arb_poly_clear(modK);
    arb_poly_clear(sig);
}

//initialise and fill arrays for bt, at, moll. Process main loop  
void static
set_bound_constants(slong res, const arb_t t, const arb_t y, slong N, slong N0, const slong type, const slong tail,
          slong nprimes, slong digits, slong n, slong prec)
{

    slong D;

    slong i, j;

    arb_poly_ptr moll;
    slong divs;
    ulong *primes;
    slong *d;

    divs=0; D=0; 
 
    arb_t logp, coeff;
    arb_init(logp);
    arb_init(coeff);

    arb_poly_t tpoly, ypoly, btpoly, bound;
    arb_poly_init(tpoly);
    arb_poly_init(ypoly);
    arb_poly_init(btpoly);
    arb_poly_init(bound);
 
    arb_poly_set_arb(tpoly, t);
    arb_poly_set_arb(ypoly, y);
 
    divs = 1 << nprimes;
    primes = flint_malloc(nprimes * sizeof(*primes));
    _first_n_primes(primes, nprimes);
 
    D = 1;
    for (i = 0; i < nprimes; i++)
    {
        D *= primes[i];
    }
 
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

    if (type == 0)
    {
       arb_poly_abbeff_tbound(bound, tpoly, ypoly, moll, nprimes, primes, d, D, N, tail, 1, prec);
       printf("\n"); 
       printf("Triangle bound ");
       if (tail == 1)
       {
          printf("(with tail sum) ");
       }
       else
       {
          printf("(no tail sum) ");
       }
       printf("for t=");
       arb_printn(t, 3, ARB_STR_NO_RADIUS);
       printf(", y=");
       arb_printn(y, 3, ARB_STR_NO_RADIUS);
       printf(", N=%ld, ", N);
       printf("using %ld primes:", nprimes);
       printf("\n");
       printf("\n");
       arb_poly_get_coeff_arb(coeff, bound, 0);
       arf_printd(arb_midref(coeff), digits);
       printf("\n");
       printf("\n");
    }
    if (type == 1)
    {
       arb_poly_abbeff_lbound(bound, tpoly, ypoly, moll, nprimes, primes, d, D, N, tail, 1, prec);
       printf("\n");
       printf("Lemma bound ");
       if (tail == 1)
       {
          printf("(with tail sum) ");
       }
       else
       {
          printf("(no tail sum) ");
       }
       printf("for t=");
       arb_printn(t, 3, ARB_STR_NO_RADIUS);
       printf(", y=");
       arb_printn(y, 3, ARB_STR_NO_RADIUS);
       printf(", N=%ld, ", N);
       printf("using %ld primes:", nprimes);
       printf("\n");
       printf("\n");
       arb_poly_get_coeff_arb(coeff, bound, 0);
       arf_printd(arb_midref(coeff), digits);
       printf("\n");
       printf("\n");
    }

    if (type == 2)
    {
       nprimes = 0;
       arb_poly_abbeff_approx_tbound_moll0(bound, tpoly, ypoly, moll, nprimes, primes, d, D, N, N0, 1, prec);
       printf("\n");
       printf("Approximate Triangle bound (with tail sum) moll(0) ");
       printf("for t=");
       arb_printn(t, 3, ARB_STR_NO_RADIUS);
       printf(", y=");
       arb_printn(y, 3, ARB_STR_NO_RADIUS);
       printf(", N=%ld, ", N);
       printf("N0=%ld, ", N0);
       printf("using %ld primes:", nprimes);
       printf("\n");
       printf("\n");
       arb_poly_get_coeff_arb(coeff, bound, 0);
       arf_printd(arb_midref(coeff), digits);
       printf("\n");
       printf("\n");
    }

    if (type == 3)
    {
       arb_poly_abbeff_approx_tbound_moll2(bound, tpoly, ypoly, N, N0, 1, prec);
       printf("\n");
       printf("Approximate Triangle bound (with tail integral) moll(2) ");
       printf("for t=");
       arb_printn(t, 3, ARB_STR_NO_RADIUS);
       printf(", y=");
       arb_printn(y, 3, ARB_STR_NO_RADIUS);
       printf(", N=%ld, ", N);
       printf("N0=%ld, ", N0);
       printf("using 1 prime:");
       printf("\n");
       printf("\n");
       arb_poly_get_coeff_arb(coeff, bound, 0);
       arf_printd(arb_midref(coeff), digits);
       printf("\n");
       printf("\n");
    }

    //clear all variables 
    for (i = 0; i < divs; i++)
    {
        arb_poly_clear(moll + i + 1);
    }
    flint_free(d);
    flint_free(moll);
    flint_free(primes);

    arb_clear(logp);
    arb_clear(coeff);

    arb_poly_clear(tpoly);
    arb_poly_clear(ypoly);
    arb_poly_clear(btpoly);
    arb_poly_clear(bound);
}
 
//main initialisation - processing and validating input values 
int main(int argc, char *argv[])
{
    const char *t_str, *y_str, *N_str, *N0_str, *type_str, *tail_str, *m_str, *d_str;

    slong N, N0, type, tail, m, d, prec, res;

    int result = EXIT_SUCCESS;
    int usage_error = 0;
    res = 0;
 
    arb_t t, y;
    arb_init(t);
    arb_init(y);
 
    if (argc != 9)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }
 
    t_str = argv[1];
    y_str = argv[2];
    N_str = argv[3];
    N0_str = argv[4];
    type_str = argv[5];
    tail_str = argv[6];
    m_str = argv[7];
    d_str = argv[8];
 
    N = atol(N_str);
    N0 = atol(N0_str);
    type = atol(type_str);
    tail = atol(tail_str);
    m = atol(m_str);
    d = atol(d_str);

    if (N < 1 || N0 < 0 || N0 > N || m < 0 || m > 20 || d < 1 || tail < 0 || tail > 1 || type < 0 || type > 3)
    {
        usage_error = 1;
        result = EXIT_FAILURE;
        goto finish;
    }
 
    prec = d * 3.32192809488736 + 30;

    arb_set_str(t, t_str, prec);
    arb_set_str(y, y_str, prec);

TIMEIT_ONCE_START

        set_bound_constants(res, t, y, N, N0, type, tail, m, d, 1, prec);

TIMEIT_ONCE_STOP 
 
finish:
 
    if (usage_error && result == EXIT_FAILURE)
    {
        flint_printf("Usage:\n");
        flint_printf("%s t y N N0 type tail m d\n\n", argv[0]);
        flint_printf(
    "This script computes the lower bound of abs(Aeff + Beff)/abs(Beff0)\n"
    "for t,y,N. It uses an Euler mollification that includes m primes,\n"
    "so that for example m=0 corresponds to no mollification\n"
    "and m=3 corresponds to an Euler mollification that uses\n"
    "the first three primes {2, 3, 5}.\n"
    "Calculations can be done with a tailsum (tail=1) or without (tail=0).\n"
    "The 'type' of bound can be selected:\n"
    " 0 = Triangle lower bound (N0 is ignored);\n" 
    " 1 = Lemma lower bound (N0 is ignored);\n"
    " 2 = Approximate Triangle lower bound for N, N0 (no mollifier).\n"
    " 3 = Approximate Triangle lower bound for N, N0 (moll(2)).\n"
    "The number of decimal digits to be displayed can be set by d."
    "\n");
    }
 
    arb_clear(t);
    arb_clear(y);;
 
    flint_cleanup();
    return result;
}
