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
#include "acb_calc.h"
#include "flint/profiler.h"
#include "pthread.h"

//define some global parameters to be accessed by individual threads
acb_mat_t ests;
arb_mat_t minmodabb;

arb_t x;
arb_t y;
arb_t t;
arb_t n0;
arb_t logtn0;
acb_t ssexpo;
acb_poly_t finpoly;

struct ThreadData {
slong start, stop, id, prec, num, prt;
}; 

// prepare the function as integrand for the integral
int
f_ddzbound(acb_ptr res, const acb_t z, void *param, slong order, slong prec)
{
    arb_t a, b, c, u, pi, Nu, logNu;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(u);
    arb_init(Nu);
    arb_init(logNu);
    arb_init(pi);
    arb_const_pi(pi, prec);    

    arb_ptr(Narb);
    arb_ptr(y);
    arb_ptr(x);
    arb_ptr(t);
    arb_ptr(tdiv4divxmin6);
    arb_ptr(afac);
    arb_ptr(const1);
    arb_ptr(estinit);

    Narb = ((arb_ptr *)(param))[0];
    y = ((arb_ptr *)(param))[1];
    x = ((arb_ptr *)(param))[2];
    t = ((arb_ptr *)(param))[3];
    tdiv4divxmin6 = ((arb_ptr *)(param))[4];
    afac = ((arb_ptr *)(param))[5];
    const1 = ((arb_ptr *)(param))[6];
    estinit = ((arb_ptr *)(param))[7];

    acb_get_real(u, z);
    arb_mul(Nu, u, Narb, prec);
    arb_log(logNu, Nu, prec);

    //establish the integrand
    arb_pow(a, Nu, y, prec);
    arb_mul(b, afac, a, prec); 
    arb_add_si(b, b, 1, prec);
    arb_mul(b, b, tdiv4divxmin6, prec);
    arb_one(c);
    arb_mul_2exp_si(c, c, -1); 
    arb_add(b, b, c, prec);
    arb_mul(b, logNu, b, prec);
    arb_mul(c, estinit, a, prec);
    arb_add(c, b, c, prec);
    arb_neg(a, y);
    arb_sub_si(a, a, 1, prec);
    arb_mul_2exp_si(a, a, -1);                
    arb_mul_2exp_si(b, logNu, -2);
    arb_sub(b, b, const1, prec);
    arb_mul(b, b, t, prec);
    arb_add(a, a, b, prec);
    arb_pow(a, Nu, a, prec);       
    arb_mul(a, c, a, prec); 
    acb_set_arb(res, a);

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(u);
    arb_clear(Nu);
    arb_clear(pi);
    arb_clear(logNu);

    return 0;
}

void
generate_ddzbound(arb_t esta, const arb_t x, const arb_t y, const arb_t t, arb_t Narb, arb_t afac, arb_t const1,
				arb_t tdiv4divxmin6, slong prec)

{  
    acb_t tmp, ai, bi, est;
    acb_init(ai);
    acb_init(bi);
    acb_init(tmp);
    acb_init(est);

    arb_t a, b, c, d, e, estinit, pi;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(e);
    arb_init(estinit);
    arb_init(pi);
    arb_const_pi(pi, prec);    

    slong goal;

    //estinit
    arb_add_si(b, y, 1, prec);
    acb_set_arb_arb(tmp, b, x);
    acb_abs(b, tmp, prec);
    arb_div(b, b, pi, prec);
    arb_mul_2exp_si(b, b, -2);
    arb_log(b, b, prec); 
    arb_add(b, b, pi, prec); 
    arb_set_si(c, 3);
    arb_div(c, c, x, prec);
    arb_add(b, b, c, prec);
    arb_mul(b, b, afac, prec);
    arb_one(d);
    arb_mul_2exp_si(d, d, -1);
    arb_add(d, tdiv4divxmin6, d, prec);
    arb_mul(estinit, b, d, prec);

    //evaluate the integral
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

    void *param[8];

    param[0] = (void *) Narb;
    param[1] = (void *) y;
    param[2] = (void *) x;
    param[3] = (void *) t;
    param[4] = (void *) tdiv4divxmin6;
    param[5] = (void *) afac;
    param[6] = (void *) const1;
    param[7] = (void *) estinit;

    acb_calc_integrate(est, f_ddzbound, param, ai, bi, goal, tol, NULL, prec);
    acb_abs(esta, est, prec);
    arb_mul(esta, esta, Narb, prec);
    arb_add(esta, esta, estinit, prec);

    mag_clear(tol);

    acb_clear(tmp);
    acb_clear(ai);
    acb_clear(bi);
    acb_clear(est);

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(e);
    arb_clear(estinit);
    arb_clear(pi);
}

// prepare the function as integrand for the integral
int
f_ddtbound(acb_ptr res, const acb_t z, void *param, slong order, slong prec)
{
    arb_t a, b, c, d, u, pi, Nu, logNu;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(u);
    arb_init(Nu);
    arb_init(logNu);
    arb_init(pi);
    arb_const_pi(pi, prec);    

    arb_ptr(Narb);
    arb_ptr(y);
    arb_ptr(t);
    arb_ptr(afac);
    arb_ptr(const1);
    arb_ptr(const2);
    arb_ptr(estinit);

    Narb = ((arb_ptr *)(param))[0];
    y = ((arb_ptr *)(param))[1];
    t = ((arb_ptr *)(param))[2];
    afac = ((arb_ptr *)(param))[3];
    const1 = ((arb_ptr *)(param))[4];
    const2 = ((arb_ptr *)(param))[5];
    estinit = ((arb_ptr *)(param))[6];

    acb_get_real(u, z);

    arb_mul(Nu, u, Narb, prec);
    arb_log(logNu, Nu, prec);

    //establish the integrand
    arb_mul_2exp_si(a, logNu, -2); 
    arb_sub(a, const2, a, prec);
    arb_mul(a, a, logNu, prec); 
    arb_pow(d, Nu, y, prec);
    arb_mul(b, d, afac, prec); 
    arb_add_si(b, b, 1, prec); 
    arb_mul(a, b, a, prec); 
    arb_mul(c, estinit, d, prec);
    arb_add(c, c, a, prec);
    arb_neg(a, y);
    arb_sub_si(a, a, 1, prec);
    arb_mul_2exp_si(a, a, -1);  
    arb_mul_2exp_si(b, logNu, -2);
    arb_sub(b, b, const1, prec);
    arb_mul(b, b, t, prec);
    arb_add(a, a, b, prec);
    arb_pow(a, Nu, a, prec);       
    arb_mul(a, c, a, prec); 
    acb_set_arb(res, a);

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(u);
    arb_clear(Nu);
    arb_clear(pi);
    arb_clear(logNu);

    return 0;
}

void
generate_ddtbound(arb_t esta, const arb_t x, const arb_t y, const arb_t t, arb_t Narb, arb_t afac, arb_t const1,
                  arb_t xdiv4pi, arb_t logxdiv4pi, arb_t onedivxmin6, slong prec)

{  
    arb_t a, b, c, d, const2, estinit, pi;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(const2);
    arb_init(estinit);
    arb_init(pi);
    arb_const_pi(pi, prec);    

    acb_t ai, bi, est;
    acb_init(ai);
    acb_init(bi);
    acb_init(est);

    slong goal;

    //const2
    arb_mul_2exp_si(a, logxdiv4pi, -2);
    arb_mul_2exp_si(b, pi, -3);
    arb_add(a, a, b, prec);
    arb_mul_2exp_si(b, onedivxmin6, 1);
    arb_add(const2, a, b, prec);

    //estinit
    arb_mul_2exp_si(a, onedivxmin6, 3);
    arb_mul_2exp_si(b, pi, -1);
    arb_add(b, a, b, prec);
    arb_add(c, logxdiv4pi, a, prec);
    arb_mul(c, b, c, prec);
    arb_mul(c, c, afac, prec);
    arb_mul_2exp_si(estinit, c, -2);

    //evaluate the integral
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

    void *param[7];

    param[0] = (void *) Narb;
    param[1] = (void *) y;
    param[2] = (void *) t;
    param[3] = (void *) afac;
    param[4] = (void *) const1;
    param[5] = (void *) const2;
    param[6] = (void *) estinit;

    acb_calc_integrate(est, f_ddtbound, param, ai, bi, goal, tol, NULL, prec);
    acb_abs(esta, est, prec);
    arb_mul(esta, esta, Narb, prec);
    arb_add(esta, esta, estinit, prec);

    mag_clear(tol);

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(const2);
    arb_clear(estinit);
    arb_clear(pi);

    acb_clear(ai);
    acb_clear(bi);
    acb_clear(est);
} 

void H01(acb_t z, const acb_t s, slong prec)
{
    acb_t a, b, c, d, half, piacb;
    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(d);
    acb_init(half);
    acb_init(piacb);

    arb_t ab, pi;
    arb_init(ab);
    arb_init(pi);
    arb_const_pi(pi, prec);  
 
    acb_set_arb(piacb, pi);
    acb_one(a);
    acb_mul_2exp_si(half, a, -1);
 
    acb_mul_2exp_si(d, s, -1);
    acb_sub_si(b, s, 1, prec);
    acb_mul(b, b, d, prec);   

    acb_neg(a, s);
    acb_mul_2exp_si(a, a, -1);
    acb_pow(a, piacb, a, prec);
    acb_mul(b, b, a, prec); 

    arb_mul_2exp_si(ab, pi, 1);
    arb_sqrt(ab, ab, prec);
    acb_mul_arb(b, b, ab, prec); 

    acb_log(a, d, prec);
    acb_sub(c, d, half, prec);
    acb_mul(c, c, a, prec);
    acb_sub(c, c, d, prec);
    acb_exp(c, c, prec);
    acb_mul(z, c, b, prec); 
 
    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(d);
    acb_clear(half);
    acb_clear(piacb);

    arb_clear(ab);
    arb_clear(pi);
}

void alpha1(acb_t z, const acb_t s, slong prec)
{
    acb_t a, b;
    acb_init(a);
    acb_init(b);

    arb_t ab, pi;
    arb_init(ab);
    arb_init(pi);
    arb_const_pi(pi, prec);  
 
    acb_mul_2exp_si(a, s, 1);
    acb_inv(a, a, prec);
 
    acb_sub_si(b, s, 1, prec);
    acb_inv(b, b, prec);   

    acb_add(a, a, b, prec);

    arb_mul_2exp_si(ab, pi , 1);
    acb_div_arb(b, s, ab, prec);
    acb_log(b, b, prec);   
    acb_mul_2exp_si(b, b, -1);

    acb_add(z, a, b, prec);
 
    acb_clear(a);
    acb_clear(b);

    arb_clear(ab);
    arb_clear(pi);
}

static void
bexpo_aexpo_afac_bsums_asums(acb_mat_t ests, acb_t sarr, const acb_t ssexpo, const acb_poly_t finpoly, const arb_t t, 
                             const arb_t logtn0, const arb_t n0, slong id, const slong k, const slong prec)

{
    acb_t a, b, c, d, s, afac, bexpo, aexpo, bsums, asums, negs, conjs, onemins, alphas, alphaconjs, alpha1mins, one;
    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(d);
    acb_init(s);
    acb_init(afac);
    acb_init(bexpo);
    acb_init(aexpo);
    acb_init(bsums);
    acb_init(asums);
    acb_init(negs);
    acb_init(conjs);
    acb_init(onemins);
    acb_init(alphas);
    acb_init(alpha1mins);
    acb_init(alphaconjs);
    acb_init(one);
    acb_one(one);

    arb_t ab, negt;
    arb_init(ab);
    arb_init(negt);

    arb_neg(negt, t);

    //set alpha1(sarr(k))/2, alpha1(1-sarr(k))/2, alpha1conj(sarr(k))/2
    acb_set(s, sarr);
    alpha1(alphas, s, prec);
    acb_mul_2exp_si(alphas, alphas, -1);

    acb_conj(conjs, s);
    alpha1(alphaconjs, conjs, prec);
    acb_mul_2exp_si(alphaconjs, alphaconjs, -1);

    acb_neg(negs, s);
    acb_add_si(onemins, negs, 1, prec);
    alpha1(alpha1mins, onemins, prec);
    acb_mul_2exp_si(alpha1mins, alpha1mins, -1);

    //bexpo
    acb_mul_arb(b, alpha1mins, negt, prec);
    acb_sub(b, b, onemins, prec);
    acb_add(bexpo, b, ssexpo, prec);

    //aexpo
    acb_mul_arb(b, alphaconjs, negt, prec);
    acb_sub(b, b, conjs, prec);
    acb_add(aexpo, b, ssexpo, prec);

    //afac
    acb_pow_ui(b, alphas, 2, prec);
    acb_pow_ui(d, alpha1mins, 2, prec);
    acb_sub(b, b, d, prec);
    acb_mul_arb(b, b, t, prec);
    acb_exp(a, b, prec);
    H01(b, s, prec);            
    H01(d, onemins, prec);   
    acb_div(b, b, d, prec);
    acb_mul(afac, b, a, prec);

    //bsums
    arb_mul_2exp_si(ab, logtn0, -1);
    acb_add_arb(a, bexpo, ab, prec);
    acb_set_arb(b, n0);
    acb_pow(a, b, a, prec);
    acb_add_arb(b, bexpo, logtn0, prec);
    acb_poly_evaluate(c, finpoly, b, prec);
    acb_mul(bsums, c, a, prec);

    //asums
    acb_add_arb(a, aexpo, ab, prec);
    acb_set_arb(b, n0);
    acb_pow(a, b, a, prec);
    acb_add_arb(b, aexpo, logtn0, prec);
    acb_poly_evaluate(c, finpoly, b, prec);
    acb_mul(asums, c, a, prec);

    //ests
    acb_conj(b, asums);
    acb_mul(a, afac, b, prec);
    acb_add(a, a, bsums, prec);
    acb_set(acb_mat_entry(ests, k, 0), a);

    //establish minimum value of modabb
    acb_abs(ab, acb_mat_entry(ests, k, 0), prec);
    if (arb_lt(ab, arb_mat_entry(minmodabb, id, 0)))
    {
          arb_set(arb_mat_entry(minmodabb, id, 0), ab);
    } 

    arb_clear(ab);
    arb_clear(negt);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(d);
    acb_clear(s);
    acb_clear(negs);
    acb_clear(conjs);
    acb_clear(onemins);
    acb_clear(alphas);
    acb_clear(alpha1mins);
    acb_clear(alphaconjs);
    acb_clear(one);
}

static void
print_details(slong res, const arb_t t, const arb_t a, const arb_t b, const acb_t c)

{
    arf_printd(arb_midref(t), 20);
    printf(", ");
    arf_printd(arb_midref(a), 30);
    printf(", ");
    arf_printd(arb_midref(b), 30);
    printf(", ");
    acb_printn(c, 30, ARB_STR_NO_RADIUS);
    printf("\n");
}

//process a rectangle for each t

void* abbeff_symmetric_rectangle(void *voidData)
{
    //recover the data passed to this specific thread
    struct ThreadData* data=voidData;
    slong start=data->start;
    slong stop=data->stop;
    slong id=data->id;
    slong prec=data->prec;
    slong num=data->num;
    slong prt=data->prt;

    acb_t a, b, c, sarr, one, onei, Idiv4, IXdiv2;
    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(sarr);
    acb_init(Idiv4);
    acb_init(IXdiv2);
    acb_init(one);
    acb_one(one);
    acb_init(onei);
    acb_onei(onei);

    arb_t ab, ac, ad, thtarrv, thtarrv1, zarrv, zarrv1, half, oneminy, oneminydiv2, ydiv2, numarb, nummin1, varb;
    arb_init(ab);
    arb_init(ac);
    arb_init(ad);
    arb_init(thtarrv);
    arb_init(thtarrv1);
    arb_init(zarrv);
    arb_init(zarrv1);
    arb_init(half);
    arb_init(numarb);
    arb_init(nummin1);
    arb_init(varb);
    arb_init(oneminy);
    arb_init(oneminydiv2);
    arb_init(ydiv2);

    slong k, v, result;
    result = 0;

    //shared values     
    acb_mul_2exp_si(Idiv4, onei, -2); 
    arb_neg(ab, y); 
    arb_add_si(oneminy, ab, 1, prec); 
    arb_mul_2exp_si(oneminydiv2, oneminy, -1);
    acb_set_arb(b, x); 
    acb_mul_onei(a, b); 
    acb_mul_2exp_si(IXdiv2, a, -1);
    arb_mul_2exp_si(ydiv2, y, -1);
    arb_one(ab);
    arb_mul_2exp_si(half, ab, -1);	
    arb_set_ui(numarb, num);
    arb_sub_si(nummin1, numarb, 1, prec); 

    arb_mul_2exp_si(ab, half, -1);
    arb_neg(thtarrv, ab);  
    arb_neg(zarrv, oneminydiv2);

    for (v = start; v < stop; v++)
    { 
       arb_set_ui(varb, v); 
       arb_add_si(varb, varb, 1, prec); 
       arb_div(ad, varb, nummin1, prec);

       //thtarr
       arb_sub(ac, ad, half, prec);
       arb_mul_2exp_si(thtarrv1, ac, -1);

       //zarr 
       arb_mul(ac, oneminydiv2, ad, prec);
       arb_mul_2exp_si(ac, ac, 1);
       arb_sub(zarrv1, ac, oneminydiv2, prec);

       //x lower constant
       k = v;
       acb_add_arb(a, IXdiv2, oneminydiv2, prec);
       acb_sub_arb(a, a, zarrv, prec);
       acb_sub(sarr, a, Idiv4, prec);
       bexpo_aexpo_afac_bsums_asums(ests, sarr, ssexpo, finpoly, t, logtn0, n0, id, k, prec);

       if (prt==1)
       {      
          arb_mul_2exp_si(ab, zarrv, 1);
          arb_add(ab, y, ab, prec); 
          acb_set_arb(acb_mat_entry(ests, k, 1), ab);
          arb_sub(ad, x, half, prec);
          acb_set_arb(acb_mat_entry(ests, k, 2), ad);
        }

       //y upper constant
       k = num+v-1;
       acb_set_arb(a, thtarrv);
       acb_mul_onei(a, a);
       acb_add(sarr, IXdiv2, a, prec);
       bexpo_aexpo_afac_bsums_asums(ests, sarr, ssexpo, finpoly, t, logtn0, n0, id, k, prec);

       if (prt==1)
       {      
          arb_set_si(ab, 1);
          acb_set_arb(acb_mat_entry(ests, k, 1), ab);
          arb_mul_2exp_si(ad, thtarrv, 1);
          arb_add(ad, x, ad, prec);
          acb_set_arb(acb_mat_entry(ests, k, 2), ad);
       }

       //x upper and output to be attached in reverse order
       k = 3*num-(v+1)-3;
       acb_add_arb(a, IXdiv2, oneminydiv2, prec);
       acb_sub_arb(a, a, zarrv1, prec);
       acb_add(sarr, a, Idiv4, prec);
       bexpo_aexpo_afac_bsums_asums(ests, sarr, ssexpo, finpoly, t, logtn0, n0, id, k, prec);

       if (prt==1)
       {  
          arb_mul_2exp_si(ab, zarrv1, 1);    
          arb_add(ab, y, ab, prec);
          acb_set_arb(acb_mat_entry(ests, k, 1), ab);
          arb_add(ad, x, half, prec);
          acb_set_arb(acb_mat_entry(ests, k, 2), ad);
       }

       //y lower and output to be attached in reverse order
       k = 4*num-(v+1)-4;
       acb_add_arb(a, IXdiv2, oneminydiv2, prec);
       acb_set_arb(b, thtarrv1);
       acb_mul_onei(b, b);
       acb_add(a, a, b, prec);
       acb_add_arb(sarr, a, oneminydiv2, prec);
       bexpo_aexpo_afac_bsums_asums(ests, sarr, ssexpo, finpoly, t, logtn0, n0, id, k, prec);

       if (prt==1)
       {      
          arb_add_si(ab, y, -1, prec);
          arb_add(ab, ab, y, prec);
          acb_set_arb(acb_mat_entry(ests, k, 1), ab);
          arb_mul_2exp_si(ad, thtarrv1, 1);
          arb_add(ad, x, ad, prec);
          acb_set_arb(acb_mat_entry(ests, k, 2), ad);
       }
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(sarr);
    acb_clear(one);
    acb_clear(onei);
    acb_clear(Idiv4);
    acb_clear(IXdiv2);

    arb_clear(ab);
    arb_clear(ac);
    arb_clear(ad);
    arb_clear(thtarrv);
    arb_clear(thtarrv1);
    arb_clear(zarrv);
    arb_clear(zarrv1);
    arb_clear(half);
    arb_clear(oneminy);
    arb_clear(oneminydiv2);
    arb_clear(ydiv2);
    arb_clear(numarb);
    arb_clear(nummin1);
    arb_clear(varb);

    return(NULL);
}

//procedure to transform the matrix into a poly
static void
mattopoly(acb_poly_t polyres, const arb_t tval, const slong expterms, const slong taylorterms, 
           const acb_mat_t mattmp, const slong prec)
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
abbeff_t_loop(slong res, const arb_t X, const arb_t y0, const arb_t ts, const arb_t te, const slong N, 
              const slong taylorterms, const slong expterms, const acb_mat_t finalmat,
              const slong prt, slong numthreads, const slong digits, const slong prec)

{
    arb_t a, b, c, d, pi, ab, bb, Narb, xdiv4pi, logxdiv4pi, onedivxmin6, tdiv4divxmin6;
    arb_t minmod, windnum, windtot, dzabb, dtabb, afac, const1;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(x);
    arb_init(y);
    arb_init(t);
    arb_init(ab);
    arb_init(bb);
    arb_init(Narb);
    arb_init(xdiv4pi);
    arb_init(logxdiv4pi);
    arb_init(onedivxmin6);
    arb_init(tdiv4divxmin6);
    arb_init(n0);
    arb_init(logtn0);
    arb_init(minmod);
    arb_init(windnum);
    arb_init(windtot);
    arb_init(dzabb);
    arb_init(dtabb);
    arb_init(afac);
    arb_init(const1);
    arb_init(pi);
    arb_const_pi(pi, prec);

    acb_t ca, argdiv, one;
    acb_init(ca);
    acb_init(argdiv);
    acb_init(ssexpo);
    acb_init(one);
    acb_one(one);

    acb_poly_init(finpoly);
   
    arb_set_si(Narb, N);

    printf("\n");
    printf("Processing the barrier for X= ");
    arf_printd(arb_midref(X), 20);
    printf("...");
    arb_add_si(a, X, 1, prec);
    arf_printd(arb_midref(a), 20);
    printf(" (N = %ld)", N);
    printf(", y0 = ");
    arf_printd(arb_midref(y0), 10);
    printf("...1 ");
    printf(", t = ");
    arf_printd(arb_midref(ts), 10);
    printf("...");
    arf_printd(arb_midref(te), 10);
    printf("\n");
    printf("\n");

    slong i, idx, count, rectmesh, num, prtresult;
    prtresult = 0;

    //change X and y to midpoints
    arb_one(b);
    arb_mul_2exp_si(a, b, -1);  
    arb_add(x, X, a, prec); 
    arb_add_si(a, y0, 1, prec);
    arb_mul_2exp_si(y, a, -1);

    //n0
    arb_mul_2exp_si(n0, Narb, -1);

    //ssexpo
    acb_set_arb(ca, x);
    acb_mul_onei(ca, ca);   
    acb_neg(ca, ca);
    acb_add_si(ca, ca, 1, prec);
    acb_mul_2exp_si(ssexpo, ca, -1);   
    
    //shared values for the ddz and ddtbounds
    arb_set_si(Narb, N);
    arb_mul_2exp_si(a, pi, 2);
    arb_div(xdiv4pi, X, a, prec);
    arb_log(logxdiv4pi, xdiv4pi, prec);
    arb_sub_si(b, X, 6, prec);
    arb_inv(onedivxmin6, b, prec);

    //const1 required for the ddz and ddtbounds
    arb_mul_2exp_si(a, logxdiv4pi, -2);
    arb_set_si(b, -3);
    arb_mul(b, b, y0, prec);
    arb_add_si(b, b, 1, prec);
    arb_mul_2exp_si(c, y0, 2);
    arb_add_si(d, y0, 1, prec);
    arb_mul(c, c, d, prec);
    arb_mul(d, X, X, prec);
    arb_div(c, c, d, prec);
    arb_add(c, c, b, prec);
    arb_zero(b);
    arb_max(c, b, c, prec);
    arb_div(c, c, d, prec);
    arb_mul_2exp_si(c, c, -1);
    arb_sub(const1, a, c, prec);

    //perform the t-loop over all the X..X+1, y0..1 rectangles
    arb_set(t, ts);
    arb_zero(windtot);
    count = 0;

    while(arb_le(t, te))
    {
        count=count+1;

        //logtn0
        arb_mul_2exp_si(a, t, -1);
        arb_log(b, n0, prec);
        arb_mul(logtn0, b, a, prec);

        mattopoly(finpoly, t, expterms, taylorterms, finalmat, prec);

        //calculate ddz and ddt bounds
        arb_mul_2exp_si(a, t, -2);
        arb_mul(tdiv4divxmin6, a, onedivxmin6, prec);

        //afac-term
        arb_set_si(a, 2);
        arb_set_si(b, 100);
        arb_div(a, a, b, prec);
        arb_mul(a, y0, a, prec);
        arb_exp(a, a, prec);
        arb_neg(b, y0);
        arb_mul_2exp_si(b, b, -1); 
        arb_pow(b, xdiv4pi, b, prec);
        arb_mul(d, a, b, prec);
        arb_mul(a, t, y0, prec);
        arb_mul(a, a, onedivxmin6, prec);
        arb_mul_2exp_si(a, a, -1);
        arb_pow(b, Narb, a, prec); 
        arb_mul(afac, d, b, prec);

        generate_ddzbound(dzabb, X, y0, t, Narb, afac, const1, tdiv4divxmin6, prec);
        generate_ddtbound(dtabb, X, y0, t, Narb, afac, const1, xdiv4pi, logxdiv4pi, onedivxmin6, prec);  
        arb_ceil(dzabb, dzabb, prec);
        arb_ceil(dtabb, dtabb, prec);

        //establish number of mesh points
        num = arf_get_d(arb_midref(dzabb), ARF_RND_NEAR);
        rectmesh = 4 * num - 4;

        acb_mat_init(ests, rectmesh, 3);

        //for very small num it doesn't make sense to split the loop into threads
        //if (num < numthreads*20)
        //   numthreads = 1;
        //if (rectmesh < 3000)
          // numthreads = 1;

        //starting values for the minmodabbs per thread
        arb_mat_init(minmodabb, numthreads, 1);
        for (i = 0; i < numthreads; i++)
        { 
           arb_set_si(arb_mat_entry(minmodabb, i, 0), 1000);
        }

        //prep the threads
        pthread_t thread[numthreads];

        struct ThreadData data[numthreads];

        //prep all the thread data (divide up ranges in the overall loop with start and stop of range)
        slong threadtasks = (num-1)/numthreads;
        for (i = 0; i < numthreads; i++)
        {
            data[i].start=i*threadtasks;
            data[i].stop= (i+1)*threadtasks;
            data[i].id= i;
            data[i].prec= prec;
            data[i].num= num;
            data[i].prt= prt;
        }
        data[numthreads-1].stop = num - 1;

        //start the threads with an indexed (array) of a data-structure (with pointers to it from the threads).
        for (i = 0; i < numthreads; i++)
        {
            pthread_create(&thread[i], 0, abbeff_symmetric_rectangle, &data[i]);
        }

        //wait for all threads to complete
        for (i = 0; i < numthreads; i++)
        {
            pthread_join(thread[i], NULL);
        }

        //calculate and print the winding number for this x,y rectangle
        arb_zero(windnum);
        for (idx = 0; idx < rectmesh-1; idx++)
        { 
            acb_div(argdiv, acb_mat_entry(ests, idx, 0), acb_mat_entry(ests, idx+1, 0), prec);
            acb_arg(a, argdiv, prec);
            arb_add(windnum, windnum, a, prec);
        }   

        acb_div(argdiv, acb_mat_entry(ests, rectmesh-1, 0), acb_mat_entry(ests, 0, 0), prec);
        acb_arg(a, argdiv, prec);
        arb_add(windnum, windnum, a, prec);
        arb_div(windnum, windnum, pi, prec);
        arb_mul_2exp_si(windnum, windnum, -1);

        arb_add(windtot, windtot, windnum, prec);

        //establish minmodabb across all threads
        arb_set_si(minmod, 1000);
        for (i = 0; i < numthreads; i++)
        { 
           if (arb_lt(arb_mat_entry(minmodabb, i, 0), minmod))
               arb_set(minmod, arb_mat_entry(minmodabb, i, 0));
        }

        //print rectangle summary
        printf("Rectangle(%ld) : ", count);
        arf_printd(arb_midref(t), 20);
        printf(", ");
        arf_printd(arb_midref(dtabb), 20);
        printf(", ");
        arf_printd(arb_midref(dzabb), 20);
        printf(", %f, ", arf_get_d(arb_midref(windtot), ARF_RND_NEAR));
        arf_printd(arb_midref(minmod), 20);
        printf(", ");
        printf("%ld\n", rectmesh);

        if (prt==1)
        {      
           for (idx = 0; idx < rectmesh; idx++)
           { 
               acb_get_real(a, acb_mat_entry(ests, idx, 1));
               acb_get_real(b, acb_mat_entry(ests, idx, 2));
               acb_set(ca, acb_mat_entry(ests, idx, 0));
               print_details(prtresult, t, a, b, ca);
           }
        }

        arb_set_si(a, 1);
        if (arb_lt(minmod, a))
        {
            printf("\n");
            printf("!!!   Run has been aborted !!!\n");
            printf("!!! 'modabb' did become < 1 !!!\n");
            printf("\n");
            goto end;
        }

        //establish the next t
        arb_set_si(b, 1);
        arb_mul_2exp_si(a, b, -1);
        arb_sub(a, minmod, a, prec);
        arb_div(a, a, dtabb, prec);
        arb_add(t, t, a, prec); 

        //prevent drift in the error term of t
        arb_get_mid_arb(t,t);

        arb_mat_clear(minmodabb);
        acb_mat_clear(ests);
    }

    flint_printf("\n");
    printf("Overall winding number: %f \n", arf_get_d(arb_midref(windtot), ARF_RND_NEAR));
    printf("\n");

end:

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(x);
    arb_clear(y);
    arb_clear(t);
    arb_clear(Narb);
    arb_clear(xdiv4pi);
    arb_clear(logxdiv4pi);
    arb_clear(onedivxmin6);
    arb_clear(tdiv4divxmin6);
    arb_clear(n0);
    arb_clear(logtn0);
    arb_clear(minmod);
    arb_clear(windnum);
    arb_clear(windtot);
    arb_clear(dzabb);
    arb_clear(dtabb);
    arb_clear(afac);
    arb_clear(const1);
    arb_clear(pi);

    acb_clear(ca);
    acb_clear(argdiv);
    acb_clear(ssexpo);
    acb_clear(one);

    acb_poly_clear(finpoly);
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
    arb_init(X);
    arb_init(y0);
    arb_init(ts);
    arb_init(te);
    arb_init(re);
    arb_init(im);

    acb_t tmp;
    acb_init(tmp);

    const char *ts_str, *te_str, *y0_str, *thread_str, *prt_str;

    slong e, t, N, prec, expterms, taylorterms, prt, res, linesize, numthreads, digits;
    res = 0;
    expterms = 0;
    taylorterms = 0;
    N = 0;
    prec = 0;

    int result = EXIT_SUCCESS;

    linesize=100000;
    char line[linesize];

    if (argc != 7)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    ts_str = argv[1];
    te_str = argv[2];
    y0_str = argv[3];
    prt_str = argv[4];
    thread_str = argv[5];

    prec = 128;

    arb_set_str(ts, ts_str, prec);
    arb_set_str(te, te_str, prec);
    arb_set_str(y0, y0_str, prec);

    prt = atol(prt_str);
    numthreads = atol(thread_str);

    //process finalmata and finalmatb file
    FILE *f = fopen(argv[6], "r");

    //recover X, y0, expterms, taylorterms, digits
    fgets(line, linesize, f);
    char* str = strdup(line);
    arb_set_str(X, getfield(str, 1), prec);
    str = strdup(line);
    expterms = atol(getfield(str, 2));
    str = strdup(line);
    taylorterms = atol(getfield(str, 3));
    str = strdup(line);
    digits = atol(getfield(str, 4));
    free(str);

    acb_mat_t finalmat;
    acb_mat_init(finalmat, expterms, taylorterms);

    //precision
    prec = digits * 3.32192809488736 + 60;

    //fill finalmat
    printf("\n");
    printf("Filling stored sums matrix with %ld X %ld terms guaranteeing %ld digits accuracy...\n", expterms, taylorterms, digits);
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
            acb_set_arb_arb(tmp, re, im);
            acb_set(acb_mat_entry(finalmat, e-1, (t-1)/2), tmp);
            free(str);
            t=t+2;
        }
        e=e+1;
    }

    fclose(f);

    N = get_N(ts, X, prec);

TIMEIT_ONCE_START

    abbeff_t_loop(res, X, y0, ts, te, N, taylorterms, expterms, finalmat, prt, numthreads, digits, prec);

TIMEIT_ONCE_STOP

finish: 
 
    if (result == EXIT_FAILURE)
    {
        flint_printf("Required inputs:\n");
        flint_printf("%s ts te y0 Prt numthreads singlematfile\n\n", argv[0]);
        flint_printf(
    " This script computes the winding number for a '3D-Barrier',\n"
    " that runs along rectangle: [X <= x <= X+1] + i[y0 <= y <= 1],\n"
    " and along: [ts <= t <= te]. It takes X and the required number\n"
    " of Taylor expansion terms from the singlematfile, that contains\n" 
    " a matrix with polynomial coefficients for precalculated sums.\n"
	" With parameter Prt the output can be controlled:\n"
    " 0 = prints rectangle summary only, 1 = prints full details.\n"
    " With numthreads the number of threads csn be selected.\n");
    }

    arb_clear(X);
    arb_clear(y0);
    arb_clear(ts);
    arb_clear(te);
    arb_clear(re);
    arb_clear(im);

    acb_clear(tmp);

    acb_mat_clear(finalmat);
 
    flint_cleanup();

    return result;
}
