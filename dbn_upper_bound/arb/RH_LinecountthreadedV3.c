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
#include "acb_dirichlet.h"
#include "flint/profiler.h"
#include "pthread.h"

//define some global parameters to be accessed from individual threads
    arb_mat_t clientworksplit;
    arb_mat_t threadworksplit;
    arb_mat_t gramblockcounts;
    arb_mat_t Rootsfound;
    arb_mat_t prevgpdata;
    arb_mat_t currgpdata;
    arb_mat_t Zstat;
    arb_mat_t Rosserruleviol;
    arb_mat_t Rosserruletype;

    slong zprec[257];

    struct ThreadData {
    slong id, prec;
    }; 

//count all zeros in the critical strip using the N(t) formula for the contour integral
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

//establish the number of zeros in the strip using the N(t) formula (only at good Gram points)
void
Establish_number_of_roots_in_strip(arb_t res, arb_t start, arb_t end, slong prec)
{
    arb_t a, b;
    arb_init(a);
    arb_init(b);

    N_count(a, start, prec);
    N_count(b, end, prec);
    arb_sub(res, end, start, prec);

    arb_clear(a);
    arb_clear(b);
}

//Establish the required accuracy for Z(t)
slong
z_accuracy(slong zpreccnt, arb_t xa, slong targetprec)
{
    acb_t s, t;
    acb_init(s);
    acb_init(t);

    //prevent drift
    arb_get_mid_arb(xa, xa);
    acb_set_arb(t, xa);

    zpreccnt = targetprec; 
    acb_dirichlet_hardy_z(s, t, NULL, NULL, 1, zpreccnt);
    while(arb_rel_accuracy_bits(acb_realref(s)) < targetprec)
    {
        zpreccnt = zpreccnt + 1;
        acb_dirichlet_hardy_z(s, t, NULL, NULL, 1, zpreccnt);
    }

    acb_clear(s);
    acb_clear(t);

    return(zpreccnt);
}

//Detailed root search to fix violations of Rosser's rule  
slong
CountAllRoots(slong roots, arb_t gvb, arb_t gve, slong depth, slong id, slong prec)
{
    arb_t a, b, c, d, currg, prevg, endz, currz, prevz, currzd, prevzd, stepsize, diffba, accrt, Zdgexact;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(currg);
    arb_init(prevg);
    arb_init(endz);
    arb_init(currz);
    arb_init(prevz);
    arb_init(currzd);
    arb_init(prevzd);
    arb_init(stepsize);
    arb_init(diffba);
    arb_init(accrt);
    arb_init(Zdgexact);

    acb_t s, t;
    acb_init(s);
    acb_init(t);

    acb_struct(z[2]);
    acb_init(z);
    acb_init(z+1);

    slong sectioncnt, zprec1; zprec1 = 0;

    //set required accuracy for evaluating Z'(t) at 10 digits
    zprec1 = z_accuracy(zprec1, gve, 10 * 3.32192809488736 + 5);

    arb_set_si(a, 1000000000);
    arb_inv(accrt, a, prec);

    //divide the difference between the end and start Gram point values by the number of search sections chosen 
    arb_sub(a, gve, gvb, prec);
    arb_set_si(b, depth);
    arb_div(stepsize, a, b, prec);

    arb_set(prevg, gvb);
    acb_set_arb(t, gvb);
    acb_dirichlet_hardy_z(z, t, NULL, NULL, 2, zprec1);
    acb_get_real(prevz, z);
    acb_get_real(prevzd, z + 1);
    arb_add_si(arb_mat_entry(Zstat, id, 2), arb_mat_entry(Zstat, id, 2), 1, prec);

    acb_set_arb(t, gve);
    acb_dirichlet_hardy_z(z, t, NULL, NULL, 1, zprec1);
    acb_get_real(endz, z);
    arb_add_si(arb_mat_entry(Zstat, id, 2), arb_mat_entry(Zstat, id, 2), 1, prec);

    //loop through all sections of the interval gvb to gve
    sectioncnt = 0; roots = 0;
    while((sectioncnt < depth))
    {
        arb_mul_si(a, stepsize, sectioncnt + 1, prec);
        arb_add(currg, gvb, a, prec);

        acb_set_arb(t, currg);
        acb_dirichlet_hardy_z(z, t, NULL, NULL, 2, zprec1);
        acb_get_real(currzd, z + 1);
        arb_add_si(arb_mat_entry(Zstat, id, 2), arb_mat_entry(Zstat, id, 2), 1, prec);

        //check for a sign change of Z'(t)
        arb_mul(a, currzd, prevzd, prec);
        if (arb_is_negative(a))
        {
           //if sign change of Z'(t), then find the exact (at 'accrt' precision) value where Z(t) = 0.
           arb_set(a, prevg);
           arb_set(b, currg);
           arb_set(Zdgexact, a);
           arb_sub(diffba, b, a, prec);

           while(arb_ge(diffba, accrt))
           {
                arb_add(c, a, b , prec);
                arb_mul_2exp_si(Zdgexact, c, -1);
                acb_set_arb(t, Zdgexact);
                acb_dirichlet_hardy_z(z, t, NULL, NULL, 2, zprec1);
                acb_get_real(c, z + 1);
                acb_set_arb(t, a);
                acb_dirichlet_hardy_z(z, t, NULL, NULL, 2, zprec1);
                acb_get_real(d, z + 1);
                arb_mul(d, c, d, prec);
                arb_add_si(arb_mat_entry(Zstat, id, 2), arb_mat_entry(Zstat, id, 2), 2, prec);
                if (arb_is_negative(d))
                   arb_set(b, Zdgexact);
                else
                   arb_set(a, Zdgexact);

                arb_sub(diffba, b, a, prec);
           }  

           //check for sign change between Z(t) at the current Z'(t) = 0 and the previous Z'(t) = 0.
           acb_get_real(currz, z);
           arb_mul(a, currz, prevz, prec);
           if (arb_is_negative(a))
              roots++;

           arb_set(prevz, currz);
        }

        arb_set(prevg, currg);
        arb_set(prevzd, currzd);
        sectioncnt++;
    }

    //count an additional root when Z(previous)*Z(end) < 0
    arb_mul(a, prevz, endz, prec);
    if (arb_is_negative(a))
    roots++;

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(currg);
    arb_clear(prevg);
    arb_clear(endz);
    arb_clear(currz);
    arb_clear(prevz);
    arb_clear(currzd);
    arb_clear(prevzd);
    arb_clear(stepsize);
    arb_clear(diffba);
    arb_clear(accrt);
    arb_clear(Zdgexact);

    acb_clear(s);
    acb_clear(t);

    acb_clear(z);
    acb_clear(z+1);

    return(roots);
}

//Find all roots (sign changes) of Z(t) on the critical line  
slong
CountRootsIntervalZ(slong introots, arb_mat_t intsections, slong intnum, arb_t gvb, arb_t gve, arb_t gzb, arb_t gze,
                    slong numsections, slong target, slong maxdepthpow, slong id, slong prec)
{
    arb_t a, b, currz, prevz, stepsize;
    arb_init(a);
    arb_init(b);
    arb_init(currz);
    arb_init(prevz);
    arb_init(stepsize);

    acb_t s, t;
    acb_init(s);
    acb_init(t);

    slong sectioncnt, sectionendp;

    arb_set(prevz, gzb);
    arb_set(arb_mat_entry(intsections, intnum, 0), gzb);
    arb_set(arb_mat_entry(intsections, intnum, maxdepthpow), gze);

    //divide the difference between the end and start Gram point values by the number of search sections chosen 
    arb_sub(a, gve, gvb, prec);
    arb_set_si(b, numsections);
    arb_div(stepsize, a, b, prec);

    //loop through all search sections and stop once the target number of roots of Z(t) has been found. 
    sectioncnt = 0; sectionendp = maxdepthpow/numsections; introots = 0;
    while((sectioncnt < numsections) && (introots < target))
    {

        //only evaluate Z(t) when it concerns a new calculations (older values are stored in the matrix 'intsections')
        if(arb_is_zero(arb_mat_entry(intsections, intnum, sectionendp)))
        {
            arb_mul_si(a, stepsize, sectioncnt + 1, prec);
            arb_add(a, a, gvb, prec);
            acb_set_arb(t, a);
            acb_dirichlet_hardy_z(s, t, NULL, NULL, 1, zprec[id]);
            acb_get_real(currz, s);
            arb_set(arb_mat_entry(intsections, intnum, sectionendp), currz);
            arb_add_si(arb_mat_entry(Zstat, id, 1), arb_mat_entry(Zstat, id, 1), 1, prec);
        }
        else
            arb_set(currz, arb_mat_entry(intsections, intnum, sectionendp)); 

        //check for a sign change in this subsection of the interval
        arb_mul(a, currz, prevz, prec);
        if (arb_is_negative(a))
           introots++;

        arb_set(prevz, currz);

        sectioncnt++; sectionendp = (sectioncnt + 1)*maxdepthpow/numsections;     
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(currz);
    arb_clear(prevz);
    arb_clear(stepsize);

    acb_clear(s);
    acb_clear(t);

    return(introots);
}

//Establish the optimal searching approach for this Gramblock
slong
EstablishSearchPattern(slong rootsfound, arb_mat_t gramblock, slong blocklength, slong id, slong prec)
{
    arb_t gvb, gve, gzb, gze, gbsz, gbez;
    arb_init(gvb);
    arb_init(gve);
    arb_init(gzb);
    arb_init(gze);
    arb_init(gbsz);
    arb_init(gbez);

    slong intcnt, numsections, target, maxdepth;
    rootsfound = 0; numsections = 1; maxdepth = 256;

    //recover Z(t)-values at the beginning and end of the Gram block.
    arb_abs(gbsz, arb_mat_entry(gramblock, 0, 3)); arb_abs(gbez, arb_mat_entry(gramblock, blocklength, 3));

    //all calculated values of Z(t) for each subsection of an interval will be stored here.
    arb_mat_t(intsections);
    arb_mat_init(intsections, blocklength + 1, maxdepth);

    //run through all intervals of this Gram block with an increasing subdivision into numsections.
    while (numsections <= maxdepth && rootsfound < blocklength)
    {
       numsections = numsections * 2;

       //search Gram block from left to right if right outer Z-value is larger than left outer Z-value
       if (arb_le(gbsz, gbez))
       {
          //check first outer left interval (target 1 zero) then the middle one (target 2) and then the outer right (target 1)
          intcnt = 0;
          while(intcnt < blocklength)
          {
              target = 2;
              if (intcnt == 0 || intcnt == blocklength - 1)
                 target = 1;
              arb_set(gvb, arb_mat_entry(gramblock, intcnt, 2)); arb_set(gve, arb_mat_entry(gramblock, intcnt + 1, 2)); 
              arb_set(gzb, arb_mat_entry(gramblock, intcnt, 3)); arb_set(gze, arb_mat_entry(gramblock, intcnt + 1, 3));  
              rootsfound = CountRootsIntervalZ(rootsfound, intsections, intcnt, gvb, gve, gzb, gze, numsections, target, maxdepth, id, prec);

              if (rootsfound != target)
                 intcnt++;
              else
              {
                 rootsfound = blocklength; 
                 break;
              }
          } 
       }

       //search Gram block from right to left
       else 
       {
          //check first outer right interval (target 1 zero) then the preceding middle ones (target 2) and then the left outer (target 1)
          intcnt = blocklength - 1;
          while(intcnt >= 0)
          {
              target = 2;
              if (intcnt == 0 || intcnt == blocklength - 1)
                 target = 1;
              arb_set(gvb, arb_mat_entry(gramblock, intcnt, 2)); arb_set(gve, arb_mat_entry(gramblock, intcnt + 1, 2)); 
              arb_set(gzb, arb_mat_entry(gramblock, intcnt, 3)); arb_set(gze, arb_mat_entry(gramblock, intcnt + 1, 3));  
   
              rootsfound = CountRootsIntervalZ(rootsfound, intsections, intcnt, gvb, gve, gzb, gze, numsections, target, maxdepth, id, prec);

              if (rootsfound != target)
                 intcnt--;  
              else
              {
                 rootsfound = blocklength; 
                 break;
              }
          } 
       }
    }

    arb_clear(gvb);
    arb_clear(gve);
    arb_clear(gzb);
    arb_clear(gze);
    arb_clear(gbsz);
    arb_clear(gbez);

    arb_mat_clear(intsections);

    return(rootsfound);
}

//test procedure for gram function abs(tn-tn1) > gram target precision. 0 = false, 1 = true.
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

//Establish good approximation of the n-th gram point (gn)
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
    arb_set_si(b, 10000000000);
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

//Test procedure to establish whether a gram point is good or bad. 0 = bad, 1 = good.
int
GoodGram(slong gramquality, arb_t gi, slong id, slong prec)
{
    arb_t a, b, gq, gv, gz;
    arb_init(a);
    arb_init(b);
    arb_init(gq);
    arb_init(gv);
    arb_init(gz);

    acb_t s, t;
    acb_init(s);
    acb_init(t);

    slong i;

    gramquality = 0;

    //Establish Gram value and prevent drift
    gram(gv, gi, prec);
    arb_get_mid_arb(gv, gv);

    //Establish the value of Z(g_n)
    acb_set_arb(t, gv);
    acb_dirichlet_hardy_z(s, t, NULL, NULL, 1, zprec[id]);
    acb_get_real(gz, s);
    arb_add_si(arb_mat_entry(Zstat, id, 0), arb_mat_entry(Zstat, id, 0), 1, prec);

    //In case of doubt, when |Z(gp_m)| < 10e-8 then just assume the Gram block is Bad
    arb_set_si(a, 100000000); arb_inv(a, a, prec); arb_abs(b, gz);
    if (arb_lt(b, a))
    {
        gramquality = 0;
    }
    else
    {
        //Establish (-1)^gi
        arb_set_si(b, -1); arb_pow(b, b, gi, prec);

        //Determine quality (-1)^m*(Z(gp_m)) > 0 = Good (i.e. gramquality = 1)
        arb_mul(a, b, gz, prec);
        arb_zero(b);
        if (arb_gt(a, b))
           gramquality = 1;
    }

    //keep all Gram point data from current and previous points
    for (i=0; i < 4; i++)
        arb_set(arb_mat_entry(prevgpdata, id, i), arb_mat_entry(currgpdata, id, i));

    arb_set_si(gq, gramquality);
    arb_set(arb_mat_entry(currgpdata, id, 0), gi);
    arb_set(arb_mat_entry(currgpdata, id, 1), gq);
    arb_set(arb_mat_entry(currgpdata, id, 2), gv);
    arb_set(arb_mat_entry(currgpdata, id, 3), gz);

    arb_clear(a);
    arb_clear(b);
    arb_clear(gq);
    arb_clear(gv);
    arb_clear(gz);

    acb_clear(s);
    acb_clear(t);

    return(gramquality);
}

//Determine the (middle) index of a set of three-in-a-row Good Gram points starting from searchvalue.
void
Find_prev_three_good_grams(arb_t result, arb_t searchvalue, slong prec)
{
    arb_t a, b, c, d;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);

    slong gramquality, targetprec;
    gramquality = 0; targetprec = 0;

    arb_add_si(d, searchvalue, -1, prec);

    //determine the required decimal precision (3 digits) for the search value Z(gram(gi))
    targetprec = 3 * 3.32192809488736 + 10;
    gram(a, d, prec);
    zprec[256] = z_accuracy(zprec[256], a, targetprec);

    //loop until a consecutive pattern of G G G Grampoints has been found and pick the middle one.
    arb_set(a, d);
    arb_add_si(b, a, 1, prec);
    arb_add_si(c, b, 1, prec);
    while ((GoodGram(gramquality, a, 256, prec) == 0) || (GoodGram(gramquality, b, 256, prec) == 0)
	        || (GoodGram(gramquality, c, 256, prec) == 0))
    {
        arb_add_si(a, a, -1, prec);
        arb_add_si(b, b, -1, prec);
        arb_add_si(c, c, -1, prec);
    }

    arb_set(result, b);

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
}

//Determine the (middle) index of a set of three-in-a-row Good Gram points starting from searchvalue.
void
Find_next_three_good_grams(arb_t result, arb_t searchvalue, slong prec)
{
    arb_t a, b, c, d;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);

    slong gramquality, targetprec;
    gramquality = 0;

    arb_add_si(d, searchvalue, -1, prec);

    //determine the required decimal precision (3 digits) for the search value Z(gram(gi))
    targetprec = 3 * 3.32192809488736 + 10;
    gram(a, d, prec);
    zprec[256] = z_accuracy(zprec[256], a, targetprec);

    //loop until a consecutive pattern of G G G Grampoints has been found and pick the middle one.
    arb_set(a, d);
    arb_add_si(b, a, 1, prec);
    arb_add_si(c, b, 1, prec);
    while ((GoodGram(gramquality, a, 256, prec) == 0) || (GoodGram(gramquality, b, 256, prec) == 0)
	        || (GoodGram(gramquality, c, 256, prec) == 0))
    {
        arb_add_si(a, a, 1, prec);
        arb_add_si(b, b, 1, prec);
        arb_add_si(c, c, 1, prec);
    }

    arb_set(result, b);

    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
}

//determine the chunks of work for each Thread. Ensure each 'chunk' starts and ends as the midpoint of G G G.
void
Allocate_work_to_Thread(arb_t rcount, const arb_t start, const arb_t end, const slong numthreads, slong prec)

{  
    arb_t a, gpcnt, numthreadarb, threadsize, blockstart, blockend;
    arb_init(a);
    arb_init(gpcnt);
    arb_init(numthreadarb);
    arb_init(threadsize);
    arb_init(blockstart);
    arb_init(blockend);

    slong n;

    //determine worksize per thread
    arb_set_si(numthreadarb, numthreads);
    arb_sub(a, end, start, prec);
    arb_div(threadsize, a, numthreadarb, prec);
    arb_floor(threadsize, threadsize, prec);

    //loop through all threads and establish beginning and end points (minding G G G).
    n = 0;
    arb_set(blockstart, start);
    arb_add(blockend, start, threadsize, prec);
    while (n < numthreads)
    {
       Find_next_three_good_grams(gpcnt, blockend, prec);

       arb_set(arb_mat_entry(threadworksplit, n, 0), blockstart);
       arb_set(arb_mat_entry(threadworksplit, n, 1), gpcnt);

       arb_add_si(blockstart, gpcnt, 0, prec);
       arb_add(blockend, gpcnt, threadsize, prec);
       n = n + 1;
    }

    arb_set(arb_mat_entry(threadworksplit, n - 1, 1), end);

    arb_clear(a);
    arb_clear(gpcnt);
    arb_clear(numthreadarb);
    arb_clear(threadsize);
    arb_clear(blockstart);
    arb_clear(blockend);
}

//Determine the chunks of work per client and establish start and end values (ensuring midpoint of G G G) for clientid.
void
Allocate_work_to_Client(arb_t rcount, const arb_t start, const arb_t end, const slong numclients, slong prec)
{  
    arb_t a, gpcnt, numclientarb, clientsize, blockstart, blockend;
    arb_init(a);
    arb_init(gpcnt);
    arb_init(numclientarb);
    arb_init(clientsize);
    arb_init(blockstart);
    arb_init(blockend);

    slong n;

    //determine worksize per client
    arb_set_si(numclientarb, numclients);
    arb_sub(a, end, start, prec);
    arb_div(clientsize, a, numclientarb, prec);
    arb_floor(clientsize, clientsize, prec);

    //loop through all client work until clientid and establish beginning and end points (minding G G G).
    n = 0;
    arb_set(blockstart, start);
    arb_add(blockend, start, clientsize, prec);
    while (n < numclients)
    {
       Find_next_three_good_grams(gpcnt, blockend, prec);

       arb_set(arb_mat_entry(clientworksplit, n, 0), blockstart);
       arb_set(arb_mat_entry(clientworksplit, n, 1), gpcnt);

       arb_add_si(blockstart, gpcnt, 0, prec);
       arb_add(blockend, gpcnt, clientsize, prec);
       n = n + 1;
    }

    arb_set(arb_mat_entry(clientworksplit, n - 1, 1), end);

    arb_clear(a);
    arb_clear(gpcnt);
    arb_clear(numclientarb);
    arb_clear(clientsize);
    arb_clear(blockstart);
    arb_clear(blockend);
}

//Thread-calculations to assess whether Gram points are 'good' or 'bad' and fix the bad ones.
void* isolate_bad_Gram_points(void *voidData)
{
    //Recover the data passed to this specific thread
    struct ThreadData* data=voidData;
    slong id=data->id;
    slong prec=data->prec;

    arb_t a, gi, gf;
    arb_init(a);
    arb_init(gi);
    arb_init(gf);

    //Data of previous and current Gramblocks (or G-G intervals) are kept here
    arb_mat_t gramblock, prevgramblock;
    arb_mat_init(gramblock, 30, 4);
    arb_mat_init(prevgramblock, 30, 4);

    slong i, j, blocklength, prevblklen, targetprec, gramquality, prevqual, roots, missingzeros, surpluszeros, rossercnt;

    rossercnt = 0;

    //Get the allocated start and end values of the range for this thread
    arb_set(gi, arb_mat_entry(threadworksplit, id, 0));
    arb_set(gf, arb_mat_entry(threadworksplit, id, 1));

    gramquality = 0; prevqual = 1; blocklength = 0; missingzeros = 0; surpluszeros = 0;

    //Determine the required decimal precision (3 digits) for the start value Z(gram(gi))
    targetprec = 3 * 3.32192809488736 + 5;
    gram(a, gf, prec);
    zprec[id] = z_accuracy(prec, a, targetprec);

    while (arb_lt(gi, gf))
    {
       arb_add_si(gi, gi, 1, prec);

       //Establish Gram point quality
       gramquality = GoodGram(gramquality, gi, id, prec);

       //Good Gram point followed by Good (a good 'Gramblock').
       if(prevqual && gramquality)
       {    
           arb_add_si(arb_mat_entry(gramblockcounts, id, 0), arb_mat_entry(gramblockcounts, id, 0), 1, prec);
           arb_add_si(arb_mat_entry(Rootsfound, id, 0), arb_mat_entry(Rootsfound, id, 0), 1, prec);

           //Type I Rosser's rule violation in previous block? Check if any missing zeros can be found now. 
           if(missingzeros > 0)
           {
               roots = CountAllRoots(roots, arb_mat_entry(prevgpdata, id, 2), arb_mat_entry(currgpdata, id, 2), 16, id, prec);
               if(missingzeros + 1 - roots <= 0)
               {
                   arb_add_si(arb_mat_entry(Rootsfound, id, 2), arb_mat_entry(Rootsfound, id, 2), prevblklen, prec);
                   rossercnt++;
                   arb_set(arb_mat_entry(Rosserruleviol, id, rossercnt), arb_mat_entry(prevgramblock, 0, 0));
                   arb_set_si(arb_mat_entry(Rosserruletype, id, rossercnt), 1);
                   surpluszeros = roots - missingzeros - 1;
                   missingzeros = 0;
               }
           }
           //Save this interval as a 'gram block' for reference in case a Rosser type II violation occurs
           for (i = 0; i < 4; i++)
           {
               arb_set(arb_mat_entry(prevgramblock, 0, i), arb_mat_entry(prevgpdata, id, i));
               arb_set(arb_mat_entry(prevgramblock, 1, i), arb_mat_entry(currgpdata, id, i));
           }
           prevblklen = 1;
       }

       //Good Gram point followed by Bad (start of a new Gram block)
       if(prevqual && gramquality == 0)
       {
           blocklength = 1;
           for (i = 0; i < 4; i++)
           {
               arb_set(arb_mat_entry(gramblock, 0, i), arb_mat_entry(prevgpdata, id, i));
               arb_set(arb_mat_entry(gramblock, blocklength, i), arb_mat_entry(currgpdata, id, i));
           }
       }

       //Bad Gram point followed by Bad (inner part of a Gram Block)
       if(prevqual== 0 && gramquality == 0)
       {
           blocklength = blocklength + 1;
           for (i = 0; i < 4; i++)
           {
               arb_set(arb_mat_entry(gramblock, blocklength, i), arb_mat_entry(currgpdata, id, i));
           }
       }

       //Bad Gram point followed by Good (end of a Gram block)
       if(prevqual== 0 && gramquality)
       {
           arb_add_si(arb_mat_entry(gramblockcounts, id, blocklength), arb_mat_entry(gramblockcounts, id, blocklength), 1, prec);
           blocklength = blocklength + 1;
           for (i = 0; i < 4; i++)
           {
               arb_set(arb_mat_entry(gramblock, blocklength, i), arb_mat_entry(currgpdata, id, i));
           }

           //Type III Rosser's rule violation in previous block? Check wether any missing zeros can be found now. 
           if(missingzeros > 0)
           {
               roots = CountAllRoots(roots, arb_mat_entry(gramblock, 0, 2), arb_mat_entry(gramblock, blocklength, 2), blocklength*16, id, prec);
               if(missingzeros + blocklength - roots <= 0)
               {
                   arb_add_si(arb_mat_entry(Rootsfound, id, 4), arb_mat_entry(Rootsfound, id, 4), prevblklen, prec);
                   rossercnt++;
                   arb_set(arb_mat_entry(Rosserruleviol, id, rossercnt), arb_mat_entry(prevgramblock, 0, 0));
                   arb_set_si(arb_mat_entry(Rosserruletype, id, rossercnt), 3);
                   surpluszeros = roots - missingzeros - blocklength;
                   missingzeros = 0;
               }
               else
               {
                   printf("!!! Failed to find the missing zeros of Rosser's rule violation of Gram point :");
                   arb_printd(arb_mat_entry(prevgramblock, 0, 0), 20);
                   printf(" !!! \n Please inspect manually! \n");
                   missingzeros = 0;
               }
           }

           //Check the current bad Gram block or establish any missing zeros (indicating a violation of Rosser's rule)
          roots = EstablishSearchPattern(missingzeros, gramblock, blocklength, id, prec);
          
           if(roots == blocklength)
                 arb_add_si(arb_mat_entry(Rootsfound, id, 1), arb_mat_entry(Rootsfound, id, 1), blocklength, prec);
           else
           //If not all zeros can be found perform a rigorous search in the current Gramblock to establish all its roots
           {
                roots = CountAllRoots(roots, arb_mat_entry(gramblock, 0, 2), arb_mat_entry(gramblock, blocklength, 2), blocklength*16, id, prec);
                missingzeros = blocklength - roots;

                if(missingzeros == 0)
                      arb_add_si(arb_mat_entry(Rootsfound, id, 1), arb_mat_entry(Rootsfound, id, 1), blocklength, prec);
                else
                //Type II Rosser's rule violation in current block? Check if missing zeros can be found in the preceding block or equal an outstanding surplus.
                { 
                   if(missingzeros == surpluszeros)
                   {
                       arb_add_si(arb_mat_entry(Rootsfound, id, 3), arb_mat_entry(Rootsfound, id, 3), missingzeros, prec);
                       rossercnt++;
                       arb_set(arb_mat_entry(Rosserruleviol, id, rossercnt), arb_mat_entry(gramblock, 0, 0));
                       arb_set_si(arb_mat_entry(Rosserruletype, id, rossercnt), 2);
                       surpluszeros = 0;
                       missingzeros = 0;
                   }
                   else
                   {
                       roots = CountAllRoots(roots, arb_mat_entry(prevgramblock, 0, 2), arb_mat_entry(prevgramblock, prevblklen, 2), prevblklen*16, id, prec);
                       if(missingzeros + prevblklen - roots == 0)
                       {
                            arb_add_si(arb_mat_entry(Rootsfound, id, 3), arb_mat_entry(Rootsfound, id, 3), blocklength, prec);
                            rossercnt++;
                            arb_set(arb_mat_entry(Rosserruleviol, id, rossercnt), arb_mat_entry(gramblock, 0, 0));
                            arb_set_si(arb_mat_entry(Rosserruletype, id, rossercnt), 2);
                            missingzeros = 0;
                       }
                    }
                }
           }

           //Save the current Gram block data for reference when a Type II Rosser rule violation occurs
           for (j = 0; j <= blocklength; j++)
           {
              for (i = 0; i < 4; i++)
              {
                  arb_set(arb_mat_entry(prevgramblock, j, i), arb_mat_entry(gramblock, j, i));
              }
           }
           prevblklen = blocklength;
       }

       prevqual = gramquality;
    }

    arb_clear(a);
    arb_clear(gi);
    arb_clear(gf);

    arb_mat_clear(gramblock);
    arb_mat_clear(prevgramblock);

    return(NULL);
}

//main program that issues and rejoins all threads
void
Count_number_of_roots_on_line(arb_t totcount, slong numthreads, slong prec)
{
    arf_t u;
    arf_init(u);

    arb_t subcount;
    arb_init(subcount);

    slong i, j, tmp, tmp1;

    double d1, d2;
    
    //prep the threads
    pthread_t thread[numthreads];

    struct ThreadData data[numthreads];

    //prep the thread data
    for (i = 0; i < numthreads; i++)
    {
        data[i].id=i;
        data[i].prec=prec;
    }

    //start the various threads.
    for (i = 0; i < numthreads; i++)
    {
        pthread_create(&thread[i], 0, isolate_bad_Gram_points, &data[i]); 
    }

    //wait for all threads to complete
    for (i = 0; i < numthreads; i++)
    {
        pthread_join(thread[i], NULL);
    }

    //print all Rosser's rule violations
    printf("\nSummary of Rosser rule violations:\n");
    for (j = 0; j < 10000; j++)
    {    
       for (i = 0; i < numthreads; i++)
       { 
          arb_get_ubound_arf(u, arb_mat_entry(Rosserruleviol, i, j), prec);
          tmp=arf_get_si(u, ARF_RND_NEAR);   
          arb_get_ubound_arf(u, arb_mat_entry(Rosserruletype, i, j), prec);
          tmp1=arf_get_si(u, ARF_RND_NEAR);   
          if (tmp > 0)
          {
              printf("%ld of type: %ld\n", tmp, tmp1);
          }
       }
    }
    printf("\n\n");

    //print the number of Gram points per Gramblock length
    arb_zero(totcount);
    for (j = 0; j < 30; j++)
    {    
       arb_zero(subcount);
       for (i = 0; i < numthreads; i++)
       { 
           arb_add(subcount, subcount, arb_mat_entry(gramblockcounts, i, j), prec);
       }
       if (arb_is_positive(subcount))
       {
           arb_get_ubound_arf(u, subcount, prec);
           tmp=arf_get_si(u, ARF_RND_NEAR);
           printf("Total gramblocks of length %ld : %ld \n", j+1, tmp);
           arb_add(totcount, totcount, subcount, prec);
       }
    }
    arb_get_ubound_arf(u, totcount, prec);
    tmp=arf_get_si(u, ARF_RND_NEAR);
    printf("Total gramblocks  ========== : %ld \n", tmp);

    //print the statistics on evaluating the Z(t) function
	printf("\n");
    arb_zero(totcount);
    for (j = 0; j < 3; j++)
    {    
       arb_zero(subcount);
       for (i = 0; i < numthreads; i++)
       { 
           arb_add(subcount, subcount, arb_mat_entry(Zstat, i, j), prec);
       }
       if (arb_is_positive(subcount))
       {
           arb_get_ubound_arf(u, subcount, prec);
           tmp=arf_get_si(u, ARF_RND_NEAR);
           printf("Evaluations of Z(t) required ");
           switch (j)
           {
              case 0:
                 printf("to assess a Gram point   : %ld \n", tmp); tmp1 = tmp; break;
              case 1:
                 printf("to fix bad Gram points   : %ld \n", tmp); break;
              case 2:
                 printf("to fix Rosser violations : %ld \n", tmp); break;
              default:
                 break;
           }
           arb_add(totcount, totcount, subcount, prec);
       }
    }
    arb_get_ubound_arf(u, totcount, prec);
    tmp=arf_get_si(u, ARF_RND_NEAR);
    printf("Total evaluations of Z(t)  ========================== : %ld \n\n", tmp);

    d1 = tmp;
    d2 = tmp1;
    printf("Average number of evaluations of Z(t) required for this range = %3.3f \n\n", d1/d2);

    //Print the roots found per category
    arb_zero(totcount);
    for (j = 0; j < 5; j++)
    {    
       arb_zero(subcount);
       for (i = 0; i < numthreads; i++)
       { 
           arb_add(subcount, subcount, arb_mat_entry(Rootsfound, i, j), prec);
       }
       if (arb_is_positive(subcount))
       {
           arb_get_ubound_arf(u, subcount, prec);
           tmp=arf_get_si(u, ARF_RND_NEAR);
           printf("Roots of Z(t) found in ");
           switch (j)
           {
              case 0:
                 printf("Good to Good intervals      : %ld \n", tmp); break;
              case 1:
                 printf("Gramblocks G B...B G        : %ld \n", tmp); break;
              case 2:
                 printf("current & next G-G interval : %ld (Rosser violation Type 1)\n", tmp); break;
              case 3:
                 printf("current & previous interval : %ld (Rosser violation Type 2)\n", tmp); break;
              case 4:
                 printf("current & next Gramblock    : %ld (Rosser violation Type 3)\n", tmp); break;
              default:
                 break;
           }
           arb_add(totcount, totcount, subcount, prec);
       }
    }
    arb_get_ubound_arf(u, totcount, prec);
    tmp=arf_get_si(u, ARF_RND_NEAR);
    printf("Total roots of Z(t)  ============================= : %ld \n\n", tmp);

    arb_clear(subcount);

    arf_clear(u);
}

int main(int argc, char *argv[])
{
    arb_t xa, xb, res, rootsinstrip, rootsonline;
    arb_init(xa);
    arb_init(xb);
    arb_init(res);
    arb_init(rootsinstrip);
    arb_init(rootsonline);

    arf_t u;
    arf_init(u);

    const char *xa_str, *xb_str, *numclients_str, *clientid_str, *numthread_str;
    int result = EXIT_SUCCESS;
    slong numclients, clientid, numthreads, prec, tmp, tmp1;

    if (argc != 6)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    //Initial precision 
    prec = 90;

    xa_str = argv[1];
    arb_set_str(xa, xa_str, prec);
	
    xb_str = argv[2];
    arb_set_str(xb, xb_str, prec);

    if(arb_le(xb, xa))
    {
        result = EXIT_FAILURE; 
        goto finish;
    }

    numclients_str = argv[3];
    numclients = atol(numclients_str);

    clientid_str = argv[4];
    clientid = atol(clientid_str);

    numthread_str = argv[5];
    numthreads = atol(numthread_str);

    if(clientid < 0 || clientid > numclients - 1 || numclients < 1 || numthreads < 1 || numthreads > 254)
    {
        result = EXIT_FAILURE; 
        goto finish;
    }

    arb_mat_init(clientworksplit, numclients, 2);
    arb_mat_init(threadworksplit, numthreads, 2);
    arb_mat_init(gramblockcounts, numthreads, 30);
    arb_mat_init(Rootsfound, numthreads, 5);
    arb_mat_init(prevgpdata, 257, 4);
    arb_mat_init(currgpdata, 257, 4);
    arb_mat_init(Zstat, 257, 3);
    arb_mat_init(Rosserruleviol, numthreads, 10000);
    arb_mat_init(Rosserruletype, numthreads, 10000);

//tmp=0;
//tmp = CountAllRoots(tmp, xa, xb, 16, 256, prec);
//printf("Roots: %ld \n", tmp);
//goto finish;

    arb_set_si(res, 4);
    if(arb_gt(xa, res))
    {
       Find_prev_three_good_grams(res, xa, prec);
       if (arb_lt(res, xa))
       {
          arb_get_ubound_arf(u, res, prec);
          slong tmp=arf_get_si(u, ARF_RND_NEAR);
          printf("\nWarning: the starting point you have selected is not the mid point of three Good Gram points,\n");
          printf("         it has therefore been automatically decreased to: %ld. \n\n", tmp);
          arb_set(xa, res);
       }
    }

    Find_next_three_good_grams(res, xb, prec);
    if (arb_gt(res, xb))
    {
       arb_get_ubound_arf(u, res, prec);
       slong tmp=arf_get_si(u, ARF_RND_NEAR);
       printf("\nWarning: the end point you have selected is not the mid point of three Good Gram points,\n");
       printf("         it has therefore been automatically augmented to: %ld. \n\n", tmp);
       arb_set(xb, res);
    }

    Allocate_work_to_Client(res, xa, xb, numclients, prec);
    if(arb_ge(arb_mat_entry(clientworksplit, numclients-1, 0), arb_mat_entry(clientworksplit, numclients-1, 1)))
    {
         printf("! The chosen number of clients does not allow to adequately split the selected range !\n");
         printf("Either reduce the number of clients or widen the range of Gram points to be verified.\n\n");
         goto finish;
    }

    Allocate_work_to_Thread(res, arb_mat_entry(clientworksplit, clientid, 0), arb_mat_entry(clientworksplit, clientid, 1), numthreads, prec);
    if(arb_ge(arb_mat_entry(threadworksplit, numthreads-1, 0), arb_mat_entry(threadworksplit, numthreads-1, 1)))
    {
         printf("! The chosen number of threads does not allow to adequately split the selected range !\n");
         printf("Either reduce the number of threads or widen the range of Gram points to be verified.\n\n");
         goto finish;
    }

    Establish_number_of_roots_in_strip(rootsinstrip, arb_mat_entry(clientworksplit, clientid, 0), arb_mat_entry(clientworksplit, clientid, 1), prec);

TIMEIT_ONCE_START

    Count_number_of_roots_on_line(rootsonline, numthreads, prec);

    arb_get_ubound_arf(u, arb_mat_entry(clientworksplit, clientid, 0), prec);
    tmp=arf_get_si(u, ARF_RND_NEAR);
    arb_get_ubound_arf(u, arb_mat_entry(clientworksplit, clientid, 1), prec);
    tmp1=arf_get_si(u, ARF_RND_NEAR);
    printf("In the range of Gram points from %ld to %ld, \n", tmp, tmp1);
    arb_get_ubound_arf(u, rootsonline, prec);
    tmp=arf_get_si(u, ARF_RND_NEAR);
    arb_get_ubound_arf(u, rootsinstrip, prec);
    tmp1=arf_get_si(u, ARF_RND_NEAR);
    printf("%ld roots of Z(t) have been found on the line and %ld in the strip.\n\n", tmp, tmp1);

    if (arb_eq(rootsinstrip, rootsonline))
       printf("SUCCESS");
    else
       printf("FAILURE");
    printf("\n\n");

TIMEIT_ONCE_STOP

finish:
    if (result == EXIT_FAILURE)
    {
        flint_printf("Required inputs:\n");
        flint_printf("%s Gp_b Gp_e numclients clientid numthreads \n\n", argv[0]);
        flint_printf(
    "This script aims to count all zeros of Z(t) on the critical line \n"
    "in the range: Grampoint Gp_b to Grampoint Gp-e. The work can be \n"
    "parallellized through grid computing by chosing a desired number of\n" 
    "clients (1..N) and a clientid (0..N-1) for this run. When numclients = 1 \n"
    "and clientid = 0 the script is processed in stand-alone mode.\n"
    "Further parallellization per client can be done through threading (1..255).\n"
    "The accuracy of Z(t) is automatically kept at >= 3 digits.\n");
    }

    arb_mat_clear(clientworksplit);
    arb_mat_clear(gramblockcounts);
    arb_mat_clear(threadworksplit);
    arb_mat_clear(Rootsfound);
    arb_mat_clear(prevgpdata);
    arb_mat_clear(currgpdata);
    arb_mat_clear(Zstat);
    arb_mat_clear(Rosserruleviol);
    arb_mat_clear(Rosserruletype);

    arb_clear(xa);
    arb_clear(xb);
    arb_clear(res);
    arb_clear(rootsinstrip);
    arb_clear(rootsonline);

    arf_clear(u);
 
    flint_cleanup();

    return result;
}
