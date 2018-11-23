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

//define some global parameters to be accessed by individual threads
    arb_mat_t clientworksplit;
    arb_mat_t threadworksplit;
    arb_mat_t gramblockcounts;
    arb_mat_t Rootsfound;
    arb_mat_t currgi;
    arb_mat_t currgv;
    arb_mat_t currgz;
    arb_mat_t Zstat;

    slong zprec[257];

    struct ThreadData {
    slong id, prec;
    }; 

//count all zeros in the critical strip
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

//Find all roots (sign changes) of Z(t) on the critical line  (depth=1 means 2 intervals, etc.)
int
CountRootZ(slong roots, arb_t begp, arb_t endp, arb_t begpz, arb_t endpz, arb_t depth, slong id, slong prec)
{
    arb_t a, b, w, currz, prevz, stepsize;
    arb_init(a);
    arb_init(b);
    arb_init(w);
    arb_init(currz);
    arb_init(prevz);
    arb_init(stepsize);

    acb_t s, t;
    acb_init(s);
    acb_init(t);

    roots = 0;
    arb_set(prevz, begpz);
    arb_sub(a, endp, begp, prec);
    arb_add_si(b, depth, 1, prec);
    arb_div(stepsize, a, b, prec);

    arb_one(w);
    while(arb_le(w, depth))
    {
        arb_mul(a, w, stepsize, prec);
        arb_add(a, a, begp, prec);

        acb_set_arb(t, a);
        acb_dirichlet_hardy_z(s, t, NULL, NULL, 1, zprec[id]);
        acb_get_real(currz,s);
        arb_add_si(arb_mat_entry(Zstat, id, 1), arb_mat_entry(Zstat, id, 1), 1, prec);

        arb_mul(a, currz, prevz, prec);
        if (arb_is_negative(a))
           roots++;

        arb_set(prevz, currz);

        arb_add_si(w, w, 1, prec);        
    }

    arb_mul(a, currz, endpz, prec);
    if (arb_is_negative(a))
       roots++;

    arb_clear(a);
    arb_clear(b);
    arb_clear(w);
    arb_clear(currz);
    arb_clear(prevz);
    arb_clear(stepsize);

    acb_clear(s);
    acb_clear(t);

    return(roots);
}

//Check if the expected zeros of Z(t) could be been found in the interval gvb..gve (0 = all roots have been found).
slong
CountAllRootZ(slong res, slong gvb, slong gve, slong gzb, slong gze, slong target, slong id, slong prec)
{
    arb_t depth, maxdepth;
    arb_init(depth);
    arb_init(maxdepth);

    slong blkroots;

    //count all zeros in this Gram interval
    arb_one(depth); blkroots = 0; arb_set_si(maxdepth, 100); 
    while(blkroots < target)
    {
        arb_add_si(depth, depth, 1, prec);
        blkroots = CountRootZ(blkroots, arb_mat_entry(currgv, id, gvb), arb_mat_entry(currgv, id, gve), arb_mat_entry(currgz, id, gzb), arb_mat_entry(currgz, id, gze), depth, id, prec);
        if(arb_gt(depth, maxdepth))
                break;
    }
    res = target - blkroots;

    arb_clear(depth);
    arb_clear(maxdepth);

    return(res);
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
    arb_set_si(b, 100000000000000);
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
    arb_t a, b, gv, gz;
    arb_init(a);
    arb_init(b);
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

    //maintain history of previous two Good Gram point values and associated Z-values (currx + 2 is the latest)
    if(gramquality == 1)
    {
       for(i = 0; i < 2; i++)
       {
            arb_set(arb_mat_entry(currgi, id, i), arb_mat_entry(currgi, id, i + 1));
            arb_set(arb_mat_entry(currgv, id, i), arb_mat_entry(currgv, id, i + 1));
            arb_set(arb_mat_entry(currgz, id, i), arb_mat_entry(currgz, id, i + 1));
        }
        arb_set(arb_mat_entry(currgi, id, 2), gi);
        arb_set(arb_mat_entry(currgv, id, 2), gv);
        arb_set(arb_mat_entry(currgz, id, 2), gz);
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(gv);
    arb_clear(gz);

    acb_clear(s);
    acb_clear(t);

    return(gramquality);
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
    zprec[256] = z_accuracy(prec, a, targetprec);

    //loop until a consecutive pattern of G G G Grampoints has been found and pick the middle one.
    arb_set(a, d);
    arb_add_si(b, a, 1, prec);
    arb_add_si(c, b, 1, prec);
    while ((GoodGram(gramquality, a, 256, prec)==0) || (GoodGram(gramquality, b, 256, prec)==0)
	        || (GoodGram(gramquality, c, 256, prec)==0))
    {
        arb_add_si(a, a, 1, prec);
        arb_add_si(b, a, 1, prec);
        arb_add_si(c, b, 1, prec);
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
    //recover the data passed to this specific thread
    struct ThreadData* data=voidData;
    slong id=data->id;
    slong prec=data->prec;

    arb_t a, gi, gf;
    arb_init(a);
    arb_init(gi);
    arb_init(gf);

    slong tmp, blocklength, prevblklen, targetprec, gramquality, prevqual, roots, missingzeros;

    //get start and end of the range for this thread
    arb_set(gi, arb_mat_entry(threadworksplit, id, 0));
    arb_set(gf, arb_mat_entry(threadworksplit, id, 1));

    gramquality = 0; prevqual = 1; blocklength = 0; missingzeros = 0;

    //determine the required decimal precision (3 digits) for the start value Z(gram(gi))
    targetprec = 3 * 3.32192809488736 + 5;
    gram(a, gf, prec);
    zprec[id] = z_accuracy(prec, a, targetprec);

    while (arb_lt(gi, gf))
    {
       arb_add_si(gi, gi, 1, prec);

       //establish Gram quality
       gramquality = GoodGram(gramquality, gi, id, prec);

       //Good Gram point followed by Good
       if(prevqual && gramquality)
       {    
           blocklength = 0;
           prevblklen = blocklength;
           arb_add_si(arb_mat_entry(gramblockcounts, id, blocklength), arb_mat_entry(gramblockcounts, id, blocklength), 1, prec);
           arb_add_si(arb_mat_entry(Rootsfound, id, 0), arb_mat_entry(Rootsfound, id, 0), 1, prec);

           //Type I Rosser's rule violation in previous block? Check if any missing zeros can be found now. 
           if(missingzeros > 0)
           {
               tmp = missingzeros;
               missingzeros = CountAllRootZ(missingzeros, 1, 2, 1, 2, missingzeros + 1, id, prec);
               if(missingzeros == 0)
               {
                   arb_add_si(arb_mat_entry(Rootsfound, id, 2), arb_mat_entry(Rootsfound, id, 2), tmp, prec);
               }
           }
       }

       //Good Gram point followed by Bad (start of a new Gram block)
       if(prevqual && gramquality == 0)
           blocklength = 1;

       //Bad Gram point followed by Bad (inner part of a Gram Block)
       if(prevqual== 0 && gramquality == 0)
           blocklength = blocklength + 1;

       //Bad Gram point followed by Good (end of a Gram block)
       if(prevqual== 0 && gramquality)
       {
           arb_add_si(arb_mat_entry(gramblockcounts, id, blocklength), arb_mat_entry(gramblockcounts, id, blocklength), 1, prec);

           //Type III Rosser's rule violation in previous block? Check wether any missing zeros can be found now. 
           if(missingzeros > 0)
           {
               tmp = missingzeros;
               missingzeros = CountAllRootZ(missingzeros, 1, 2, 1, 2, blocklength + missingzeros + 1, id, prec);
               if(missingzeros == 0)
               {
                   arb_add_si(arb_mat_entry(Rootsfound, id, 4), arb_mat_entry(Rootsfound, id, 4), tmp, prec);
               }
               else
               {
                   printf("!!! Failed to find the missing zeros of Rosser's rule violation of Gram point :");
                   arb_printd(arb_mat_entry(currgi, id, 0), 20);
                   printf(" !!! \n Please inspect manually! \n");
                   missingzeros = 0;
               }
           }

           //Check the current bad Gram block or establish any missing zeros (indicating a violation of Rosser's rule)
           missingzeros = CountAllRootZ(missingzeros, 1, 2, 1, 2, blocklength + 1, id, prec);

           //Type II Rosser's rule violation in current block? Check if missing zeros can be found in the preceding block. 
           if(missingzeros == 0)
                 arb_add_si(arb_mat_entry(Rootsfound, id, 1), arb_mat_entry(Rootsfound, id, 1), blocklength + 1, prec);
           else
           {
                tmp = missingzeros;
                missingzeros = CountAllRootZ(missingzeros, 0, 1, 0, 1, prevblklen + missingzeros + 1, id, prec);
                if(missingzeros == 0)
                {
                   arb_add_si(arb_mat_entry(Rootsfound, id, 3), arb_mat_entry(Rootsfound, id, 3), tmp, prec);
                }
           }
 
           //printf("Number of roots fixed: %ld requiring for blocklength %ld: \n", blocklength + 1, blocklength + 1);

           roots = roots + blocklength + 1;
           prevblklen = blocklength;
       }

       prevqual = gramquality;
    }

    arb_clear(a);
    arb_clear(gi);
    arb_clear(gf);

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
    for (j = 0; j < 2; j++)
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
                 printf("to assess a Gram point : %ld \n", tmp); tmp1 = tmp; break;
              case 1:
                 printf("to fix bad Gram points : %ld \n", tmp); break;
              default:
                 break;
           }
           arb_add(totcount, totcount, subcount, prec);
       }
    }
    arb_get_ubound_arf(u, totcount, prec);
    tmp=arf_get_si(u, ARF_RND_NEAR);
    printf("Total evaluations of Z(t)  ======================== : %ld \n\n", tmp);

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
                 printf("Good to Good interval : %ld \n", tmp); break;
              case 1:
                 printf("Gramblocks G B...B G  : %ld \n", tmp); break;
              case 2:
                 printf("next Good to Good int.: %ld (Rosser violation Type I)\n", tmp); break;
              case 3:
                 printf("previous interval     : %ld (Rosser violation Type II)\n", tmp); break;
              case 4:
                 printf("next Gramblock        : %ld (Rosser violation Type III)\n", tmp); break;
              default:
                 break;
           }
           arb_add(totcount, totcount, subcount, prec);
       }
    }
    arb_get_ubound_arf(u, totcount, prec);
    tmp=arf_get_si(u, ARF_RND_NEAR);
    printf("Total roots of Z(t)  ======================= : %ld \n\n", tmp);

    arb_clear(subcount);

    arf_clear(u);
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
    slong gramqual, numclients, clientid, numthreads, prec, tmp, tmp1; gramqual = 0;

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
    arb_mat_init(currgi, 257, 3);
    arb_mat_init(currgv, 257, 3);
    arb_mat_init(currgz, 257, 3);
    arb_mat_init(Zstat, 257, 2);

TIMEIT_ONCE_START

    zprec[256] = prec;
    if(GoodGram(gramqual, xa, 256, prec) == 0)
    {
         printf("\n ! The starting value chosen is a bad Gram point. Please select a good one !\n\n");
         goto finish;
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

    Count_number_of_roots_on_line(rootsonline, numthreads, prec);

    arb_get_ubound_arf(u, arb_mat_entry(clientworksplit, clientid, 0), prec);
    tmp=arf_get_si(u, ARF_RND_NEAR);
    arb_get_ubound_arf(u, arb_mat_entry(clientworksplit, clientid, 1), prec);
    tmp1=arf_get_si(u, ARF_RND_NEAR);
    printf("In the range of Gram points from %ld to %ld, \n", tmp, tmp1);
    arb_get_ubound_arf(u, rootsinstrip, prec);
    tmp=arf_get_si(u, ARF_RND_NEAR);
    arb_get_ubound_arf(u, rootsonline, prec);
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
    arb_mat_clear(currgi);
    arb_mat_clear(currgv);
    arb_mat_clear(currgz);
    arb_mat_clear(Zstat);

    arb_clear(xa);
    arb_clear(xb);
    arb_clear(res);
    arb_clear(rootsinstrip);
    arb_clear(rootsonline);

    arf_clear(u);
 
    flint_cleanup();

    return result;
}
