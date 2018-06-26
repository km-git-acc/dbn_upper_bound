The scripts in this folder are written in Pari/GP. To install it, please follow the steps <a href='https://pari.math.u-bordeaux.fr/'>here</a>

It is recommended to first understand the theory and derivations on which the below guide and the scripts are based. Towards that, please read <a href='https://terrytao.wordpress.com/'>Prof. Tao's blog</a>, and the <a href='https://github.com/km-git-acc/dbn_upper_bound/blob/master/Writeup/debruijn.pdf'>writeup</a>. 

General steps followed to prove a certain dbn bound 
-----------------------------------------------------------------------------------------------
1) For a given t0 and y0, choose an N_min where some euler bound, mollified or otherwise is reasonably positive
2) For the N_min determined in 1), choose a barrier region [X,X+1] where |(A+B)/B0| seems large
3) Calculate the winding number of |A+B|/B0 for the rectangle [[X,X+1],[y0,1]], moving in small increments of t from 0 to t0
4) If all the winding numbers in 3) are zero, choose an Nmax (with y0 and t0) where the non mollified triangle bound, or its more conservative tail integral version is positive.
5) Compute the euler bounds or suitably modified versions for each N in [Nmin,Nmax]
6) If all the computed bounds in 5) are reasonably positive, calculate the dbn bound using y0 and t0 (which will be conditional to RH verified atleast upto height X/2)

(In 3), 4) and 5), also conduct a check whether the error bounds for |(H-(A+B))/B0| are smaller than the smallest |(A+B)/B0| and the smallest euler bounds)

Guide to using the scripts
-------------------------------------------
-> Euler bounds can be calculated through the functions in the abbeff_largex_bounds file. Within that file, 
for 1) and 4), the function abbeff_largex_ep_bound can be used, 
while for 5) where [Nmin,Nmax] can be a large range, functions abbeff_largex_ep_sawtooth_base_lbound and abbeff_largex_ep_sawtooth_incremental_lbounds can be used
-> For barrier strip related calculations, the file barrier_multieval_t_agnostic can be used
-> To compute error bounds, the file error bounds can be used

The remaining scripts contain definitions and work related to earlier approaches to establishing dbn bounds, which one is encouraged to explore as well. 
