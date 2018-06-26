The scripts in this folder are written in Pari/GP. To install it, please follow the steps <a href='https://pari.math.u-bordeaux.fr/'>here</a>

General steps followed to prove a certain dbn bound 
-----------------------------------------------------------------------------------------------
1) For a given t0 and y0, choose an N_min where some euler bound, mollified or otherwise is reasonably positive
2) For the N_min determined in 1), choose a barrier region [X,X+1] where |A+B|/B0 seems large
3) Calculate the winding number of |A+B|/B0 for the rectangle [[X,X+1],[y0,1]], moving in small increments of t from 0 to t0
4) If all the winding numbers in 3) are zero, choose an Nmax (with y0 and t0) where the non mollified triangle bound, or its more conservative tail integral version is positive.
5) Compute the euler bounds or suitably modified versions for each N in [Nmin,Nmax]
6) If all the computed bounds in 5) are reasonably positive, calculate the dbn bound using y0 and t0 (which will be conditional to RH verified atleast upto height X/2)

