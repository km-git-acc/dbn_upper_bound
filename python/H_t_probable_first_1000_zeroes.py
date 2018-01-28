import sys

orig_stdout = sys.stdout
f = open('out_1000_01.txt', 'w')
sys.stdout = f

import math
from math import *
from utility import *
from scipy.optimize import *
from first_1000_zeta import *

first_1000_H0_zero = [x*2 for x in first_1000_zeta]

#k is a stopping parameter and is the number of zeros to calculate up to set to 1000 by default
k=1000
for j in range(1,50):
 t=.01*j
 for i,nearby_root  in enumerate(first_1000_H0_zero[1:k]): 
       Ht_root = fsolve(Ht_real,nearby_root,args=(t,))[0]
       print (t,i,Ht_root)
	   
	   
sys.stdout = orig_stdout
f.close()