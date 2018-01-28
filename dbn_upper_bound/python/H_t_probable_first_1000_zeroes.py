import sys

orig_stdout = sys.stdout
f = open('out_600_01.txt', 'w')
sys.stdout = f

import math
from math import *
from utility import *
from Ht_real_compute import *
from scipy.optimize import *
from first_1000_zeta import *

first_1000_H0_zero = [x*2 for x in first_1000_zeta]

for j in range(0,50):
 t=.01*1
 for i,nearby_root  in enumerate(first_1000_H0_zero[600:700]): 
       Ht_root = fsolve(Ht_real,nearby_root,args=(t,))[0]
       print (t,i,Ht_root)
	   
	   
sys.stdout = orig_stdout
f.close()