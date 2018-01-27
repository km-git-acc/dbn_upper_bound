import math
from math import *
from utility import *
from Ht_real_compute import *
from scipy.optimize import *
from first_1000_zeta import *

first_1000_H0_zero = [x*2 for x in first_1000_zeta]

t=0.2
for i,nearby_root  in enumerate(first_1000_H0_zero): 
      Ht_root = fsolve(Ht_real,nearby_root,args=(t,))[0]
      print (t,i,Ht_root)
