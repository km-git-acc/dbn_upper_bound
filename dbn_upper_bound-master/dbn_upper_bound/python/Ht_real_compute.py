'''
from scipy.integrate import quad
import numpy as np
import math
from math import *
from utility import *

def Ht_real_integrand(u, z, t):
     return exp(t*u*u)*phi_decay(u)*cos(z*u)

def Ht_real(z,t):
     #return quad(Ht_real_integrand, 0, np.inf, args=(z,t)) #causing overflow errors so np.inf replaced with 10
     return quad(Ht_real_integrand, 0, 10, args=(z,t))
'''
'''
The above part is no longer necessary since the functions are present in utility.py. 
It could be that the math functions are faster than the cmath functions in which case 
we may have to make further changes to speed up Ht_real computations'''

'''The snippet below are simple loops for evaluating Ht_real on multiple points. This can be reused if required'''

#check values of Ht_real
#for t in range(0,6):
 #     t=float(t)/10.0
  #    print(t)
   #   for i in range(0,10000):
    #        z=float(i)/100
     #       print(t,z, Ht_real(z,t))
