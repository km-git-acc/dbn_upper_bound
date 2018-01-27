import math
from math import *

#pi upto 20 decimal places. add further precision if necessary
PI = 3.14159265358979323846
PI_sq = PI*PI

#print PI
#print PI_sq

def Nt(t,T):
      T=float(T)
      T_by_4PI = T/(4*PI)
      return T_by_4PI*ln(T_by_4PI) - T_by_4PI + t*ln(T) 

def phi_decay(u,n_max=100):
    running_sum=0
    for n in range(1,n_max+1):
        term1=2*PI_sq*pow(n,4)*exp(9*u) - 3*PI*pow(n,2)*exp(5*u)
        term2=exp(-1*PI*pow(n,2)*exp(4*u))
        running_sum += term1*term2
        #print n,term1, term2, running_sum
    return running_sum

#check phi_decay values
print phi_decay(0.001)
print phi_decay(0.01)
print phi_decay(0.1)
print phi_decay(0.5)
print phi_decay(1)
