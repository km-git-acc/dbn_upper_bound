import scipy
from scipy.integrate import *
import numpy as np
from cmath import *

#pi upto 20 decimal places. add further precision if necessary
PI = 3.14159265358979323846
PI_sq = PI*PI


def Nt(t, T):
      """
      Evaluates equation (4) in Terry's blog
      https://terrytao.wordpress.com/2018/01/27/polymath15-first-thread-computing-h_t-asymptotics-and-dynamics-of-zeroes/
      :param t: the "time" parameter
      :param T: height
      :return: right side of (4) in the blog link
      """
      T=float(T)
      T_by_4PI = T/(4*PI)
      return T_by_4PI*log(T_by_4PI) - T_by_4PI + (t/16.0)*log(T_by_4PI)

def phi_decay(u,n_max=100):
    running_sum=0
    for n in range(1,n_max+1):
        term1=2*PI_sq*pow(n,4)*exp(9*u) - 3*PI*pow(n,2)*exp(5*u)
        term2=exp(-1*PI*pow(n,2)*exp(4*u))
        running_sum += term1*term2
        #print n,term1, term2, running_sum
    return running_sum

def Ht_complex_integrand(u, z, t):
     return exp(t*u*u)*phi_decay(u)*cos(z*u)

def Ht_complex(z,t):
    #may work well only for small to medium values of z  
    def real_func(a,b,c): return scipy.real(Ht_complex_integrand(a,b,c))
    def imag_func(a,b,c): return scipy.imag(Ht_complex_integrand(a,b,c))   
    real_part = quad(real_func, 0, 10, args=(z,t))
    imag_part = quad(imag_func, 0, 10, args=(z,t))
    return (real_part[0] + 1j*imag_part[0], real_part[1], imag_part[1])

def Ht_real_integrand(u, z, t):
     if(abs(scipy.imag(z))>0 or abs(scipy.imag(t))>0 or abs(scipy.imag(u))>0): 
        print ("complex values not allowed for this function")
        return "error"
     u,z,t = scipy.real(u),scipy.real(z),scipy.real(t)
     return scipy.real(exp(t*u*u)*phi_decay(u)*cos(z*u))

def Ht_real(z,t):
     if(abs(scipy.imag(z))>0 or abs(scipy.imag(t))>0): 
        print ("complex values not allowed for this function")
        return "error"
     z,t = scipy.real(z),scipy.real(t)
     #return quad(Ht_real_integrand, 0, np.inf, args=(z,t)) #causing overflow errors so np.inf replaced with 10
     return quad(Ht_real_integrand, 0, 10, args=(z,t))

#check phi_decay values
'''print (phi_decay(0.001))
print (phi_decay(0.01))
print (phi_decay(0.1))
print (phi_decay(0.5))
print (phi_decay(1))'''
