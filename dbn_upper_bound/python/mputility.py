'''This module has commonly used functions computed using the mpmath library, hence arbitray precision, 
instead of using scipy as in utility.py. For explanations of the first few functions, please check utility.py. The new ones are explained below.
Default precision is set to 30 significant digits below, but this can be overriden in your code.
If you can get the gmpy2 library installed, that could provide a major speed boost for mpmath'''

'''The new functions (not present in utility.py) are 
1) The main and helper functions for computing the Riemann Siegel Z function used to compute zeta zeroes (for benchmarking purposes),
2) The approx functional eqn Ht_AFE and helper functions for computing H_t as derived here 
http://michaelnielsen.org/polymath1/index.php?title=Asymptotics_of_H_t''' 

from mpmath import mp
from dbn_upper_bound.python.constants import PI

mp.dps=30

PI = mp.mpf(PI)
PI_sq = PI*PI


def phi_decay(u,n_max=100):
    running_sum=0
    u=mp.mpc(u)
    for n in range(1,n_max+1):
        term1=2*PI_sq*mp.power(n,4)*mp.exp(9*u) - 3*PI*mp.power(n,2)*mp.exp(5*u)
        term2=mp.exp(-1*PI*mp.power(n,2)*mp.exp(4*u))
        running_sum += term1*term2
        #print n,term1, term2, running_sum
    return running_sum

def Ht_complex_integrand(u, z, t):
     u,z,t=mp.mpc(u),mp.mpc(z),mp.mpc(t)
     return mp.exp(t*u*u)*phi_decay(u)*mp.cos(z*u)

def Ht_complex(z,t):
    #may work well only for small to medium values of z  
    part = mp.quad(lambda u: Ht_complex_integrand(u,z,t), [0, 10])
    return part

def Ht_complex_root_finder(complex_guess,t):
    result = mp.findroot(lambda z: Ht_complex(z,t),mp.mpc(complex_guess))
    return result

'''Begin Riemann Siegel Z block'''

def RStheta(x):
    return (x/2)*(mp.log(x/(2*PI)) - 1) - PI/8 + 1/(48*x) + 7/(5760*mp.power(x,3)) 

def c0(p):
   return mp.cos(2*PI*(p*p - p - 1/16))/mp.cos(2*PI*p)

def c1(p):
   return (-1/(96*PI*PI))*mp.diff(lambda x: c0(x), p, 3)

def RSZ_plain(x):
    x = mp.mpf(x.real)
    tau = mp.sqrt(x/(2*PI))
    N = int(tau.real)
    z = 2*(x-N)-1
    running_sum=0
    for n in range(1,N+1):
        running_sum += mp.cos(RStheta(x) - x*mp.log(n))/mp.sqrt(n)
    return (2*running_sum).real

def RSZ_upto_c0(x):
    x = mp.mpf(x.real)
    tau = mp.sqrt(x/(2*PI))
    N = int(tau.real)
    p = tau-N
    running_sum=0
    for n in range(1,N+1):
        running_sum += mp.cos(RStheta(x) - x*mp.log(n))/mp.sqrt(n)
    return (2*running_sum + mp.power(-1,N-1)*mp.power(tau,-0.5)*c0(p)).real

def RSZ_upto_c1(x):
    x = mp.mpf(x.real)
    tau = mp.sqrt(x/(2*PI))
    N = int(tau.real)
    p = tau-N
    running_sum=0
    for n in range(1,N+1):
        running_sum += mp.cos(RStheta(x) - x*mp.log(n))/mp.sqrt(n)
    return (2*running_sum + mp.power(-1,N-1)*mp.power(tau,-0.5)*(c0(p) - (1/tau)*c1(p))).real

'''End Riemann Siegel Z block'''

'''Begin H_t approx functional eqn (Ht_AFE) block'''

def F0(s):
    # works only for x > 14 approx due to the sqrt(x/4PI) factor
    term1 = mp.power(PI,-1*s/2)*mp.gamma((s+4)/2)
    x=2*s.imag
    N = int(mp.sqrt(x/(4*PI)))
    running_sum=0
    for n in range(1,N+1):
        running_sum += 1/mp.power(n,s)
    return term1*running_sum

def H0_AFE(z):
    z = mp.mpc(z)
    f_arg1 = (1 + 1j*z.real - z.imag)/2
    f_arg2 = (1 + 1j*z.real + z.imag)/2
    if z.imag==0:
        result= 0.5*(F0(f_arg1).real)
    else: 
        result = 0.25*(F0(f_arg1) + F0(f_arg2).conjugate())
    return result

def Ft(s,t):
    # works only for x > 14 approx due to the sqrt(x/4PI) factor
    term1 = mp.power(PI,-1*s/2)*mp.gamma((s+4)/2)
    x=2*s.imag
    N = int(mp.sqrt(x/(4*PI)))
    running_sum=0
    for n in range(1,N+1):
        running_sum += mp.exp((t/16)*mp.power(mp.log((s+4)/(2*PI*n*n)),2))/mp.power(n,s) # main eqn 
    return term1*running_sum

def Ht_AFE(z,t):
    z,t = mp.mpc(z),mp.mpc(t)
    f_arg1 = (1 + 1j*z.real - z.imag)/2
    f_arg2 = (1 + 1j*z.real + z.imag)/2
    if z.imag==0:
        result= 0.5*(Ft(f_arg1,t).real)
    else: 
        result = 0.25*(Ft(f_arg1,t) + Ft(f_arg2,t).conjugate())
    return result

'''End H_t approx functional eqn (Ht_AFE) block'''
