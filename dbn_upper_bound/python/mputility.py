'''Using mpmath instead of numpy and scipy to arbitray precision
Could be quite slow to compute but much more accurate'''

from dbn_upper_bound.python.constants import PI
from mpmath import mp, mpf, mpc, exp, power, log, cos, sqrt, quad, findroot

mp.dps=20

PI = mpf(PI)
PI_sq = PI*PI


def phi_decay(u,n_max=100):
    running_sum=0
    u=mpc(u)
    for n in range(1,n_max+1):
        term1=2*PI_sq*power(n,4)*exp(9*u) - 3*PI*power(n,2)*exp(5*u)
        term2=exp(-1*PI*power(n,2)*exp(4*u))
        running_sum += term1*term2
    return running_sum

def Ht_complex_integrand(u, z, t):
     u,z,t=mpc(u),mpc(z),mpc(t)
     return exp(t*u*u)*phi_decay(u)*cos(z*u)

def Ht_complex(z,t):
    #may work well only for small to medium values of z  
    part = quad(lambda u: Ht_complex_integrand(u,z,t), [0, 10])
    return part

def Ht_complex_root_finder(complex_guess,t):
    result = findroot(lambda z: Ht_complex(z,t),mpc(complex_guess))
    return result

