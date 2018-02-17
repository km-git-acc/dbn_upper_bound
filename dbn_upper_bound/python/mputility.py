"""
This module has commonly used functions computed using the
mpmath library, hence arbitray precision,
instead of using scipy as in utility.py. For explanations of the
first few functions, please check utility.py. The new ones are explained below.
Default precision is set to 30 significant digits below, but
this can be overriden in your code. If you can get the gmpy2 library
installed, that could provide a major speed boost for mpmath
The new functions (not present in utility.py) are
1) The main and helper functions for computing the Riemann Siegel
    Z function used to compute zeta zeroes (for benchmarking purposes),
2) The approx functional eqn Ht_AFE and helper functions for
    computing H_t as derived here
http://michaelnielsen.org/polymath1/index.php?title=Asymptotics_of_H_t
"""

import csv
from mpmath import mp

mp.dps = 30


def phi_decay(u, n_max=100):
    """
    Computes Phi(u) in Terry' blog at
    https://terrytao.wordpress.com/2018/02/02/polymath15-second-thread-generalising-the-riemann-siegel-approximate-functional-equation/
    :param u: input complex number
    :param n_max: upper limit for summation. It has to be a positive
    integer
    :return: Phi(u)
    """
    running_sum = 0
    u = mp.mpc(u)
    for n in range(1, n_max+1):
        term1 = 2 * mp.pi() * mp.pi() * mp.power(n, 4) * mp.exp(9*u) \
                - 3 * mp.pi() * mp.power(n, 2) * mp.exp(5*u)
        term2 = mp.exp(-1*mp.pi() * mp.power(n, 2) * mp.exp(4*u))
        running_sum += term1 * term2

    return running_sum


def Ht_complex_integrand(u, z, t):
    """
    Computes the integrand of \H_t(u) in Terry' blog at
    https://terrytao.wordpress.com/2018/02/02/polymath15-second-thread-generalising-the-riemann-siegel-approximate-functional-equation/
    :param u: integration parameter
    :param z: point at which H_t is computed
    :param t: the "time" parameter
    :return: integrand of H_t(z)
    """
    u, z, t = mp.mpc(u), mp.mpc(z), mp.mpc(t)
    return mp.exp(t * u * u) * phi_decay(u) * mp.cos(z * u)


def Ht_complex(z, t):
    #may work well only for small to medium values of z  
    part = mp.quad(lambda u: Ht_complex_integrand(u, z, t), [0, 10])
    return part

def Ht_complex_root_finder(complex_guess, t):
    result = mp.findroot(lambda z: Ht_complex(z, t), mp.mpc(complex_guess))
    return result

'''Begin Riemann Siegel Z block'''

def RStheta(x):
    return (x/2)*(mp.log(x/(2*mp.pi())) - 1) - mp.pi()/8 + 1/(48*x) + 7/(5760*mp.power(x,3)) 

def c0(p):
    return mp.cos(2*mp.pi()*(p*p - p - 1/16))/mp.cos(2*mp.pi()*p)

def c1(p):
    return (-1/(96*mp.pi()*mp.pi()))*mp.diff(lambda x: c0(x), p, 3)

def RSZ_plain(x):
    x = mp.mpf(x.real)
    tau = mp.sqrt(x/(2*mp.pi()))
    N = int(tau.real)
    z = 2*(x-N)-1
    running_sum = 0
    for n in range(1, N+1):
        running_sum += mp.cos(RStheta(x) - x*mp.log(n))/mp.sqrt(n)
    return (2 * running_sum).real

def RSZ_upto_c0(x):
    x = mp.mpf(x.real)
    tau = mp.sqrt(x / (2 * mp.pi()))
    N = int(tau.real)
    p = tau - N
    running_sum = 0
    for n in range(1, N+1):
        running_sum += mp.cos(RStheta(x) - x * mp.log(n)) / mp.sqrt(n)
    return (2 * running_sum + mp.power(-1,N-1) * mp.power(tau,-0.5) * c0(p)).real

def RSZ_upto_c1(x):
    x = mp.mpf(x.real)
    tau = mp.sqrt(x / (2 * mp.pi()))
    N = int(tau.real)
    p = tau - N
    running_sum = 0
    for n in range(1, N+1):
        running_sum += mp.cos(RStheta(x) - x * mp.log(n)) / mp.sqrt(n)
    return (2 * running_sum + mp.power(-1,N-1) * mp.power(tau,-0.5) * (c0(p) - (1/tau) * c1(p))).real

'''End Riemann Siegel Z block'''

'''Begin H_t approx functional eqn (Ht_AFE) block'''

def F0(s):
    # works only for x > 14 approx due to the sqrt(x/4PI) factor
    term1 = mp.power(mp.pi(),-1*s/2)*mp.gamma((s+4)/2)
    x = 2 * s.imag
    N = int(mp.sqrt(x/(4*mp.pi())))
    running_sum = 0
    for n in range(1, N+1):
        running_sum += 1/mp.power(n,s)
    return term1 * running_sum

def H0_AFE(z):
    z = mp.mpc(z)
    f_arg1 = (1 + 1j*z.real - z.imag)/2
    f_arg2 = (1 + 1j*z.real + z.imag)/2
    if z.imag == 0:
        result = 0.5 * (F0(f_arg1).real)
    else: 
        result = 0.25 * (F0(f_arg1) + F0(f_arg2).conjugate())
    return result

def Ft(s, t):
    # works only for x > 14 approx due to the sqrt(x/4PI) factor
    term1 = mp.power(mp.pi(),-1*s/2) * mp.gamma((s+4)/2)
    x = 2 * s.imag
    N = int(mp.sqrt(x/(4 * mp.pi())))
    running_sum = 0
    for n in range(1, N+1):
        running_sum += mp.exp((t/16) * mp.power(mp.log((s+4) / (2 * mp.pi() * n * n)),2)) / mp.power(n,s) # main eqn 
    return term1 * running_sum

def Ht_AFE(z, t):
    z, t = mp.mpc(z), mp.mpc(t)
    f_arg1 = (1 + 1j*z.real - z.imag)/2
    f_arg2 = (1 + 1j*z.real + z.imag)/2
    if z.imag == 0:
        result = 0.5 * (Ft(f_arg1,t).real)
    else: 
        result = 0.25 * (Ft(f_arg1,t) + Ft(f_arg2,t).conjugate())
    return result

def psi(p):
    term1 = 2 * mp.pi()
    term2 = mp.cos(mp.pi() * (p*p/2 - p - 1/8))
    term3 = mp.exp(0.5j * mp.pi() * p * p - 5j * mp.pi()/8)
    term4 = mp.cos(mp.pi() * p)
    return term1 * term2 * term3 / term4


def Ht_AFE_A(z, t):
    """
    This is the much more accurate approx functional eqn posted by Terry at
    https://terrytao.wordpress.com/2018/02/02/polymath15-second-thread-generalising-the-riemann-siegel-approximate-functional-equation/#comment-492182
    :param z: point at which H_t is computed
    :param t: the "time" parameter
    :return: the A part in Ht
    """
    z, t = mp.mpc(z), mp.mpc(t)
    s = (1 + 1j * z.real - z.imag)/2
    tau = mp.sqrt(s.imag/(2 * mp.pi()))
    N = int(tau)

    A_pre = (1/16) * s * (s-1) \
            * mp.power(mp.pi(), -1*s/2) * mp.gamma(s/2)
    A_sum = 0.0
    for n in range(1, N+1):
        if t.real > 0:
            A_sum += mp.exp((t/16) * mp.power(mp.log((s+4)/(2*mp.pi()*n*n)), 2))/mp.power(n, s)
        else:
            A_sum += 1/mp.power(n, s)

    return A_pre * A_sum


def Ht_AFE_B(z, t):
    """
    This is the much more accurate approx functional eqn posted by Terry at
    https://terrytao.wordpress.com/2018/02/02/polymath15-second-thread-generalising-the-riemann-siegel-approximate-functional-equation/#comment-492182
    :param z: point at which H_t is computed
    :param t: the "time" parameter
    :return: the B part in Ht
    """
    z, t = mp.mpc(z), mp.mpc(t)
    s = (1 + 1j * z.real - z.imag) / 2
    tau = mp.sqrt(s.imag/(2*mp.pi()))
    M = int(tau)

    B_pre = (1/16) * s * (s-1) * mp.power(mp.pi(), (s-1)/2) * mp.gamma((1-s)/2)
    B_sum = 0.0
    for m in range(1, M+1):
        if t.real > 0:
            B_sum += mp.exp((t/16) * mp.power(mp.log((5-s)/(2 * mp.pi() * m * m)), 2))/mp.power(m, 1-s)
        else:
            B_sum += 1/mp.power(m, 1-s)

    return B_pre * B_sum


def Ht_AFE_C(z, t):
    """
    This is the much more accurate approx functional eqn posted by Terry at
    https://terrytao.wordpress.com/2018/02/02/polymath15-second-thread-generalising-the-riemann-siegel-approximate-functional-equation/#comment-492182
    :param z: point at which H_t is computed
    :param t: the "time" parameter
    :return: the C part in Ht
    """
    z, t = mp.mpc(z), mp.mpc(t)
    s = (1 + 1j * z.real - z.imag) / 2
    tau = mp.sqrt(s.imag/(2 * mp.pi()))
    N = int(tau)
    M = int(tau)

    C_pre1 = -(1/16) * s * (s-1) * mp.power(mp.pi(), -1 * s/2) * mp.gamma(s/2)
    C_pre2 = mp.gamma(1-s) * mp.exp(-1j * mp.pi() * s)/(2j * mp.pi())
    C_pre3 = mp.power(2j * mp.pi() * M, s-1) * mp.exp(-1 * t * mp.pi() * mp.pi()/64)
    C_psi = psi((s/(2j * M * mp.pi())) - N)

    return C_pre1 * C_pre2 * C_pre3 * C_psi


def Ht_AFE_ABC(z, t):
    """
    This is the much more accurate approx functional eqn posted by Terry at
    https://terrytao.wordpress.com/2018/02/02/polymath15-second-thread-generalising-the-riemann-siegel-approximate-functional-equation/#comment-492182
    :param z: point at which H_t is computed
    :param t: the "time" parameter
    :return: Ht as a sum of A, B, C
    """
    A = Ht_AFE_A(z, t)
    B = Ht_AFE_B(z, t)
    C = Ht_AFE_C(z, t)

    H = A+B+C
    if z.imag == 0:
        return H.real
    else:
        return H

def AB_analysis(z,t):
    z,t = mp.mpc(z),mp.mpc(t)
    s = (1 + 1j*z.real - z.imag)/2
    tau = mp.sqrt(s.imag/(2*mp.pi()))
    N = int(tau)
    M = int(tau)
    
    A_pre = (1/16) * s * (s-1) * mp.power(mp.pi(),-1*s/2) * mp.gamma(s/2)
    A_sum = 0.0
    for n in range(1,N+1):
        if t.real > 0: 
            A_sum += mp.exp((t/16)*mp.power(mp.log((s+4)/(2*mp.pi()*n*n)),2)) / mp.power(n,s)
        else: 
            A_sum += 1/mp.power(n,s) 
    A = A_pre*A_sum
    
    B_pre = (1/16) * s * (s-1) * mp.power(mp.pi(),(s-1)/2) * mp.gamma((1-s)/2)
    B_sum = 0.0
    for m in range(1,M+1):
        if t.real > 0: 
            B_sum += mp.exp((t/16)*mp.power(mp.log((5-s)/(2*mp.pi()*m*m)),2)) / mp.power(m,1-s)
        else: 
            B_sum += 1/mp.power(m,1-s) 
    B = B_pre * B_sum
    
    B0 = B_pre * mp.exp((t/16)*mp.power(mp.log((5-s)/(2*mp.pi())),2))
    
    AplusB = A + B
    ABB0 = AplusB/B0
    return (AplusB, B0, ABB0, abs(ABB0)) 

def adjusted_AB_analysis(z,t):
    z,t = mp.mpc(z),mp.mpc(t)
    s = (1 + 1j*z.real - z.imag)/2
    tau = mp.sqrt(s.imag/(2*mp.pi()))
    N = int(tau)
    M = int(tau)
    
    s_like1 = (s+4)/2
    s_like2 = (5-s)/2
    initial_term = (1/4) * mp.sqrt(2*mp.pi())
    
    A_pre = initial_term * mp.power(mp.pi(),-1*s/2) * mp.exp((s_like1-0.5) * mp.log(s_like1)-s_like1)
    A_sum = 0.0
    for n in range(1,N+1):
        if t.real > 0: 
            A_sum += mp.exp((t/16) * mp.power(mp.log(s_like1/(mp.pi()*n*n)),2)) / mp.power(n,s)
        else: 
            A_sum += 1/mp.power(n,s) 
    A = A_pre * A_sum
    
    B_pre = initial_term * mp.power(mp.pi(),(s-1)/2) * mp.exp((s_like2-0.5) * mp.log(s_like2)-s_like2)
    B_sum = 0.0
    ddxB_sum = 0.0
    for m in range(1,M+1):
        if t.real > 0: 
            mth_term = mp.exp((t/16) * mp.power(mp.log(s_like2/(mp.pi()*m*m)),2)) / mp.power(m,1-s) 
        else: 
            mth_term = 1/mp.power(m,1-s)
        B_sum += mth_term
        ddxB_sum += mth_term * mp.log(m) 
    B = B_pre * B_sum
    B_common_term = mp.exp((t/16) * mp.power(mp.log(s_like2/mp.pi()),2))    
    B0 = B_pre * B_common_term
    ddxBB0 = 0.5j * ddxB_sum / B_common_term
    
    AplusB = A + B
    ABB0 = AplusB/B0
    return (AplusB, B0, ABB0, abs(ABB0), ddxBB0, abs(ddxBB0)) 

def Ht_AFE_ADJ_AB(z,t):
    '''This uses the adjusted A'+B' estimate to compute H_t for moderate to large T, where the C term is small.
    Also, there are some rearrangements to the A' and B' formulas to make them faster to compute''' 
    z,t = mp.mpc(z),mp.mpc(t)
    s = (1 + 1j*z.real - z.imag)/2
    tau = mp.sqrt(s.imag/(2*mp.pi()))
    N = int(tau)
    M = int(tau)
    
    s_like1 = (s+4)/2
    s_like2 = (5-s)/2
    log_sp1 = mp.log(s_like1/mp.pi())
    log_sp2 = mp.log(s_like2/mp.pi())
    initial_term = (1/4) * mp.sqrt(2*mp.pi())
    A0 = initial_term * mp.power(mp.pi(),-1*s/2) * mp.exp((s_like1-0.5)*mp.log(s_like1) - s_like1 + (t/16)*mp.power(log_sp1,2)) 
    B0 = initial_term * mp.power(mp.pi(),(s-1)/2) * mp.exp((s_like2-0.5)*mp.log(s_like2) - s_like2 + (t/16)*mp.power(log_sp2,2))   
    A_sum_exponent = s + (t/4) * log_sp1
    B_sum_exponent = 1 - s + (t/4) * log_sp2
    
    A_sum = 0.0
    for n in range(1,N+1): 
        A_sum += mp.power(n,(t/4)*mp.log(n) - A_sum_exponent)
    A = A0*A_sum
    
    B_sum = 0.0
    for m in range(1,M+1): 
        B_sum += mp.power(m,(t/4)*mp.log(m) - B_sum_exponent) 
    B = B0*B_sum
    
    H = A + B
    if z.imag==0: return H.real
    else: return H

'''End H_t approx functional eqn (Ht_AFE) block'''

'''Miscellanous useful functions'''
def Nt(t, T):
    t, T = mp.mpf(t), mp.mpf(T)
    Tsmall = T/(4 * mp.pi())
    N0 = Tsmall * mp.log(Tsmall) - Tsmall
    extra = (t/16.0) * mp.log(Tsmall)
    return (N0 + extra).real

def expected_zero_gap(t, T):
    t, T = mp.mpf(t), mp.mpf(T)
    return (T - 0.9*T)/(Nt(t,T) - Nt(t,0.9*T))

def append_data(filename, rows):
    with open(filename, 'a') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(rows)

def sign_change(x, y):
    if x*y < 0: return 1
    else: return 0
