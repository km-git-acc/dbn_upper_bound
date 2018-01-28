"""
This module contains various math utilities for the main project
"""

import scipy
from scipy.integrate import quad
from scipy.optimize import fsolve
from cmath import *
from dbn_upper_bound.python.constants import PI, PI_sq


def Nt(t, T):
    """
    Evaluates equation (4) in Terry's blog
    https://terrytao.wordpress.com/2018/01/27/polymath15-first-thread-computing-h_t-asymptotics-and-dynamics-of-zeroes/
    :param t: the "time" parameter
    :param T: height
    :return: right side of (4) in the blog link
    """
    T = float(T)
    T_by_4PI = T/(4*PI)

    return T_by_4PI*log(T_by_4PI) - T_by_4PI + (t/16.0)*log(T_by_4PI)


def phi_decay(u, n_max=100):
    """
    Add docstring
    :param u:
    :param n_max:
    :return:
    """
    running_sum = 0
    for n in range(1, n_max+1):
        term1 = 2*PI_sq*pow(n, 4)*exp(9*u) - 3*PI*pow(n, 2)*exp(5*u)
        term2 = exp(-1*PI*pow(n, 2)*exp(4*u))
        running_sum += term1*term2

    return running_sum


def Ht_complex_integrand(u, z, t):
    """
    Add docstring
    :param u:
    :param z:
    :param t:
    :return:
    """
    return exp(t*u*u)*phi_decay(u)*cos(z*u)


def Ht_complex(z, t):
    """
    Add docstring
    :param z:
    :param t:
    :return:
    """
    #may work well only for small to medium values of z
    def real_func(a, b, c): return scipy.real(Ht_complex_integrand(a, b, c))
    def imag_func(a, b, c): return scipy.imag(Ht_complex_integrand(a, b, c))
    real_part = quad(real_func, 0, 10, args=(z, t))
    imag_part = quad(imag_func, 0, 10, args=(z, t))
    return (real_part[0] + 1j*imag_part[0], real_part[1], imag_part[1])


def Ht_complex_root_finding_helper(z_as_array, t):
    """
    Add docstring
    :param z_as_array:
    :param t:
    :return:
    """
    z = float(z_as_array[0]) + 1j*float(z_as_array[1])
    Ht = Ht_complex(z, t)[0]
    return (scipy.real(Ht), scipy.imag(Ht))


def Ht_complex_root_finder(complex_guess, t):
    """
    Add docstring
    :param complex_guess:
    :param t:
    :return:
    """
    result = fsolve(Ht_complex_root_finding_helper,
                    [scipy.real(complex_guess), scipy.imag(complex_guess)],
                    args=(t,))
    return result[0]+1j*result[1]


def Ht_complex_zlarge(z, t):
    """
    Approx formula for Ht for large z values
    check last section of
    https://terrytao.wordpress.com/2018/01/27/polymath15-first-thread-computing-h_t-asymptotics-and-dynamics-of-zeroes/
    :param z:
    :param t:
    :return:
    """
    x = float(scipy.real(z))
    y = float(scipy.imag(z))
    t = float(t)
    x_by_4PI = x/(4*PI)
    B = ((PI*t/16)+(x/4))*log(x_by_4PI) - (x/4) + PI*(9+y)/8
    A1 = PI_sq*sqrt(PI/(2*1j*x))*pow(x_by_4PI, (9+y)/4)
    A2 = exp(-1*PI*x/8)
    A3 = exp((t/16)*pow(log(x_by_4PI), 2) - (t*PI_sq/64))
    Ht_without_PIby8_amplitude = A1*A3*exp(-1j*B)
    Ht = Ht_without_PIby8_amplitude*A2
    Ht_magnitude = abs(Ht)
    return (Ht, Ht_magnitude, Ht_without_PIby8_amplitude)


def Ht_real_integrand(u, z, t):
    """
    Add docstring
    :param u:
    :param z:
    :param t:
    :return:
    """
    if abs(scipy.imag(z)) > 0 or abs(scipy.imag(t)) > 0 \
            or abs(scipy.imag(u)) > 0:
        print("complex values not allowed for this function")
        return "error"

    u, z, t = scipy.real(u), scipy.real(z), scipy.real(t)
    return scipy.real(exp(t*u*u)*phi_decay(u)*cos(z*u))


def Ht_real(z, t):
    """
    Add docstring
    :param z:
    :param t:
    :return:
    """
    if abs(scipy.imag(z)) > 0 or abs(scipy.imag(t)) > 0:
        print("complex values not allowed for this function")
        return "error"
    z, t = scipy.real(z), scipy.real(t)
    # return quad(Ht_real_integrand, 0, np.inf, args=(z,t))
    # causing overflow errors so np.inf replaced with 10
    return quad(Ht_real_integrand, 0, 10, args=(z, t))

#check phi_decay values
'''print(phi_decay(0.001))
print(phi_decay(0.01))
print(phi_decay(0.1))
print(phi_decay(0.5))
print(phi_decay(1))'''
