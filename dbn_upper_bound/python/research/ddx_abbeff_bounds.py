from mpmath import mp
mp.dps=30
mp.pretty=True

def log(x): return mp.log(x)

def exp(x): return mp.exp(x)

def pi(): return mp.pi()

def sqrt(x): return mp.sqrt(x)

def power(a,x): return mp.power(a,x)

def mpf(x): return mp.mpf(x)

pi=pi()
I = 1j

def alpha1(s):
    return 0.5/s + 1/(s-1) + 0.5*log(0.5*s/pi)

def alpha1prime(s):
    return -0.5/s**2 - 1/(s-1)**2 + 0.5/s

def H01(s):
    return 0.5*s*(s-1)*power(pi, -s/2)*sqrt(2*pi)*exp((0.5*s- 0.5)*log(0.5*s)-0.5*s)

#N bound for (d/dx)|B_eff/B0_eff|
def ddx_bbeff_bound(N,y=0.4,t=0.4):
    y,t = mpf(y), mpf(t)
    xN = 4.0*pi*N**2 - pi*t/4
    s = 0.5*(1+y-I*xN)
    deriv_sum = 0.0
    for n in range(1,N+1): deriv_sum += abs(0.5*I*(1+0.5*t*alpha1prime(s))*log(n)/n**(s + 0.5*t*alpha1(s)-0.25*t*log(n)))
    return (N, deriv_sum)

#N bound for (d/dx)|A_eff/A0_eff|
def ddx_aaeff_bound(N,y=0.4,t=0.4):
    y,t = mpf(y), mpf(t)
    xN = 4.0*pi*N**2 - pi*t/4
    s = 0.5*(1-y+I*xN)
    deriv_sum = 0.0
    for n in range(1,N+1): deriv_sum += abs(-0.5*I*(1+0.5*t*alpha1prime(s))*log(n)/n**(s + 0.5*t*alpha1(s)-0.25*t*log(n)))
    return (N, deriv_sum)

#N bound for (d/dx)|(A_eff + B_eff)/B0_eff|
def ddx_abbeff_bound(N,y=0.4,t=0.4):
    y,t = mpf(y), mpf(t)
    xN = 4*pi*N**2 - pi*t/4.0
    xNp1 = 4*pi*(N+1)**2 - pi*t/4.0
    delta = pi*y/(2*(xN - 6 - (14 + 2*y)/pi)) + 2*y*(7+y)*log(abs(1+y+I*xNp1)/(4*pi))/(xN*xN)
    ddx_b_pre = 1 + 0.5*t/(xN-6)
    ddx_a_pre = exp(delta)/N**y
    ddxsum_b, ddxsum_a = 0.0, 0.0
    for n in range(1,N+1):
        bn = exp((t/4.0)*log(n)**2)
        expo_b = 0.5*(1+y) + 0.5*t*((3*y-1)/(xNp1**2+1) + log(N))
        expo_a = 0.5*(1-y) + 0.5*t*((2-3*y)/xN**2 + log(N))
        b_numerator = bn*log(n)/2.0
        a_numerator = bn*(0.25*t*log(n)/(xN-6) + 0.25*log(abs(1-y+I*xNp1)*abs(1+y-I*xNp1)/(4*pi*n**2)) + 0.25*(3*t+1)/(xN*(xN-6)))
        ddxsum_b += b_numerator/n**expo_b  
        ddxsum_a += a_numerator/n**expo_a
    ddxsum = ddx_b_pre*ddxsum_b + ddx_a_pre*ddxsum_a
    return (N,ddxsum)
