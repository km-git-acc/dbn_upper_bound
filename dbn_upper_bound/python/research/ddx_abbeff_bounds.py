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

#N dependent upper bound for |(d/dx)(A_eff + B_eff)/B0_eff|
def ddx_abbeff_bound(N,y=0.4,t=0.4):
    y,t = mpf(y), mpf(t)
    xN = 4*pi*N**2 - pi*t/4.0
    xNp1 = 4*pi*(N+1)**2 - pi*t/4.0
    delta = pi*y/(2*(xN - 6 - (14 + 2*y)/pi)) + 2*y*(7+y)*log(abs(1+y+I*xNp1)/(4*pi))/(xN*xN)
    ddx_b_pre = 1 + 0.5*t/(xN-6)
    ddx_a_pre = exp(delta)/N**y
    ddxsum_b, ddxsum_a = 0.0, 0.0
    #print(ddxsum_b)
    for n in range(1,N+1):
        #print(ddxsum_b,type(N))
        bn = exp((t/4.0)*log(n)**2)
        expo_b = 0.5*(1+y) + 0.5*t*((3*y-1)/(xNp1**2+1) + log(N))
        expo_a = 0.5*(1-y) + 0.5*t*((2-3*y)/xN**2 + log(N))
        b_numerator = bn*log(n)/2.0
        a_numerator = bn*(0.25*t*log(n)/(xN-6) + 0.25*log(abs(1-y+I*xNp1)*abs(1+y-I*xNp1)/(16*(pi**2)*(n**2))) + 0.25*(3*t+1)/(xN*(xN-6)))
        ddxsum_b += b_numerator/n**expo_b  
        ddxsum_a += a_numerator/n**expo_a
    ddxsum = ddx_b_pre*ddxsum_b + ddx_a_pre*ddxsum_a
    return (N,ddxsum)

#mod of exact derivative at x -- |(d/dx)(A_eff + B_eff)/B0_eff|
def ddx_abbeff_x(x,y=0.4,t=0.4):
    y,t = mpf(y), mpf(t)
    N = int(sqrt((x+pi*t/4)/(4*pi)))
    ddxsum_b, ddxsum_a = 0.0, 0.0
    sb, sa = 0.5*(1+y-I*x), 0.5*(1-y+I*x) 
    ddx_b_pre = 1.0
    ddx_a_pre = exp((t/4.0)*(alpha1(sa)**2 - alpha1(sb)**2))*H01(sa)/H01(sb) 
    expo_b = sb + 0.5*t*alpha1(sb)
    expo_a = sa + 0.5*t*alpha1(sa)
    for n in range(1,N+1):
        bn = exp((t/4.0)*log(n)**2)
        b_numerator = -bn*log(n)*(-0.5*I - 0.25*I*t*alpha1prime(sb))
        a_numerator1 = -bn*log(n)*(0.5*I + 0.25*I*t*alpha1prime(sa))
        ddxloglambda = (-I*t/4.0)*(alpha1(sb)*alpha1prime(sb) + alpha1(1-sb)*alpha1prime(1-sb)) - (I/4.0)*log(sb/(2*pi)) - (I/4.0)*log((1-sb)/(2*pi)) + (I/4.0)*(1/sb + 1/(1-sb))
        a_numerator2 = -bn*ddxloglambda
        ddxsum_b += b_numerator/n**expo_b  
        ddxsum_a += (a_numerator1 + a_numerator2)/n**expo_a
    ddxsum = ddx_b_pre*ddxsum_b + ddx_a_pre*ddxsum_a
    return (x,abs(ddxsum))

#complex valued (A_eff+B_eff)/B0_eff
def abbeff_x(x,y=0.4,t=0.4):    
    s1 = 0.5*(1 - y + I*x)
    s2 = 0.5*(1 + y + I*x)
    N = int(sqrt((x+pi*t/4)/(4*pi)))
    alph1 = alpha1(s1)
    alph2 = alpha1(s2).conjugate()
    A0_expo = (t/4.0)*alph1*alph1
    B0_expo = (t/4.0)*alph2*alph2
    H01_est1 = H01(s1)
    H01_est2 = H01(s2).conjugate()
    
    A0 = exp(A0_expo)*H01_est1
    B0 = exp(B0_expo)*H01_est2
    A_sum = 0.0
    B_sum = 0.0
    for n in range(1, N+1):
        A_sum += 1/power(n, s1 + (t/2.0)*alph1 - (t/4.0)*log(n))
        B_sum += 1/power(n, 1-s1 + (t/2.0)*alph2 - (t/4.0)*log(n))
    A = A0 * A_sum
    B = B0 * B_sum
    return (x,(A+B)/B0)

#newton quotient of (A_eff+B_eff)/B0_eff for double checking the derivative
def newton_quot_abbeff_x(x,y=0.4,t=0.4,h=0.00001):
    newton_quot = (abbeff_x(x+h,y=0.4,t=0.4)[1] - abbeff_x(x,y=0.4,t=0.4)[1])/h
    return (x,abs(newton_quot))

