from itertools import chain, combinations
from operator import mul
from functools import reduce
import csv
from mpmath import mp
mp.dps=30
mp.pretty=True

def read_data(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        data = list(reader)    
    return data

def append_data(filename, rows):
    with open(filename, 'a') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(rows)

def periodic_saving(i,N,filename,rows): 
    if(i%N==0): 
       append_data(filename,rows)
       rows[:] = []    
  
def mpf(x): return mp.mpf(x);

def mpc(x,y): return mp.mpc(x,y);

def stringtonum(z): return mp.mpmathify(z.replace('*I','j'));

def floor(n): return mp.floor(n);

def log(n): return mp.log(n);

def exp(n): return mp.exp(n);

def sqrt(x): return mp.sqrt(x);

def power(a,x): return mp.power(a,x);

def conj(a): return a.conjugate();

def gamma(z): return mp.gamma(z);

def zeta(z): return mp.zeta(z);

def cos(z): return mp.cos(z);

def sin(z): return mp.sin(z);

def sum(n,N,summand): return mp.nsum(summand,[n,N]);

def intnum(u_lowlim,u_uplim,integrand): return mp.quad(integrand,[u_lowlim,u_uplim]);

Pi=mp.pi();

I = 1j;

def alpha1(s): return(1/(2*s) + 1/(s-1) + (1/2)*log(s/(2*Pi)));

def alpha1prime(s): return(-1/(2*s**2) - 1/(s-1)**2 + 1/(2*s));

def H01(s): return((1/2)*s*(s-1)*Pi**(-s/2)*sqrt(2*Pi)*exp((s/2-1/2)*log(s/2)-s/2));

def S(N,sigmavar,t): return(sum(1,N,lambda n: n**(sigmavar + (t/4.0)*log(N**2/n))));

def C0(p): return((exp(Pi*I*(p**2/2 + 3/8)) - I*sqrt(2)*cos(Pi*p/2))/(2*cos(Pi*p)));

def B0_eff(x,y=0.4,t=0.4): return((1/8)*exp((t/4)*alpha1((1+y-I*x)/2)**2)*H01((1+y-I*x)/2));

def habbbound(N,y=0.4,t=0.4):
    xN = 4*Pi*N*N - t*Pi/4.0
    xNp1 = 4*Pi*(N+1)*(N+1) - t*Pi/4.0
    T1 = xN/2.0
    T2 = xNp1/2.0
    Tdash1 = T1 + t*Pi/8.0
    Tdash2 = T2 + t*Pi/8.0
    sNminus = 0.5*(1-y + 0.5j*(xN+xNp1))
    sNplus = 0.5*(1+y + 0.5j*(xN+xNp1))
    alph_sNminus, alph_sNplus = alpha1(sNminus), alpha1(sNplus)
    K = (xNp1 - xN)/(4.0*(xN-6))
    delta = Pi*y/(2*(xN - 6 - (14 + 2*y)/Pi)) + 2*y*(7+y)*log(abs(1+y+1j*xNp1)/(4*Pi))/(xN*xN)
    lambdafac = exp(delta)/power(N,y)
    a0 = sqrt(Tdash1/(2*Pi))
    
    epsdash1, epsdash2 = 0.0, 0.0
    for n in range(1,N+1):
        nf = mpf(n)
        modterm_sNminus, modterm_sNplus = abs(alph_sNminus - log(n)), abs(alph_sNplus - log(n)) 
        cnminus = (t*t/4.0)*(modterm_sNminus**2 + 2*K*modterm_sNminus + K*K) + (1/3.0) + t
        cnplus = (t*t/4.0)*(modterm_sNplus**2 + 2*K*modterm_sNplus + K*K) + (1/3.0) + t
        epsdash1+= cnminus*exp(cnminus/(2*(T1-3.33)))/power(nf,0.5*(1-y) + 0.5*t*(alph_sNminus.real - K) - 0.25*t*log(n))
        epsdash2+= cnplus*exp(cnplus/(2*(T1-3.33)))/power(nf,0.5*(1+y) + 0.5*t*(alph_sNplus.real - K) - 0.25*t*log(n))
    normalized_E1 = lambdafac*0.5*epsdash1/(T1-3.33)
    normalized_E2 = 0.5*epsdash2/(T1-3.33)
    
    E3_by_E3_main_decay = (1/8.0)*sqrt(Pi) * exp(-1*(t/64.0)*(Pi**2))*exp(0.181/(Tdash1 - 3.33))*(1+5.15/(a0-1.25)) 
    E3_main_decay_term = power(Tdash1, 1.5) * exp(-1 * Pi * T1/4.0)
    s2 = 0.5*(1+y) + 1j*T1
    alph2 = alpha1(s2).conjugate()
    B0_expo = (t/4.0)*alph2*alph2
    H01_est2 = H01(s2).conjugate()
    B0_eff = (1/8.0)*exp(B0_expo) * H01_est2
    E3_main_decay_by_B0eff = abs(E3_main_decay_term/B0_eff)
    normalized_E3 = E3_by_E3_main_decay*E3_main_decay_by_B0eff
    normalized_E = normalized_E1 + normalized_E2 + normalized_E3
    
    return [N, normalized_E1, normalized_E2,E3_by_E3_main_decay,E3_main_decay_by_B0eff,normalized_E3,normalized_E]


def habc_sharperbound(N,y=0.4,t=0.4):
    xN = 4*Pi*N*N - t*Pi/4.0
    T = xN/2.0
    Tdash = T + t*Pi/8.0
    Tdash_small = Tdash/(2*Pi)
    a = sqrt(Tdash_small)
        
    Res = (1+y)/2.0 + (t/2.0)*log(xN/(4*Pi)) - t*max(0,1-3*y)/(2*xN**2)
    bt = lambda n: exp((t/4.0)*log(n)**2)
    cmn_summand = lambda n: (bt(n)/n**Res)*(exp(((t*t/32)*log(xN/(4*Pi*n*n))**2 + 0.313)/(T-3.33))-1) 
    e1 = exp(0.02*y)*((xN/(4*Pi))**(-y/2.0))*sum(1,N,lambda n: (n**(y+t*y/(2*(xN-6))))*cmn_summand(n))
    e2 = sum(1,N,lambda n: cmn_summand(n))
    
    e3_term1 = Tdash_small**(-1*(1+y)/4.0)
    e3_term2_expo1 = (-t/16.0)*log(xN/(4*Pi))**2
    e3_term2_expo2 = (3*abs(log(xN/(4*Pi) + I*Pi/2.0)) + 3.58)/(xN - 8.52)
    e3_term3_abc = 1.24*(3**y + 1/3**y)/(a - 0.125) + 3.46/(Tdash - 3.33)
    e3_term3_ab  = 1 + e3_term3_abc  
    e3_abc = e3_term1*exp(e3_term2_expo1 + e3_term2_expo2)*e3_term3_abc
    e3_ab  = e3_term1*exp(e3_term2_expo1 + e3_term2_expo2)*e3_term3_ab
    
    e_total_ab  = e1 + e2 + e3_ab
    e_total_abc = e1 + e2 + e3_abc
    
    return [N, e1, e2, e3_ab, e3_abc, e_total_ab, e_total_abc]



'''for N in range(3,1001):
    print(habc_sharperbound(N))'''
