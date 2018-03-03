from itertools import chain, combinations
from operator import mul
from functools import reduce
import csv
from mpmath import mp
mp.dps=30
mp.pretty=True

def append_data(filename, rows):
    with open(filename, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(rows)

def findsubsets(S,m):
    return set(combinations(S, m))

def powerset(iterable):
    xs = list(iterable)
    return list(chain.from_iterable(combinations(xs,n) for n in range(len(xs)+1)))

def deltaN(n,N):
    if n<=N: return 1.0
    else: return 0.0

def divdelta(n,div):
    if n%div==0: return 1.0
    else: return 0.0

def alpha1(s):
    return 0.5/s + 1/(s - 1) + 0.5 * mp.log(s/(2 * mp.pi()))

def H01(s):
    return 0.5 * s * (s-1) * mp.power(mp.pi(), -0.5 * s) * mp.sqrt(2 * mp.pi()) * mp.exp((0.5 * s - 0.5) * mp.log(0.5 * s)-0.5 * s)

def habbbound(N,y=0.4,t=0.4):
    xN = 4*mp.pi()*N*N - t*mp.pi()/4.0
    xNp1 = 4*mp.pi()*(N+1)*(N+1) - t*mp.pi()/4.0
    T1 = xN/2
    T2 = xNp1/2
    Tdash1 = T1 + t*mp.pi()/8.0
    Tdash2 = T2 + t*mp.pi()/8.0
    sNminus = 0.5*(1-y + 0.5j*(xN+xNp1))
    sNplus = 0.5*(1+y + 0.5j*(xN+xNp1))
    alph_sNminus, alph_sNplus = alpha1(sNminus), alpha1(sNplus)
    K = (xNp1 - xN)/(4.0*(xN-6))
    delta = mp.pi()*y/(2*(xN - 6 - (14 + 2*y)/mp.pi())) + 2*y*(7+y)*mp.log(abs(1+y+1j*xNp1)/(4*mp.pi))/(xN*xN)
    lambdafac = mp.exp(delta)/mp.power(N,y)
    a0 = mp.sqrt(Tdash1/(2*mp.pi()))
    
    epsdash1, epsdash2 = 0.0, 0.0
    for n in range(1,N+1):
        nf = mp.mpf(n)
        modterm_sNminus, modterm_sNplus = abs(alph_sNminus - mp.log(n)), abs(alph_sNplus - mp.log(n)) 
        cnminus = (t*t/4.0)*(modterm_sNminus**2 + 2*K*modterm_sNminus + K*K) + (1/3.0) + t
        cnplus = (t*t/4.0)*(modterm_sNplus**2 + 2*K*modterm_sNplus + K*K) + (1/3.0) + t
        epsdash1+= cnminus*mp.exp(cnminus/(2*(T1-3.33)))/mp.power(nf,0.5*(1+y) + 0.5*t*(alph_sNminus.real - K) - 0.25*t*mp.log(n))
        epsdash2+= cnplus*mp.exp(cnplus/(2*(T1-3.33)))/mp.power(nf,0.5*(1+y) + 0.5*t*(alph_sNplus.real - K) - 0.25*t*mp.log(n))
    normalized_E1 = lambdafac*0.5*epsdash1/(T1-3.33)
    normalized_E2 = 0.5*epsdash2/(T1-3.33)
    
    E3_C0 = (1/8.0)*mp.sqrt(mp.pi()) * mp.exp(-1*(t/64.0)*(mp.pi()**2)) * mp.power(Tdash2, 1.5) * mp.exp(-1 * mp.pi() * T/4.0)
    E3_vwf_bound = mp.exp(0.181/(Tdash1 - 3.33))*(1+5.15/(a0-1.25))
    E3 = E3_C0*E3_vwf_bound
    s2 = 0.5*(1+y) + 1j*T1
    alph2 = alpha1(s2).conjugate()
    B0_expo = (t/4.0)*alph2*alph2
    H01_est2 = H01(s2).conjugate()
    B0_eff = (1/8.0)*mp.exp(B0_expo) * H01_est2
    normalized_E3 = E3/B0_eff

    normalized_E = normalized_E1+normalized_E2+normalized_E3
    return [N, normalized_E]
    
