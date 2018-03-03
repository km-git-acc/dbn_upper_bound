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

def alpha1(s):
    return 0.5/s + 1/(s - 1) + 0.5 * mp.log(s/(2 * mp.pi()))

def H01(s):
    return 0.5 * s * (s-1) * mp.power(mp.pi(), -0.5 * s) * mp.sqrt(2 * mp.pi()) * mp.exp((0.5 * s - 0.5) * mp.log(0.5 * s)-0.5 * s)

def habbbound(N,y=0.4,t=0.4):
    xN = 4*mp.pi()*N*N - t*mp.pi()/4.0
    xNp1 = 4*mp.pi()*(N+1)*(N+1) - t*mp.pi()/4.0
    T1 = xN/2.0
    T2 = xNp1/2.0
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
        epsdash1+= cnminus*mp.exp(cnminus/(2*(T1-3.33)))/mp.power(nf,0.5*(1-y) + 0.5*t*(alph_sNminus.real - K) - 0.25*t*mp.log(n))
        epsdash2+= cnplus*mp.exp(cnplus/(2*(T1-3.33)))/mp.power(nf,0.5*(1+y) + 0.5*t*(alph_sNplus.real - K) - 0.25*t*mp.log(n))
    normalized_E1 = lambdafac*0.5*epsdash1/(T1-3.33)
    normalized_E2 = 0.5*epsdash2/(T1-3.33)
    return [N, normalized_E1, normalized_E2]

'''for N in range(1,501):
    print(habbbound(N))'''

def normalized_E3_xbound(x,y=0.4,t=0.4):
    T = x/2.0
    Tdash = T + t*mp.pi()/8.0
    a0 = mp.sqrt(Tdash/(2*mp.pi()))
            
    E3_C0 = (1/8.0)*mp.sqrt(mp.pi()) * mp.exp(-1*(t/64.0)*(mp.pi()**2)) * mp.power(Tdash, 1.5) * mp.exp(-1 * mp.pi() * T/4.0)
    E3_vwf_bound = mp.exp(0.181/(Tdash - 3.33))*(1+5.15/(a0-1.25))
    E3 = E3_C0*E3_vwf_bound
    s2 = 0.5*(1+y) + 1j*T
    alph2 = alpha1(s2).conjugate()
    B0_expo = (t/4.0)*alph2*alph2
    H01_est2 = H01(s2).conjugate()
    B0_eff = (1/8.0)*mp.exp(B0_expo) * H01_est2
    normalized_E3 = abs(E3/B0_eff)
    return [x, normalized_E3]

'''import time
s=time.time()
mp.pretty=True
y=0.4
t=0.4
xstart = 1.0
xend = xstart + 6*(10**5)
mesh_size = 0.1
csvfile="normalized_e3_data.csv"
evaldata=[]
i=0
x=xstart
while x < xend: 
    zrow = [t,y]+normalized_E3_xbound(x)
    evaldata.append(zrow)
    #print(zrow)
    if i % 5000 == 0: append_data(csvfile, evaldata); print(i,x); evaldata = []
    x += mesh_size
    i += 1

e=time.time()
e-s'''


