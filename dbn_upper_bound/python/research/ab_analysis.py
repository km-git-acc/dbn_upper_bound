"""
Details about the module goes here
"""
import csv
from mpmath import mp

mp.dps = 30

def append_data(filename, rows):
    with open(filename, 'a') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(rows)


def AB_analysis(z, t):
    z, t = mp.mpc(z), mp.mpc(t)
    s = (1 + 1j * z.real - z.imag) / 2
    tau = mp.sqrt(s.imag / (2 * mp.pi()))
    N = int(tau)
    M = int(tau)

    A_pre = (1 / 16) * s * (s - 1) * mp.power(mp.pi(), -1 * s / 2) * mp.gamma(s / 2)
    A_sum = 0.0
    for n in range(1, N + 1):
        if t.real > 0:
            A_sum += mp.exp((t / 16) * mp.power(mp.log((s + 4) / (2 * mp.pi() * n * n)), 2)) / mp.power(n, s)
        else:
            A_sum += 1 / mp.power(n, s)
    A = A_pre * A_sum

    B_pre = (1 / 16) * s * (s - 1) * mp.power(mp.pi(), (s - 1) / 2) * mp.gamma((1 - s) / 2)
    B_sum = 0.0
    for m in range(1, M + 1):
        if t.real > 0:
            B_sum += mp.exp((t / 16) * mp.power(mp.log((5 - s) / (2 * mp.pi() * m * m)), 2)) / mp.power(m, 1 - s)
        else:
            B_sum += 1 / mp.power(m, 1 - s)
    B = B_pre * B_sum

    B0 = B_pre * mp.exp((t / 16) * mp.power(mp.log((5 - s) / (2 * mp.pi())), 2))

    AplusB = A + B
    ABB0 = AplusB / B0
    return (AplusB, B0, ABB0, abs(ABB0))


def adjusted_AB_analysis(z, t):
    z, t = mp.mpc(z), mp.mpc(t)
    s = (1 + 1j * z.real - z.imag) / 2
    tau = mp.sqrt(s.imag / (2 * mp.pi()))
    N = int(tau)
    M = int(tau)

    s_like1 = (s + 4) / 2
    s_like2 = (5 - s) / 2
    initial_term = (1 / 4) * mp.sqrt(2 * mp.pi())

    A_pre = initial_term * mp.power(mp.pi(), -1 * s / 2) * mp.exp((s_like1 - 0.5) * mp.log(s_like1) - s_like1)
    A_sum = 0.0
    for n in range(1, N + 1):
        if t.real > 0:
            A_sum += mp.exp((t / 16) * mp.power(mp.log(s_like1 / (mp.pi() * n * n)), 2)) / mp.power(n, s)
        else:
            A_sum += 1 / mp.power(n, s)
    A = A_pre * A_sum

    B_pre = initial_term * mp.power(mp.pi(), (s - 1) / 2) * mp.exp((s_like2 - 0.5) * mp.log(s_like2) - s_like2)
    B_sum = 0.0
    ddxB_sum = 0.0
    for m in range(1, M + 1):
        if t.real > 0:
            mth_term = mp.exp((t / 16) * mp.power(mp.log(s_like2 / (mp.pi() * m * m)), 2)) / mp.power(m, 1 - s)
        else:
            mth_term = 1 / mp.power(m, 1 - s)
        B_sum += mth_term
        ddxB_sum += mth_term * mp.log(m)
    B = B_pre * B_sum
    B_common_term = mp.exp((t / 16) * mp.power(mp.log(s_like2 / mp.pi()), 2))
    B0 = B_pre * B_common_term
    ddxBB0 = 0.5j * ddxB_sum / B_common_term

    AplusB = A + B
    ABB0 = AplusB / B0
    return (AplusB, B0, ABB0, abs(ABB0), ddxBB0, abs(ddxBB0))

def deltaN(n,N):
    if n<=N: return 1.0
    else: return 0.0

def divdelta(n,div):
    if n%div==0: return 1.0
    else: return 0.0

def xmultibound(x):
    N=int(mp.sqrt(0.25*x/mp.pi()))
    sum1_1, sum1_2, sum12_1, sum12_2, sum123_1, sum123_2, sum1235_1, sum1235_2 = [0.0 for i in range(8)]
    factor2 = 1 - 1/mp.power(2.0,0.7+0.1*mp.log(N*N/2.0))
    factor3 = 1 - 1/mp.power(3.0,0.7+0.1*mp.log(N*N/3.0))
    factor5 = 1 - 1/mp.power(5.0,0.7+0.1*mp.log(N*N/5.0))
    factorN = 1/mp.power(N,0.4)
    pow2, pow3, pow5, pow6, pow10, pow15, pow30 = [mp.power(j,0.4) for j in [2,3,5,6,10,15,30]]
    expo2, expo3, expo5 = 0.2*mp.log(2), 0.2*mp.log(3), 0.2*mp.log(5)
    expo6, expo10, expo15, expo30 =  expo2+expo3, expo2+expo5, expo3+expo5, expo2+expo3+expo5
    exp23, exp25, exp35 = mp.exp(0.2*mp.log(2)*mp.log(3)), mp.exp(0.2*mp.log(2)*mp.log(5)), mp.exp(0.2*mp.log(3)*mp.log(5))
    exp235 = exp23*exp25*exp35
    for n in range(1,30*N+1):
        L1 = deltaN(n,N)
        L2 = deltaN(n,2*N)*divdelta(n,2);   
        if L2>0: L2 /= mp.power((n/2.0),expo2)
        L3 = deltaN(n,3*N)*divdelta(n,3)
        if L3>0 : L3 /= mp.power((n/3.0),expo3)
        L5 = deltaN(n,5*N)*divdelta(n,5)
        if L5>0 : L5 /= mp.power((n/5.0),expo5)
        L6 = deltaN(n,6*N)*divdelta(n,6)
        if L6>0 : L6 /= (mp.power((n/6.0),expo6)*exp23) 
        L10 = deltaN(n,10*N)*divdelta(n,10)
        if L10>0 : L10 /= (mp.power((n/10.0),expo10)*exp25)
        L15 = deltaN(n,15*N)*divdelta(n,15)
        if L15>0 : L15 /= (mp.power((n/15.0),expo15)*exp35)
        L30 = deltaN(n,30*N)*divdelta(n,30)
        if L30>0 : L30 /= (mp.power((n/30.0),expo30)*exp235)
        R1 = L1
        R2 = L2/pow2
        R3 = L3/pow3
        R5 = L5/pow5
        R6 = L6/pow6 
        R10 = L10/pow10
        R15 = L15/pow15
        R30 = L30/pow30
        n=float(n)
        denom1 = mp.power(n,0.7+0.1*mp.log(N*N/n))
        denom2 = mp.power(n,0.3+0.1*mp.log(N*N/n))
        
        sum1_1+=abs(L1)/denom1; sum1_2+=abs(R1)/denom2
        sum12_1+=abs(L1-L2)/denom1; sum12_2+=abs(R1-R2)/denom2
        sum123_1+=abs(L1-L2-L3+L6)/denom1; sum123_2+=abs(R1-R2-R3+R6)/denom2
        sum1235_1+=abs(L1-L2-L3-L5+L6+L10+L15-L30)/denom1; sum1235_2+=abs(R1-R2-R3-R5+R6+R10+R15-R30)/denom2
    finalsum1 = sum1_1 - 1 + sum1_2*factorN
    finalsum12 = (sum12_1 - 1 + sum12_2*factorN)/factor2
    finalsum123 = (sum123_1 - 1 + sum123_2*factorN)/(factor2*factor3)
    finalsum1235 = (sum1235_1 - 1 + sum1235_2*factorN)/(factor2*factor3*factor5)
    return [finalsum1,finalsum12,finalsum123,finalsum1235]

'''mp.pretty=True
mltp = [0.5,0.8,1,1.5,2.5,3,3.5,4]
#mltp = [10,40]
xvals = [m*(10**6) for m in mltp]
for x in xvals:
    print([x]+xmultibound(x))'''

def condcache(N):
 pow2, pow3, pow5, pow6, pow10, pow15, pow30 = [mp.power(j,0.4) for j in [2,3,5,6,10,15,30]] 
 expo2, expo3, expo5 = 0.2*mp.log(2), 0.2*mp.log(3), 0.2*mp.log(5)
 expo6, expo10, expo15, expo30 =  expo2+expo3, expo2+expo5, expo3+expo5, expo2+expo3+expo5
 exp23, exp25, exp35 = mp.exp(0.2*mp.log(2)*mp.log(3)), mp.exp(0.2*mp.log(2)*mp.log(5)), mp.exp(0.2*mp.log(3)*mp.log(5))
 exp235 = exp23*exp25*exp35
 condc=[]
 for n in range(0,30*N+1):
     try:   
        L1 = deltaN(n,N)
        L2 = deltaN(n,2*N)*divdelta(n,2);   
        if L2>0: L2 /= mp.power((n/2.0),expo2)
        L3 = deltaN(n,3*N)*divdelta(n,3)
        if L3>0 : L3 /= mp.power((n/3.0),expo3)
        L5 = deltaN(n,5*N)*divdelta(n,5)
        if L5>0 : L5 /= mp.power((n/5.0),expo5)
        L6 = deltaN(n,6*N)*divdelta(n,6)
        if L6>0 : L6 /= (mp.power((n/6.0),expo6)*exp23) 
        L10 = deltaN(n,10*N)*divdelta(n,10)
        if L10>0 : L10 /= (mp.power((n/10.0),expo10)*exp25)
        L15 = deltaN(n,15*N)*divdelta(n,15)
        if L15>0 : L15 /= (mp.power((n/15.0),expo15)*exp35)
        L30 = deltaN(n,30*N)*divdelta(n,30)
        if L30>0 : L30 /= (mp.power((n/30.0),expo30)*exp235)
        R1 = L1
        R2 = L2/pow2
        R3 = L3/pow3
        R5 = L5/pow5
        R6 = L6/pow6 
        R10 = L10/pow10
        R15 = L15/pow15
        R30 = L30/pow30 
        condc.append([0.0, L1, R1, L1-L2, R1-R2, L1-L2-L3+L6, R1-R2-R3+R6, L1-L2-L3-L5+L6+L10+L15-L30, R1-R2-R3-R5+R6+R10+R15-R30])
     except:  condc.append([0.0 for _ in range(9)])
 return condc

def abtoybound(N,y,t,cond):
    sigma1, sigma2 = 0.5*(1+y), 0.5*(1-y)
    sum1_L, sum1_R, sum12_L, sum12_R, sum123_L, sum123_R, sum1235_L, sum1235_R = [0.0 for _ in range(8)]
    ddxsum1_L, ddxsum1_R, ddxsum12_L, ddxsum12_R, ddxsum123_L, ddxsum123_R, ddxsum1235_L, ddxsum1235_R = [0.0 for _ in range(8)]
    factorN = 1/mp.power(N,0.4)
    for n in range(1,30*N+1):
        nf=float(n)
        denom1 = mp.power(nf,sigma1+(t/4.0)*mp.log(N*N/nf))
        denom2 = mp.power(nf,sigma2+(t/4.0)*mp.log(N*N/nf))        
        term1_L = abs(cond[n][1]/denom1); term1_R = abs(cond[n][2]/denom2); sum1_L+=term1_L; sum1_R+=term1_R; ddxsum1_L+=mp.log(n)*term1_L; ddxsum1_R+=mp.log(n)*term1_R
        term12_L = abs(cond[n][3]/denom1); term12_R = abs(cond[n][4]/denom2); sum12_L+=term12_L; sum12_R+=term12_R; ddxsum12_L+=mp.log(n)*term12_L; ddxsum12_R+=mp.log(n)*term12_R
        term123_L = abs(cond[n][5]/denom1); term123_R = abs(cond[n][6]/denom2); sum123_L+=term123_L; sum123_R+=term123_R; ddxsum123_L+=mp.log(n)*term123_L; ddxsum123_R+=mp.log(n)*term123_R
        term1235_L = abs(cond[n][7]/denom1); term1235_R = abs(cond[n][8]/denom2); sum1235_L+=term1235_L; sum1235_R+=term1235_R; ddxsum1235_L+=mp.log(n)*term1235_L; ddxsum1235_R+=mp.log(n)*term1235_R
    sum1_L, sum12_L, sum123_L, sum1235_L = sum1_L - 1, sum12_L - 1, sum123_L - 1, sum1235_L - 1
    sum1_R, sum12_R, sum123_R, sum1235_R = sum1_R*factorN, sum12_R*factorN, sum123_R*factorN, sum1235_R*factorN
    ddxsum1_L, ddxsum12_L, ddxsum123_L, ddxsum1235_L = 0.5*ddxsum1_L, 0.5*ddxsum12_L, 0.5*ddxsum123_L, 0.5*ddxsum1235_L
    ddxsum1_R, ddxsum12_R, ddxsum123_R, ddxsum1235_R = 0.5*ddxsum1_R*factorN, 0.5*ddxsum12_R*factorN, 0.5*ddxsum123_R*factorN, 0.5*ddxsum1235_R*factorN
    
    abdiff1, ddxsum1       = 1 - sum1_L    - sum1_R,    ddxsum1_L    + ddxsum1_R
    abdiff12, ddxsum12     = 1 - sum12_L   - sum12_R,   ddxsum12_L   + ddxsum12_R
    abdiff123, ddxsum123   = 1 - sum123_L  - sum123_R,  ddxsum123_L  + ddxsum123_R
    abdiff1235, ddxsum1235 = 1 - sum1235_L - sum1235_R, ddxsum1235_L + ddxsum1235_R
    return [sum1_L,sum1_R,abdiff1,ddxsum1_L,ddxsum1_R,ddxsum1,sum12_L,sum12_R,abdiff12,ddxsum12_L,ddxsum12_R,ddxsum12,sum123_L,sum123_R,abdiff123,ddxsum123_L,ddxsum123_R,ddxsum123,sum1235_L,sum1235_R,abdiff1235,ddxsum1235_L,ddxsum1235_R,ddxsum1235]

def abtoy_sharperbound(N,y,t,cond):
    sigma1 = 0.5*(1+y)
    sum1, sum2, sum3, sum5 = [0.0 for _ in range(4)]
    b1 = 1
    a1 = mp.power(N,-0.4)
    for n in range(2,30*N+1):
        nf=float(n)
        denom = mp.power(nf,sigma1+(t/4.0)*mp.log(N*N))
        bn = mp.exp((t/4.0)*mp.power(mp.log(nf),2))*abs(cond[n][1])
        an = mp.power(nf/N,0.4)*bn        
        bn2, bn3, bn5 = bn*abs(cond[n][3]), bn*abs(cond[n][5]), bn*abs(cond[n][7])
        an2, an3, an5 = an*abs(cond[n][4]), an*abs(cond[n][6]), an*abs(cond[n][8])
        sum1 += max(abs(bn+an)/(1+a1), abs(bn-an)/(1-a1))/denom
        sum2 += max(abs(bn2+an2)/(1+a1), abs(bn2-an2)/(1-a1))/denom
        sum3 += max(abs(bn3+an3)/(1+a1), abs(bn3-an3)/(1-a1))/denom
        sum5 += max(abs(bn5+an5)/(1+a1), abs(bn5-an5)/(1-a1))/denom
    return [sum1,sum2,sum3,sum5]

'''import time
s=time.time()
mp.pretty=True
y=0.4
t=0.4
Nrange = range(1,2001)
Nrows = []
abtoydata="abtoy_sharperbound_data.csv"
c = condcache(N)
for N in Nrange:
    c = condcache(N)
    x_lower,x_upper = mp.ceil(4*mp.pi()*N**2), mp.floor(4*mp.pi()*(N+1)**2)
    Nth_row = [t,y,x_lower,x_upper,N]+abtoy_sharperbound(N,y,t,c)
    Nrows.append(Nth_row)
    print(Nth_row)
    if N % 10 == 0: append_data(abtoydata, Nrows); Nrows = []

append_data(abtoydata, Nrows)

e=time.time()
e-s'''

def abtoyx_e3(z,t):
    x=z.real
    xdash = x + mp.pi()*t/4.0
    y=z.imag
    sigma1, sigma2 = 0.5*(1+y), 0.5*(1-y)
    s1, s2 = sigma1 + 0.5j*xdash, sigma2 + 0.5j*xdash
    N=int(mp.sqrt(0.25*x/mp.pi()))
    sum1_L, sum1_R = 0.0, 0.0
    factor2 = 1 - 1/mp.power(2.0,s1+(t/4.0)*mp.log(N*N/2.0))
    factor3 = 1 - 1/mp.power(3.0,s1+(t/4.0)*mp.log(N*N/2.0))
    factorN = mp.power(N,-0.4)
    for n in range(1,N+1):
        n=float(n)
        sum1_L+=mp.power(n,-1*s1-(t/4.0)*mp.log(N*N/n)) 
        sum1_R+=mp.power(n,-1*s2-(t/4.0)*mp.log(N*N/n))
    sum1_R = sum1_R*factorN
    sum12_L, sum12_R = sum1_L*factor2, sum1_R*factor2
    sum123_L, sum123_R = sum12_L*factor3, sum12_R*factor3
    absum1_L, absum1_R, absum12_L, absum12_R, absum123_L, absum123_R = abs(sum1_L), abs(sum1_R), abs(sum12_L), abs(sum12_R), abs(sum123_L), abs(sum123_R) 
    abdiff1, abdiff12, abdiff123 = absum1_L - absum1_R, absum12_L - absum12_R, absum123_L - absum123_R
    return [sum1_L, sum1_R, absum1_L, absum1_R, abdiff1, sum12_L, sum12_R, absum12_L, absum12_R, abdiff12, sum123_L, sum123_R, absum123_L, absum123_R, abdiff123]


'''import time
s=time.time()
mp.pretty=True
y=0.4
t=0.4
xstart = 10**4
xend = xstart+15000
zstart = xstart + 1j*y
mesh_size = 0.05
z=zstart
abtoyx_data="abtoy_x_eval_10k_to_25k_mesh_size_0.05.csv"
evaldata=[]
i=0
while z.real < xend: 
    N = int(mp.sqrt(0.25*z.real/mp.pi()))
    zrow = [t,y,z.real,N]+abtoyx_e3(z,t)
    evaldata.append(zrow)
    #print(zrow)
    if i % 2000 == 0: append_data(abtoyx_data, evaldata); print(z.real); evaldata = []
    z += mesh_size
    i += 1

append_data(abtoyx_data, evaldata); print(z.real); evaldata = []
'''

