"""
Details about the module goes here
"""
from mpmath import mp

mp.dps = 30


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

