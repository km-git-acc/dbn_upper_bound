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
