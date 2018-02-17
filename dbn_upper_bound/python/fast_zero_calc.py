"""
This example file shows how to quickly detect sign changes
and find approximate roots for t>0 near large heights T
It is inspired from the Gram's law and one of Terry's
derivations which shows that for t>0 the zeros start getting more
regularly spaced
Caution: Near smaller heights, it will miss many of the zeroes
"""
from mpmath import mp
from dbn_upper_bound.python.mputility import \
    (Ht_AFE_ADJ_AB, expected_zero_gap, sign_change, append_data)

mp.pretty = True
mp.dps = 40
t = 0.4
Htrootfilename = "Htroots_"+str(t)+".csv"
htroots = []
rootcount = 0
known_root = 684130
midpoint_estimate_flag = 0
for i in range(1, 5000001):
    avggap = expected_zero_gap(t, known_root)
    interval_min, interval_max = known_root + avggap/2, known_root + 3*avggap/2
    interval_min_eval = Ht_AFE_ADJ_AB(interval_min, t).real
    interval_max_eval = Ht_AFE_ADJ_AB(interval_max, t).real
    root_check = sign_change(interval_min_eval, interval_max_eval)
    if root_check == 1:
        try:
            approx_root = mp.findroot(lambda y: Ht_AFE_ADJ_AB(y, t).real, [interval_min, interval_max], solver="ridder")
            midpoint_estimate_flag = 0
        except:
            approx_root = (interval_min + interval_max)/2
            midpoint_estimate_flag = 1
        print(approx_root)
        rootcount += 1
        htroots.append([t, rootcount, approx_root, interval_min, interval_max, midpoint_estimate_flag])
        known_root = approx_root
    else: known_root += avggap
    if rootcount % 100 == 0:
        append_data(Htrootfilename, htroots)
        htroots = []

append_data(Htrootfilename, htroots)
