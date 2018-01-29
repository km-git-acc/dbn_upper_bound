from scipy.optimize import fsolve
from first_1000_zeta import first_1000_zeta
from utility import Ht_real

import sys

orig_stdout = sys.stdout
f = open('out_0001_1000_01_49.txt', 'w')
sys.stdout = f


first_1000_H0_zero = [x*2 for x in first_1000_zeta]

# k is a stopping parameter (1<=k<=1000) which corresponds to the kth zero of
# the H_0 fn while l is the starting parameter. Default is 1, 1000.

l = 1

k = 1000

for j in range(0, 50):
    t = .01*j
    for i, nearby_root in enumerate(first_1000_H0_zero[1:k]):
        Ht_root = fsolve(Ht_real, nearby_root, args=(t,))[0]
        print(t, i, Ht_root)

sys.stdout = orig_stdout
f.close()
