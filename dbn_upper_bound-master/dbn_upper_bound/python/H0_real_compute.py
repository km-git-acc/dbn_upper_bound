"""
This module computes the real part of H0
"""

import dbn_upper_bound.python.utility as utils


#check values of H0_real
for z in range(10, 15):
    print(2*z, utils.Ht_real(2*z, 0))

print(28.06, utils.Ht_real(28.06, 0))
print(28.12, utils.Ht_real(28.12, 0))
print(28.18, utils.Ht_real(28.18, 0))
print(28.24, utils.Ht_real(28.24, 0))
print(28.26, utils.Ht_real(28.26, 0))
print(28.28, utils.Ht_real(28.28, 0))
print(28.4, utils.Ht_real(28.4, 0))
print(29, utils.Ht_real(29, 0))
