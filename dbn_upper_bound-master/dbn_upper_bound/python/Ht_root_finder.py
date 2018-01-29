"""
Find new roots for H_t1 using known roots for H_t2, and a root finding algorithm like fsolve. t2 will often be zero since we know zeta roots. 
A similar thing has also been done in Ht_probable_roots using first 1000 zeta zeroes.
"""
from scipy.optimize import fsolve
from utility import Ht_real


H0_z1 = 2*14.1347251417346937904572519835624702707842
H0_z2 = 2*21.022039638771554992628479593896902777
H0_z3 = 2*25.010857580145688763213790992562821818659549
H0_z4 = 2*30.42487612585951321031189753058409132018
H0_z5 = 2*32.9350615877391896906623689640749034888127156
H0_z6 = 2*37.586178158825671257217763480705332821405597

#give a nearby H0 zero, not necessarily from the above
nearby_root = H0_z5

for i in range(501):
    t = float(i)/1000.0
    Ht_root = fsolve(Ht_real, nearby_root, args=(t,))[0]
    print(t, Ht_root)
    nearby_root = Ht_root
