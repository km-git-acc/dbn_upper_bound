from scipy.integrate import quad
import numpy as np
from utility import *

def H0_real_integrand(u, z):
     return phi_decay(u)*cos(z*u)

def H0_real(z):
     #return quad(H0_real_integrand, 0, np.inf, args=(z)) #causing overflow errors so np.inf replaced with 10
     return quad(H0_real_integrand, 0, 10, args=(z))

#check values of H0_real
for z in range(10,15):
      print 2*z, H0_real(2*z)

print 28.06,H0_real(28.06)
print 28.12,H0_real(28.12)
print 28.18,H0_real(28.18)
print 28.24,H0_real(28.24)
print 28.26,H0_real(28.26)
print 28.28,H0_real(28.28)
print 28.4,H0_real(28.4)
print 29,H0_real(29)

