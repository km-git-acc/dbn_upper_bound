#ABC and helper functions can be taken from a different file for now, eg. from Ht_small_x_fixed_mesh_verification.py
#This part to be fixed later

import matplotlib.pyplot as plt
import numpy as np

#H/B0 plot

C1 = [abceff_x(x,y=0.45)/B0_eff(x,y=0.45) for x in [a/100.0 for a in range(1,30000)]]
C2 = [abceff_x(300,y)/B0_eff(300,y) for y in [a/1000.0 for a in reversed(range(400,450))]]
C3 = [abceff_x(x,y=0.4)/B0_eff(x,y=0.4) for x in [a/100.0 for a in reversed(range(1,30000))]]
C4 = [abceff_x(0,y)/B0_eff(0,y) for y in [a/1000.0 for a in range(400,450)]]
C=C1+C2+C3+C4
X = [x.real for x in C]
Y = [x.imag for x in C]
maxint_XY = int(mp.ceil(max(max(X),max(Y))))
range_axis = [a/100.0 for a in range(-100*maxint_XY,100*maxint_XY)]
sz=0.01
plt.scatter(X,Y, color='red', s=sz)
plt.scatter([0 for i in range_axis],[i for i in range_axis], color='blue', s=sz)
plt.scatter([i for i in range_axis],[0 for i in range_axis], color='blue', s=sz)

plt.show()


#H'/H plot

LC1 = [newton_quot_abc(x,0.45)/abceff_x(x,0.45) for x in [a/100.0 for a in range(1,30000)]]
LC2 = [newton_quot_abc(300,y)/abceff_x(300,y) for y in [a/1000.0 for a in reversed(range(400,450))]]
LC3 = [newton_quot_abc(x,0.4)/abceff_x(x,0.4) for x in [a/100.0 for a in reversed(range(1,30000))]]
LC4 = [newton_quot_abc(0,y)/abceff_x(0,y) for y in [a/1000.0 for a in range(400,450)]]
LC=LC1+LC2+LC3+LC4
X = [x.real for x in LC]
Y = [x.imag for x in LC]
maxint_XY1 = int(mp.ceil(max(max(X),max(Y))))
maxint_XY2 = int(mp.ceil(max(max(-Xi for Xi in X),max(-Yi for Yi in Y))))
range_axis = [a/100.0 for a in range(-100*maxint_XY2,25*maxint_XY)]
sz=0.01
plt.scatter(X,Y, color='red', s=sz)
plt.scatter([0 for i in range_axis],[i for i in range_axis], color='blue', s=sz)
plt.scatter([i for i in range_axis],[0 for i in range_axis], color='blue', s=sz)

plt.show()
