"""
This module helps understand the following error
| H_t - (A + B - C)/B0 | < 0.4
| H_t - (A + B)/B0 | < 0.3
| C/B0 | < 0.1

for t = 0.4 and y = 0.4 range
and x (say) up to 10^8
"""

import mpmath as mp
import matplotlib.pyplot as plt
import numpy as np

from dbn_upper_bound.python.mputility \
    import (Ht_AFE_A, Ht_AFE_B, Ht_AFE_C, Ht_AFE_B0, Ht_complex)


mp.dps = 30

t = 0.4
y = 0.4
x_min = 10E3
x_max = 10E6
x = np.linspace(x_min, x_max)

data = {}
data['abc_b0'] = [[], [], []]
data['ab_b0'] = [[], [], []]
data['c_b0'] = [[], [], []]
data['ht_ab_b0'] = [[], [], []]

for xi in x:
    z = mp.mpc(xi, y)
    a = Ht_AFE_A(z, t)
    b = Ht_AFE_B(z, t)
    b0 = Ht_AFE_B0(z, t)
    c = Ht_AFE_C(z, t)
    abc_b0 = (a + b - c)/b0
    ab_b0 = (a + b) / b0
    c_b0 = c / b0
    ht_ab_b0 = (Ht_complex(z, t) - a - b)/b0

    data['abc_b0'][0].append(abc_b0.real)
    data['abc_b0'][1].append(abc_b0.imag)
    data['abc_b0'][2].append(mp.norm(abc_b0))

    data['ab_b0'][0].append(ab_b0.real)
    data['ab_b0'][1].append(ab_b0.imag)
    data['ab_b0'][2].append(mp.norm(ab_b0))

    data['ht_ab_b0'][0].append(ht_ab_b0.real)
    data['ht_ab_b0'][1].append(ht_ab_b0.imag)
    data['ht_ab_b0'][2].append(mp.norm(ht_ab_b0))

fig = plt.subplots(2, 3)
plt.subplot(231)
plt.title(' (A + B - C)/B0 ')
plt.plot(data['abc_b0'][0], data['abc_b0'][1], '-*')
plt.subplot(232)
plt.title(' (A + B)/B0 ')
plt.plot(data['ab_b0'][0], data['ab_b0'][1], '-*r')
plt.subplot(233)
plt.title(' (H_t - A - B)/B0 ')
plt.plot(data['ht_ab_b0'][0], data['ht_ab_b0'][1], '-*r')
# plt.subplot(233)
# plt.title(' C/B0 ')
# plt.plot(data['c_b0'][0], data['c_b0'][1], '-*r')

plt.subplot(234)
plt.title(' |(A + B - C)/B0| ')
plt.plot(data['abc_b0'][2], '-*')
plt.subplot(235)
plt.title(' |(A + B)/B0| ')
plt.plot(data['ab_b0'][2], '-*r')
plt.subplot(236)
plt.title(' |(H_t - A - B)/B0| ')
plt.plot(data['ht_ab_b0'][2], '-*r')
# plt.subplot(236)
# plt.title(' |C/B0| ')
# plt.plot(data['c_b0'][2], '-*r')

plt.suptitle('x_min: ' + '%.0e'  % x_min
             + ', x_max: ' + '%.0e'  % x_max
             + '\n' +
             't: ' + '%.1f' % t + ', y: ' + '%.1f' % y
              )

