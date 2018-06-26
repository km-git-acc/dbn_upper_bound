
Running H_t evaluations
--------------------------------------------------------------------------------------------
There are two main files in the python folder containing commonly used functions, utility.py and mputility.py. mputility uses the mpmath library which is preferred for H_t(z) for large z values.

Getting to workable code is quite simple. For eg.

from mputility import *


Ht_AFE_ABC(10000000.0,0.2)


which evaluates H_t using the approx functional eqn for z=10000000.0 and t=0.2

For a sample file showing how a large range of z values can be explored and roots found, please check sample_afe_abc_calc.py. 
