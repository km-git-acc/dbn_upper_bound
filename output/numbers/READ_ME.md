These files were generated using an early version of the code used to find
zeroes of H_t. The naming convention is as follows: "output_a_b_c_d.txt" refers
to the output of dbn_upper_bound/python/H_t_probable_first_1000_zeroes.py 

with
a = starting index of zero of H_t
b = ending index of zero of H_t
c = start of t-range divided by 100
d = end of t-range divided by 100

The t-increment is in steps of 0.01

For example, out_01_1000_20_22.txt means the following. The file contains
the first 1000 zeroes of H_t for t in (0.20, 0.21). (The range() loop is not
inclusive of the d in computation).

Edit: Files added on/after 29 January 2018 should be calculated using the
version of the python code from at least 28 January 2018.
