# dbn_upper_bound
## Computational effort to upper bound the de Bruijn-Newman constant as part of a Polymath project
-----------------------------------------------------------------------------------------------

As you may know, Prof. Terence Tao will be launching and moderating a Polymath project to upper bound the de Bruijn-Newman constant (dBN constant for brevity). This will involve both 1) an Analytic or theory part (hard math) to derive/refine several formulas and estimates, and 2) a computational part to check whether the de Bruijn family of functions H_t are zero free in the regions marked for numerical verification.

This repo is meant to facilitate the computational aspect.

For general direction, theory, blogposts, discussions on latest results, etc. please head to [Prof. Tao's blog](https://terrytao.wordpress.com/)

For the Polymath proposal, please [check this link](https://terrytao.wordpress.com/2018/01/24/polymath-proposal-upper-bounding-the-de-bruijn-newman-constant/)

For the wiki, comprehensive list of papers, permanent record of results, please head to the [Polymath webpage](http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant)


## Computational Requirements
--------------------------------------------------------------------------------------------
1. [Python 3.6](https://docs.python.org/3/library/venv.html)
2. Packages in python_requirements.txt 

For Julia, please check the Julia fork https://github.com/km-git-acc/DBNUpperBound.jl, maintained by WilCrofter.

For other languages, please update accordingly.

The algorithms can likely run on any machine, but better configs will certainly help

## Steps to setup a Python 3.6 virtual environment
1. Install the virtual environment
```bash
python3.6 -m venv venv
```
2. Activate the virtual environment
```bash
source venv/bin/activate
```
3. Make sure you have the required packages inside the virtual environment. 
```bash
pip install -r python_requirements.txt
```

## Running H_t evaluations
--------------------------------------------------------------------------------------------
There are two main files in the python folder containing commonly used functions, utility.py and mputility.py. mputility uses the mpmath library which is preferred for H_t(z) for large z values.

Getting to workable code is quite simple. For eg.

from mputility import *


Ht_AFE_ABC(10000000.0,0.2)


which evaluates H_t using the approx functional eqn for z=10000000.0 and t=0.2

For a sample file showing how a large range of z values can be explored and roots found, please check sample_afe_abc_calc.py. 

## Results
---------------------------------------------------------------------------------------------
Please express your results here. A sample format can be your name or handle, t checked, T and epsilon verified, etc. Once there are multiple results, we can arrange it as a table. Also, once you are sure about your results, please update them or ask them to be updated on the Polymath page as well.

Do double check your results, and triple check especially if you find a non real zero for any t>0, since the implications are huge! 


## Other Sections - Please add as needed
---------------------------------------------------------------------------------------------
