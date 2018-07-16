# dbn_upper_bound
Computational effort to upper bound the de Bruijn-Newman constant as part of a Polymath project
-----------------------------------------------------------------------------------------------

As you may know, Prof. Terence Tao launched and has been moderating a Polymath project to upper bound the de Bruijn-Newman constant (dBN constant for brevity). This involves both 1) an Analytic or theory part (hard math) to derive/refine several formulas and estimates, and 2) a computational part to check whether the de Bruijn family of functions H_t are zero free in the regions marked for numerical verification.

This repo is meant to facilitate the computational aspect.

For general direction, theory, blogposts, discussions on latest results, etc. please head to [Prof. Tao's blog](https://terrytao.wordpress.com/)

For the Polymath proposal, please [check this link](https://terrytao.wordpress.com/2018/01/24/polymath-proposal-upper-bounding-the-de-bruijn-newman-constant/)

For the wiki, comprehensive list of papers, permanent record of results, please head to the [Polymath webpage](http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant)

For an intuitive understanding on establishing dBN bounds, please check this <a href='https://github.com/km-git-acc/dbn_upper_bound/blob/master/De.Bruijn.Newman.Barrier.approach.visuals.V1.6.pdf'>visual guide</a> created by <a href='https://github.com/rudolph-git-acc'>Rudolph</a> 

Computational Libraries and Machine Requirements
--------------------------------------------------------------------------------------------
Most of the recent work has been done in the Pari/GP, Arb and Julia languages. For large scale runs, the Arb scripts are recommended, and for mathlike readability, the other two.

For Julia, please check the Julia fork https://github.com/km-git-acc/DBNUpperBound.jl, maintained by WilCrofter.

There was also a lot of work done in Python in the earlier phase of the project, which may be of interest.

Please refer to the README files in the respective folders on how to use the different scripts. A lot of interesting discussion has also taken place in the 'Issues' threads (which are being used as general discussion+brainstorming+issues threads), which may be worth reading for detailed understanding. 

The scripts can likely run on any machine, but better configs will certainly help


Results
---------------------------------------------------------------------------------------------
Results can be checked on this wiki page,
http://michaelnielsen.org/polymath1/index.php?title=Zero-free_regions

(Currently, a dBN bound of 0.22 has been achieved unconditionally, and several tighter bounds conditional on RH verified to appropriate heights also demonstrated)

A lot of the summarized output from the scripts is present in the Output folder, some in the python/research folder, and large files posted as links in comments on the Blog.

