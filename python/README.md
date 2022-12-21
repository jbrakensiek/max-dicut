This python code is provided to assist constructing THRESH rounding schemes. 

**Running**

- Run `python search-dicut.py` to generate a THRESH scheme for MAX DI-CUT.
- Run `python search-2and.py` to generate a THRESH scheme for MAX 2-AND.

**Files**

- `dicut.py` -- implements routines for evaluating rounding functions on instances of MAX DI-CUT / MAX 2-AND. This includes computing an optimal distribution of rounding functions for a given list of instances as well as an optimal distribution of instances for a given list of rounding functions.

- `experiments.py` -- implements the core routines for constructing rounding functions, finding new rounding functions that are optimal for a distribution of instances, and evaluating the performance of a distribution of rounding functions over an entire set of functions

- `search-dicut.py` -- iteratively builds a sequence of rounding functions which (in the limit) approach the optimal THRESH distribution for MAX DI-CUT

- `search-2and.py` -- the same for MAX 2-AND

**Dependencies**

Install the following dependencies in `pip` to get the code to run:

- `numpy`, `scipy`, `matplotlib`, `cvxpy`, `dill`

**Disclaimer**

The code provided will *not* produce the distributions found in the files `candidate_dicut.txt` and `candidate_2and.txt`, as those were found using the provided code with ad-hoc adjustments. Contact the authors if you would like more help toward getting comparable results.