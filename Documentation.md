Our code builds on the arblib library which provides real and complex interval arithmetic:

https://arblib.org/

The `src` directory contains C++ wrappers for arb functions and utility functions related to the computation of various quantities in the experiment.

The `experiment` directory contains the main file `dicut.cpp` for rigorous verification of the approximation ratio of rounding schemes. 

```c++
class Linear_scheme 
```
Contains the definition of a rounding scheme which is a weighted combination of threshold functions, each being a linear function in the input bias. We don't force the total weight to be 1. 

```c++
Arb prob_by_linear_scheme_from_rel(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s)
```
Computes the probabilty that the given rounding scheme satisfies the configuration `(b1, b2, rho)`.

```c++
Arb obj(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound) 
```
Computes the objective value, which is `soundness(b1, b2, rho) - bound * completeness(b1, b2, rho)`. Since the weights of `Linear_scheme` isn't normalized to be 1, we need to multiply the completeness by the total weight in `s`. By certifying the nonnegativity of this objective value, we show that the rounding scheme `s` achieves an approximation ratio of `bound` in the input region.

```c++
Arb obj_d_b1(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound)
Arb obj_d_b2(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound)
Arb obj_d_rho(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound)
```
Computes the partial derivatives of `obj` with respect to `b1`, `b2`, and `rho` respectively.

```c++
int check_linear_scheme_from_rel(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound) 
```
Computes the `obj` value in the input region, and returns 1 if `obj` is nonnegative, -1 if negative, and 0 if inconclusive (i.e., the computed interval for `obj` contains both negative and nonnegative values).

```c++
int check(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound) 
```
The main function that checks the region defined by `(b1, b2, rho)`, with `s` being the rounding scheme and `bound` the bound to certify.

```c++
int main(int argc, char* argv[]) 
```
Reads the input and starts the checking on each region, whose endpoints are defined by the breakpoints of the rounding scheme. If the check fails, then the assertion `assert(good)` near the end of the function will fail; otherwise, the function exits normally.



The `python` directory has code that was used to discover the rounding schemes. It is independent of interval arithmetic code.
