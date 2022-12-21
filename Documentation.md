Our code builds on the arblib library which provides real and complex interval arithmetic:

https://arblib.org/

The most common types arblib provides are `arb_t` and `acb_t` which store a real/complex interval, respectively.

Convention: In this document, we use `acb_t` variables *a1* and *a2* for the biases, *b12* for the pairwise bias, and *rho* for the relative pairwise bias. *comp* denotes the completeness of the configuration, and *sound* denotes the soundness of the rounding scheme on the configuration. The rounding functions are specified by 5 arrays of `acb_t` variables, *weight*, *c1*, *c2*, *d1*, and *d2*. *weight[i]* is the weight of the i-th rounding function, which has value $c_1[i]\times b + d_1[i]$ for the first variable and $c_2[i]\times b + d_2[i]$ for the second variable. *length* is the number of threshold functions in our rounding scheme.


### common.h

```c
int acb_is_nan(const acb_t z);
int arb_is_nan(const arb_t z);
```

Returns nonzero iff *z* equals NaN.

```c
void arb_union_d(arb_t z, double a, double b, slong prec);
void acb_union_d(acb_t z, double a, double b, slong prec);
```

Sets *z* to a ball containing both *a* and *b*, with *prec* bits of precision.

```c
void acb_norm_cdf(acb_t res, const acb_t z, slong prec);
```

Sets *res* to $\Phi(z)$, where $\Phi$ is the standard normal cumulative distribution function, with *prec* bits of precision.

```c
void acb_norm_2d_cdf_zero_mean(acb_t res, const acb_t rho, slong prec);
```

Sets *res* to $\Pr[X <= 0 \wedge Y <= 0]$ where $X$ and $Y$ are two standard normal random variables with correlation *rho*, with *prec* bits of precision.

```c
int acb_crop_d(acb_t res, const acb_t z, double lo, double hi, slong prec);
```

Sets *res* to the intersection of *z* and a ball containing both *lo* and *hi*, with *prec* bits of precision.

```c
void b12_from_rho(acb_t res, const acb_t a1, const acb_t a2, const acb_t rho, slong prec);
```

Sets *res* to the pairwise bias computed from biases *a1*, *a2*, and relative pairwise bias *rho*, with *prec* bits of precision.

```c
int thresh_eval_helper(acb_ptr res, const acb_t z, void * param, slong order, slong prec);
void thresh_eval_relative(acb_ptr res, const acb_t rho, const acb_t t1r, const acb_t t2r, slong prec);
```

The second function sets *res* to a lower bound of  *Pr[X <= t1r and Y >= t2r]* where $X$ and $Y$ are two standard normal random variables with correlation *rho*, with *prec* bits of precision. The first function is called by the second function during evaluation of integrals.

```c
int value(acb_ptr res, const acb_t a1, const acb_t a2, const acb_t b12, slong prec);
void value_from_rho(acb_t res, const acb_t a1, const acb_t a2, const acb_t rho, slong prec);
```

Sets *res* to *comp*, completeness of the instance, computed from *a1*, *a2*, and *b12*, or from *a1*, *a2*, and *rho*, with *prec* bits of precision.

### lower.h

```c
int tri_check(const acb_t a1, const acb_t a2, const acb_t b12, slong prec)
```

Checks if triangle inequalities are satisfied, with $prec$ bits of precision. Returns 0 if $(a1 + a2 + b12 + 1 < 0) \vee (-a1 - a2 + b12 + 1 < 0) \vee (-a1 + a2 - b12 + 1 < 0) \vee (a1 - a2 - b12 + 1 < 0)$, and 1 otherwise.

```c
int tri_check_rho(const acb_t a1, const acb_t a2, const acb_t rho, slong prec)
```

Checks if triangle inequalities are satisfied, with $prec$ bits of precision. Returns 0 if $(a1 + a2 + b12 + 1 < 0) \vee (-a1 - a2 + b12 + 1 < 0) \vee (-a1 + a2 - b12 + 1 < 0) \vee (a1 - a2 - b12 + 1 < 0)$, and 2 if $(a1 + a2 + b12 + 1 \geq 0) \wedge (-a1 - a2 + b12 + 1 \geq 0) \wedge (-a1 + a2 - b12 + 1 \geq 0) \wedge (a1 - a2 - b12 + 1 \geq 0)$ where $b12 = a1*a2 + rho\sqrt{(1-a1^2)(1-a2^2)}$, and 1 otherwise (note that an interval can be neither negative nor non-negative). 

```c
void partial_derivative_a1(acb_t res, arb_t bound, const ulong length, const acb_t a1, const acb_t a2, const acb_t rho, const acb_t* weights, const acb_t* c1, const acb_t* c2, const acb_t* d1, const acb_t* d2, slong prec);
void partial_derivative_a2(acb_t res, arb_t bound, const ulong length, const acb_t a1, const acb_t a2, const acb_t rho, const acb_t* weights, const acb_t* c1, const acb_t* c2, const acb_t* d1, const acb_t* d2, slong prec);
```

Sets *res* to the partial derivative of _sound - bound * comp_ with respect to *a1* or *a2*, with $prec$ bits of precision. 

```c
void partial_derivative_rho(acb_t res, arb_t bound, const ulong length, const acb_t a1, const acb_t a2, const acb_t rho, const acb_t* weights, const acb_t* c1, const acb_t* c2, const acb_t* d1, const acb_t* d2, slong prec);
```

Sets _res_ to _2pi * sqrt(1 - rho^2)_ times the partial derivative of _sound - bound * comp_ with respect to *rho*, with $prec$ bits of precision. 

```c
int rounding_linear_check_relative_pairwise(const arb_t bound, const ulong length, const acb_t a1, const acb_t a2, const acb_t rho, const acb_t* weights, const acb_t* c1, const acb_t* c2, const acb_t* d1, const acb_t* d2, slong prec);
```

Computes the value of _sound - bound * comp_ and compare it with 0. Returns 1 if it is nonnegative, -1 if it is negative, 2 if the given intervals *a1, a2, rho* do not satisfy all the triangle inequalities or have _comp_ less than COMP_LIMIT (defined in `common.h`), and returns 0 otherwise.

```c
int linear_rounding_certify_relative_pairwise(arb_t bound, const ulong length, double a1lo, double a1hi, double a2lo, double a2hi, double rholo, double rhohi, const acb_t* weights, const acb_t* c1, const acb_t* c2, const acb_t* d1, const acb_t* d2, slong prec);
```

The main verification function. It first constructs intervals _a1, a2, rho_ from the given lower and upper bounds _a1lo, a1hi, a2lo, a2hi, rholo, rhohi_. If these intervals specify a region entirely of valid configurations (satisfying the triangle inequalities and having completeness at least COMP_LIMIT), it may compute various partial derivatives and, in case of monotonicity, shrink some of these intervals to an endpoint accordingly.  This function then calls `rounding_linear_check_relative_pairwise()` function to evaluate _sound - bound * comp_. If the call returns non-zero value, then this function exits accordingly. Otherwise, it splits the longest interval among _a1, a2, rho_ into two equal-length sub-intervals and recursively checks each of the two regions.

```c
void lower_bound_dicut_engine(int ii, int jj, double dbound);
```

Reads the input, constructs the rounding functions and passes them to `linear_rounding_certify_relative_pairwise()`. Note that the rounding functions are piecewise linear, and _ii_ and _jj_ are the pieces for the first and second variables respectively. _dbound_ is the ratio to be certified.