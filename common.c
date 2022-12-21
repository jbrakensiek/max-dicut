/*
  Copyright (c) 2022 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include "common.h"

int acb_is_nan(const acb_t z) {
    arb_t t;
    arb_init(t);
    acb_get_real(t, z);
    return arf_is_nan(arb_midref(t));
}

int arb_is_nan(const arb_t z) {
    return arf_is_nan(arb_midref(z));
}

void arb_union_d(arb_t z, double a, double b, slong prec) {
    arb_t aa, bb;
    arb_init(aa); arb_init(bb);
    
    arb_set_d(aa, a);
    arb_set_d(bb, b);

    arb_union(z, aa, bb, prec);

    arb_clear(aa);
    arb_clear(bb);
}

void acb_union_d(acb_t z, double a, double b, slong prec) {
    acb_t aa, bb;
    acb_init(aa); acb_init(bb);
    
    acb_set_d(aa, a);
    acb_set_d(bb, b);

    acb_union(z, aa, bb, prec);

    acb_clear(aa);
    acb_clear(bb);
}

void acb_norm_cdf(acb_t res, const acb_t z, slong prec) {
    acb_t y;
    acb_init(y);

    acb_set_si(y, 2);
    acb_sqrt(y, y, prec);
    acb_mul_si(y, y, -1, prec);
    acb_div(y, z, y, prec);

    acb_hypgeom_erfc(res, y, prec);
    acb_div_si(res, res, 2, prec);

    acb_clear(y);
}

void acb_norm_2d_cdf_zero_mean(acb_t res, const acb_t rho, slong prec) {
    acb_t y, z;
    acb_init(y); acb_init(z);

    acb_asin(y, rho, prec);
    acb_const_pi(z, prec);
    acb_mul_si(z, z, 2, prec);
    acb_div(y, y, z, prec);

    acb_one(z);

    acb_div_si(z, z, 4, prec);

    acb_add(res, y, z, prec);

    acb_clear(y);
    acb_clear(z);
}

int acb_crop_d(acb_t res, const acb_t z, double lo, double hi, slong prec) {
    arb_t tempr;
    arb_init(tempr);
    assert(acb_is_real(z));
    acb_get_real(tempr, z);

    arb_t tempr2;
    arb_init(tempr2);
    arb_union_d(tempr2, lo, hi, prec);

    int val = arb_intersection(tempr, tempr, tempr2, prec);
    acb_set_arb(res, tempr);

    arb_clear(tempr);
    arb_clear(tempr2);

    return val;
}

config_t config_init(const ulong i1, const ulong i2, const acb_t b12)
{
    config_t m = (config_t) malloc(sizeof(config));

    m->i1 = i1;
    m->i2 = i2;

    acb_init(m->b12);
    acb_set(m->b12, b12);

    return m;
}

config_t config_init_d(const ulong i1, const ulong i2, const double b12)
{
    config_t m = (config_t) malloc(sizeof(config));
    
    m->i1 = i1;
    m->i2 = i2;

    acb_init(m->b12);
    acb_set_d(m->b12, b12);
    
    return m;
}

int config_clear(config_t c) {
    acb_clear(c -> b12);
    free(c);
    return 0;
}



void b12_from_rho(acb_t res, const acb_t a1, const acb_t a2, const acb_t rho, slong prec) {
    acb_t tmp1, tmp2, tmp3;
    
    acb_init(tmp1);
    acb_init(tmp2);
    acb_init(tmp3);
    
    // b12 = rho * sqrt(1 - a1^2) * sqrt(1 - a2^2) + a1 * a2;
    
    acb_zero(res);
    
    acb_one(tmp1);
    acb_submul(tmp1, a1, a1, prec);

    acb_sqrt(tmp1, tmp1, prec);    
    acb_set_arb(tmp1, acb_realref(tmp1));

    acb_one(tmp2);
    acb_submul(tmp2, a2, a2, prec);
    acb_sqrt(tmp2, tmp2, prec);
    acb_set_arb(tmp2, acb_realref(tmp2));

    acb_mul(tmp3, tmp1, tmp2, prec);
    acb_mul(res, rho, tmp3, prec);   
    acb_addmul(res, a1, a2, prec);
    
    acb_clear(tmp1);
    acb_clear(tmp2);
    acb_clear(tmp3);    
}

int thresh_eval_helper(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();
        
    acb_t t1r, t2r;
    acb_init(t1r), acb_init(t2r);
    acb_set(t1r, ((acb_t*)(param))[0]);
    acb_set(t2r, ((acb_t*)(param))[1]);

    // 1/sqrt(1-z^2) * exp(-(t1r^2 + 2*t1r*t2r*z + t2r^2) / (2 * (1 - z*z)))

    // a = 1 - z*z
    acb_t a;
    acb_init(a);
    acb_one(a);
    acb_submul(a, z, z, prec);

    //flint_printf("int_a: "); acb_printd(a, 10); flint_printf("\n");
    
    // b = t1r^2 + 2*t1r*t2r*z + t2r^2
    acb_t b;
    acb_init(b);
    acb_mul_si(b, t1r, 2, prec);
    acb_mul(b, b, t2r, prec);
    acb_mul(b, b, z, prec);
    acb_addmul(b, t1r, t1r, prec);
    acb_addmul(b, t2r, t2r, prec);

    //flint_printf("int_b: "); acb_printd(b, 10); flint_printf("\n");

    // c = exp(-b/2*a)
    acb_t c;
    acb_init(c);
    acb_mul_si(c, a, -2, prec);
    acb_div(c, b, c, prec);
    acb_exp(c, c, prec);

    //flint_printf("int_c: "); acb_printd(c, 10); flint_printf("\n");

    // finish
    acb_sqrt_analytic(a, a, order != 0, prec);
    acb_div(res, c, a, prec);

    //flint_printf("int_res: "); acb_printd(res, 10); flint_printf("\n");
    
    acb_clear(a);
    acb_clear(b);
    acb_clear(c);

    return 0;
}

void thresh_eval_relative(acb_ptr res, const acb_t rho, const acb_t t1r, const acb_t t2r, slong prec) { // only computes a lower bound

    acb_t one;
    acb_init(one);
    acb_one(one);
    
    if (acb_contains(rho, one)) { // compute one-dimensional integrals
        acb_clear(one);
        // When rho = 1, Pr_{rho}[X <= t1 && Y >= t2] = normalCDF(t1) - normalCDF(t2) if t2 <= t1 and 0 otherwise.
        // Assuming t1r and t2r are singleton sets (which holds if called by thresh_eval2_lower_relative_pairwise() )
        assert(acb_is_real(t1r));     
        assert(acb_is_real(t2r));
        if (arb_le(acb_realref(t1r), acb_realref(t2r))) {
            acb_zero(res);
            return;
        } else {
            acb_t val_1, val_2;
            acb_init(val_1), acb_init(val_2);
            acb_norm_cdf(val_1, t1r, prec);
            acb_norm_cdf(val_2, t2r, prec);
            acb_sub(res, val_1, val_2, prec);
            acb_clear(val_1), acb_clear(val_2);
            return;
        }
    }

    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

   // flint_printf("rho: "); acb_printd(rho, 10); flint_printf("\n");
   // flint_printf("t1r: "); acb_printd(t1r, 10); flint_printf("\n");
   // flint_printf("t2r: "); acb_printd(t2r, 10); flint_printf("\n");

    acb_t param[2];
    acb_init(param[0]);
    acb_set(param[0], t1r);
    acb_init(param[1]);
    acb_set(param[1], t2r);


    mag_t tol; mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -prec);

    slong goal = prec;

    acb_t a, b;
    acb_init(a);
    acb_init(b);
    acb_set_si(a, 0);
    acb_mul_si(b, rho, -1, prec);
    
    
    // compute a lower bound of b
    assert(acb_is_real(b));
    arf_t b_lowerb;
    arb_t b_tmp;
    arf_init(b_lowerb);
    arb_init(b_tmp);
    
    acb_get_real(b_tmp, b);
    arb_get_lbound_arf(b_lowerb, b_tmp, prec);
    arb_set_arf(b_tmp, b_lowerb);
    acb_set_arb(b, b_tmp);
    arb_clear(b_tmp);
    arf_clear(b_lowerb);
    
    

    acb_zero(res);
   // flint_printf("b: "); acb_printd(b, 10); flint_printf("\n");
    acb_calc_integrate(res, thresh_eval_helper, param, a, b, goal, tol, options, prec);
   // flint_printf("int_res: "); acb_printd(res, 10); flint_printf("\n");

    /// divide by 2*pi
    acb_t p;
    acb_init(p);
    acb_const_pi(p, prec);
    acb_mul_si(p, p, 2, prec);
    acb_div(res, res, p, prec);

    acb_t q, r;
    acb_init(q);
    acb_init(r);
    acb_mul_si(q, t1r, 1, prec);
    acb_mul_si(r, t2r, -1, prec);
    acb_norm_cdf(q, q, prec);
    acb_norm_cdf(r, r, prec);
    
    acb_addmul(res, q, r, prec);

   // flint_printf("int_res2: "); acb_printd(res, 10); flint_printf("\n");

    acb_clear(a);
    acb_clear(b);
    acb_clear(p);
    acb_clear(q);
    acb_clear(r);
    acb_clear(param[0]), acb_clear(param[1]);
}

// equivalent of value from dicut.py

int value(acb_ptr res, const acb_t a1, const acb_t a2, const acb_t b12, slong prec)
{
    acb_one(res);

    acb_addmul_si(res, a1, 1, prec);
    acb_addmul_si(res, a2, -1, prec);
    acb_addmul_si(res, b12, -1, prec);
    acb_div_si(res, res, 4, prec);

    return 0;
}



void value_from_rho(acb_t res, const acb_t a1, const acb_t a2, const acb_t rho, slong prec)
{
    acb_t b12;
    acb_init(b12);
    
    b12_from_rho(b12, a1, a2, rho, prec);
    
    acb_one(res);

    acb_addmul_si(res, a1, 1, prec);
    acb_addmul_si(res, a2, -1, prec);
    acb_addmul_si(res, b12, -1, prec);
    acb_div_si(res, res, 4, prec);
    
    acb_clear(b12);
}
