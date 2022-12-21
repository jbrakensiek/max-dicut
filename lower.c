/*
  Copyright (c) 2022 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include "lower.h"

int tri_check(const acb_t a1, const acb_t a2, const acb_t b12, slong prec) {
    acb_t temp;
    acb_init(temp);

    // a1 + a2 + b12 + 1 >= 0
    
    acb_one(temp);
    acb_add(temp, temp, a1, prec);
    acb_add(temp, temp, a2, prec);
    acb_add(temp, temp, b12, prec);
    if (arb_is_negative(acb_realref(temp))) {
        acb_clear(temp);
        return 0;
    }

    // -a1 - a2 + b12 + 1 >= 0

    acb_one(temp);
    acb_sub(temp, temp, a1, prec);
    acb_sub(temp, temp, a2, prec);
    acb_add(temp, temp, b12, prec);
    if (arb_is_negative(acb_realref(temp))) {
        acb_clear(temp);
        return 0;
    }
    
    
    // -a1 + a2 - b12 + 1 >= 0

    acb_one(temp);
    acb_sub(temp, temp, a1, prec);
    acb_add(temp, temp, a2, prec);
    acb_sub(temp, temp, b12, prec);
    if (arb_is_negative(acb_realref(temp))) {
        acb_clear(temp);
        return 0;
    }

    // a1 - a2 - b12 + 1 >= 0

    acb_one(temp);
    acb_add(temp, temp, a1, prec);
    acb_sub(temp, temp, a2, prec);
    acb_sub(temp, temp, b12, prec);
    if (arb_is_negative(acb_realref(temp))) {
        acb_clear(temp);
        return 0;
    }

    acb_clear(temp);

    return 1;
}


int tri_check_rho(const acb_t a1, const acb_t a2, const acb_t rho, slong prec) {  
    
    // returns 0 if the entire region is invalid
    // returns 2 if the entire region is valid
    // returns 1 if there are both invalid and valid configurations in the region
    
    int all_valid = 1;    
    
    acb_t temp1, temp2, temp3;
    acb_init(temp1);
    acb_init(temp2);
    acb_init(temp3);
    
    // b12 = a1*a2 + rho * sqrt(1 - a_1^2) * sqrt(1 - a_2^2)

    // a1 + a2 + b12 + 1 >= 0
    // rho * sqrt( (1 - a_1) * (1 - a_2)) + sqrt( (1 + a_1) * (1 + a_2)) >= 0
    
    acb_one(temp1);
    acb_sub(temp1, temp1, a1, prec);
    acb_one(temp2);
    acb_sub(temp2, temp2, a2, prec);
    acb_mul(temp3, temp1, temp2, prec);
    acb_sqrt(temp3, temp3, prec);
    acb_mul(temp3, rho, temp3, prec);
    
    acb_one(temp1);
    acb_add(temp1, temp1, a1, prec);
    acb_one(temp2);
    acb_add(temp2, temp2, a2, prec);
    
    acb_mul(temp1, temp1, temp2, prec);
    acb_sqrt(temp1, temp1, prec);
    acb_add(temp3, temp3, temp1, prec);
    
    if (arb_is_negative(acb_realref(temp3))) {
        acb_clear(temp1);
        acb_clear(temp2);
        acb_clear(temp3);
        return 0;
    } else if (!arb_is_nonnegative(acb_realref(temp3))) {
        all_valid = 0;
    }

    // -a1 - a2 + b12 + 1 >= 0
    // rho * sqrt( (1 + a_1) * (1 + a_2)) + sqrt( (1 - a_1) * (1 - a_2)) >= 0
    
    acb_one(temp1);
    acb_add(temp1, temp1, a1, prec);
    acb_one(temp2);
    acb_add(temp2, temp2, a2, prec);
    acb_mul(temp3, temp1, temp2, prec);
    acb_sqrt(temp3, temp3, prec);
    acb_mul(temp3, rho, temp3, prec);
    
    acb_one(temp1);
    acb_sub(temp1, temp1, a1, prec);
    acb_one(temp2);
    acb_sub(temp2, temp2, a2, prec);
    
    acb_mul(temp1, temp1, temp2, prec);
    acb_sqrt(temp1, temp1, prec);
    acb_add(temp3, temp3, temp1, prec);
    
    if (arb_is_negative(acb_realref(temp3))) {
        acb_clear(temp1);
        acb_clear(temp2);
        acb_clear(temp3);
        return 0;
    } else if (!arb_is_nonnegative(acb_realref(temp3))) {
        all_valid = 0;
    }

    
    
    // -a1 + a2 - b12 + 1 >= 0
    // - rho * sqrt( (1 + a_1) * (1 - a_2)) + sqrt( (1 - a_1) * (1 + a_2)) >= 0
    
    acb_one(temp1);
    acb_add(temp1, temp1, a1, prec);
    acb_one(temp2);
    acb_sub(temp2, temp2, a2, prec);
    acb_mul(temp3, temp1, temp2, prec);
    acb_sqrt(temp3, temp3, prec);
    acb_mul(temp3, rho, temp3, prec);
    
    acb_one(temp1);
    acb_sub(temp1, temp1, a1, prec);
    acb_one(temp2);
    acb_add(temp2, temp2, a2, prec);
    
    acb_mul(temp1, temp1, temp2, prec);
    acb_sqrt(temp1, temp1, prec);
    acb_sub(temp3, temp1, temp3, prec);
    
    if (arb_is_negative(acb_realref(temp3))) {
        acb_clear(temp1);
        acb_clear(temp2);
        acb_clear(temp3);
        return 0;
    } else if (!arb_is_nonnegative(acb_realref(temp3))) {
        all_valid = 0;
    }


    // a1 - a2 - b12 + 1 >= 0
    // - rho * sqrt( (1 - a_1) * (1 + a_2)) + sqrt( (1 + a_1) * (1 - a_2)) >= 0
    
    acb_one(temp1);
    acb_sub(temp1, temp1, a1, prec);
    acb_one(temp2);
    acb_add(temp2, temp2, a2, prec);
    acb_mul(temp3, temp1, temp2, prec);
    acb_sqrt(temp3, temp3, prec);
    acb_mul(temp3, rho, temp3, prec);
    
    acb_one(temp1);
    acb_add(temp1, temp1, a1, prec);
    acb_one(temp2);
    acb_sub(temp2, temp2, a2, prec);
    
    acb_mul(temp1, temp1, temp2, prec);
    acb_sqrt(temp1, temp1, prec);
    acb_sub(temp3, temp1, temp3, prec);
    
    if (arb_is_negative(acb_realref(temp3))) {
        acb_clear(temp1);
        acb_clear(temp2);
        acb_clear(temp3);
        return 0;
    } else if (!arb_is_nonnegative(acb_realref(temp3))) {
        all_valid = 0;
    }


    acb_clear(temp1);
    acb_clear(temp2);
    acb_clear(temp3);

    if (all_valid)
        return 2;
    else
        return 1;
}



int rounding_linear_check_relative_pairwise(const arb_t bound, const ulong length,
                          const acb_t a1, const acb_t a2, const acb_t rho,
                          const acb_t* weights, const acb_t* c1, const acb_t* c2,
                          const acb_t* d1, const acb_t* d2, slong prec) {
    
    
    if (!tri_check_rho(a1, a2, rho, prec)) {
        return 2;
    }                      
    
    acb_t b12;
    acb_init(b12);
    
    b12_from_rho(b12, a1, a2, rho, prec);
    
    acb_t comp, ws;
    acb_init(comp);
    acb_init(ws);
    value(ws, a1, a2, b12, prec);
    acb_clear(b12);
    
    arb_t c_limit;
    arb_init(c_limit);
    arb_set_d(c_limit, COMP_LIMIT);
    if (arb_le(acb_realref(ws), c_limit)) { // making sure the completeness is at least COMP_LIMIT;
        arb_clear(c_limit);
        acb_clear(comp);
        acb_clear(ws);
        return 2;
    }
    arb_clear(c_limit);
    
    
    acb_set_arb(comp, bound);
    acb_mul(comp, comp, ws, prec);
    acb_zero(ws);

    acb_t res, res_lo, res_hi, temp, t1, t2;
    acb_init(res);
    acb_init(res_lo);
    acb_init(res_hi);
    acb_zero(res_lo);
    acb_zero(res_hi);
    acb_init(temp);
    acb_init(t1), acb_init(t2);
    
    arf_t t1_lb, t1_ub, t2_lb, t2_ub;
    arb_t t1_tmp, t2_tmp;
    
    arf_init(t1_lb);
    arf_init(t1_ub);
    arf_init(t2_lb);
    arf_init(t2_ub);
    arb_init(t1_tmp);
    arb_init(t2_tmp);

    for (int i = 0; i < length; i++) {
        acb_add(ws, ws, weights[i], prec);

        acb_zero(temp);


        acb_set(t1, d1[i]);
        acb_addmul(t1, a1, c1[i], prec);

        acb_set(t2, d2[i]);
        acb_addmul(t2, a2, c2[i], prec);
        // Need a lower bound for t1 and an upper bound for t2;
        
        assert(acb_is_real(t1));
        assert(acb_is_real(t2));
        
        acb_get_real(t1_tmp, t1);
        arb_get_lbound_arf(t1_lb, t1_tmp, prec);
        arb_set_arf(t1_tmp, t1_lb);
        acb_set_arb(t1, t1_tmp);
        
        acb_get_real(t2_tmp, t2);
        arb_get_ubound_arf(t2_ub, t2_tmp, prec);
        arb_set_arf(t2_tmp, t2_ub);
        acb_set_arb(t2, t2_tmp);

        thresh_eval_relative(temp, rho, t1, t2, prec);
        if (acb_is_nan(temp)) {
            assert(0);
            return 0;
        }

        acb_addmul(res_lo, weights[i], temp, prec);
        
        //if (VERBOSE)
            //acb_printd(temp, 10);                   
        
        acb_zero(temp);

        acb_set(t1, d1[i]);
        acb_addmul(t1, a1, c1[i], prec);

        acb_set(t2, d2[i]);
        acb_addmul(t2, a2, c2[i], prec);
        
        // Need an upper bound for t1 and a lower bound for t2;
        
        assert(acb_is_real(t1));
        assert(acb_is_real(t2));
        
        acb_get_real(t1_tmp, t1);
        arb_get_ubound_arf(t1_ub, t1_tmp, prec);
        arb_set_arf(t1_tmp, t1_ub);
        acb_set_arb(t1, t1_tmp);
        
        acb_get_real(t2_tmp, t2);
        arb_get_lbound_arf(t2_lb, t2_tmp, prec);
        arb_set_arf(t2_tmp, t2_lb);
        acb_set_arb(t2, t2_tmp);

        thresh_eval_relative(temp, rho, t1, t2, prec);
        if (acb_is_nan(temp)) {
            assert(0);
            return 0;
        }
            
        acb_addmul(res_hi, weights[i], temp, prec);       
        
        //if (VERBOSE)
        //    acb_printd(temp, 10);                   
    }
    
    acb_union(res, res_lo, res_hi, prec);
    
    arf_clear(t1_lb);
    arf_clear(t1_ub);
    arf_clear(t2_lb);
    arf_clear(t2_ub);
    arb_clear(t1_tmp);
    arb_clear(t2_tmp);
    

    acb_mul(comp, comp, ws, prec);

    //flint_printf("comp: "); acb_printd(comp, 10); flint_printf("\n");
    //flint_printf("res_lo: "); acb_printd(res_lo, 10); flint_printf("\n");
    //flint_printf("res_hi: "); acb_printd(res_hi, 10); flint_printf("\n");

    acb_sub(res, res, comp, prec);

    if (VERBOSE) acb_printd(res, 10);

    assert(acb_is_real(res));

    acb_clear(t1);
    acb_clear(t2);
    acb_clear(comp);
    acb_clear(ws);

    if (arb_is_nonnegative(acb_realref(res))) {
        acb_clear(res);
        return 1;
    }
    else if(arb_is_negative(acb_realref(res))) {
        acb_clear(res);
        return -1;
    }
    else {
        acb_clear(res);
        return 0;
    }
}


int STEPS=0;

void partial_derivative_a1(acb_t res, arb_t bound, const ulong length, const acb_t a1, const acb_t a2, const acb_t rho, const acb_t* weights,
                            const acb_t* c1, const acb_t* c2,
                            const acb_t* d1, const acb_t* d2, slong prec)
{
    acb_t t1, t2, temp1, temp2, temp3;
    acb_init(t1), acb_init(t2), acb_init(temp1), acb_init(temp2), acb_init(temp3);
    
    acb_zero(res);
    
    for (int i = 0; i < length; i++) {
    
        // soundness - bound * completeness
        // partial derivative = c1 * normalPDF(t1) * normalCDF((-t2 + rho * t1) / sqrt(1 - rho * rho)) - bound / 4 * (1 - b_2 + rho * sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1)
    
        acb_set(t1, d1[i]);
        acb_addmul(t1, a1, c1[i], prec);

        acb_set(t2, d2[i]);
        acb_addmul(t2, a2, c2[i], prec);
        

        acb_set(temp1, t2); // temp1 = t2
        acb_neg(temp1, temp1); // temp1 = -t2;
        acb_addmul(temp1, rho, t1, prec); // temp1 = -t2 + rho * t1;
        acb_one(temp2); // temp2 = 1;
        acb_submul(temp2, rho, rho, prec); // temp2 = 1 - rho * rho;
        acb_sqrt(temp2, temp2, prec); // temp2 = sqrt(1 - rho * rho);
        acb_div(temp1, temp1, temp2, prec); // temp1 = (-t2 + rho * t1) / sqrt(1 - rho * rho);
        acb_norm_cdf(temp1, temp1, prec); // temp1 = normalCDF((-t2 + rho * t1) / sqrt(1 - rho * rho));
        
        acb_mul(temp2, t1, t1, prec); // temp2 = t1 * t1;
        acb_div_si(temp2, temp2, 2, prec); // temp2 = t1 * t1 / 2;
        acb_neg(temp2, temp2); // temp2 = - t1 * t1 / 2;
        acb_exp(temp2, temp2, prec); // temp2 = exp(- t1 * t1) / 2;
        acb_const_pi(temp3, prec); // temp3 = pi;
        acb_mul_si(temp3, temp3, 2, prec); // temp3 = 2pi;
        acb_sqrt(temp3, temp3, prec); // temp3 = sqrt(2pi);
        acb_div(temp2, temp2, temp3, prec); // temp2 = normalPDF(t1);
        
        acb_mul(temp1, temp1, temp2, prec); // temp1 = normalPDF(t1) * normalCDF((-t2 + rho * t1) / sqrt(1 - rho * rho));
        acb_mul(temp1, c1[i], temp1, prec); // temp1 = c1 * normalPDF(t1) * normalCDF((-t2 + rho * t1) / sqrt(1 - rho * rho));
        
        acb_one(temp2); // temp2 = 1;
        acb_submul(temp2, a2, a2, prec); // temp2 = 1 - b2 * b2;
        acb_one(temp3);
        acb_submul(temp3, a1, a1, prec); // temp3 = 1 - b1 * b1;
        acb_div(temp2, temp2, temp3, prec); // temp2 = (1 - b2 * b2) / (1 - b1 * b1);
        acb_sqrt(temp2, temp2, prec); // temp2 = sqrt( (1 - b2 * b2) / (1 - b1 * b1) );
        acb_mul(temp2, temp2, a1, prec); // temp2 =  sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1;
        acb_mul(temp2, temp2, rho, prec); // temp2 =  rho * sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1;
        acb_one(temp3);
        acb_add(temp2, temp3, temp2, prec); // temp2 = 1 + rho * sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1;
        acb_sub(temp2, temp2, a2, prec); // temp2 = 1 - b2 + rho * sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1;
        
        acb_mul_arb(temp2, temp2, bound, prec); // temp2 = bound * (1 - b2 + rho * sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1);
        acb_div_si(temp2, temp2, 4, prec);  // temp2 = bound / 4 * (1 - b2 + rho * sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1);
        
        acb_sub(temp1, temp1, temp2, prec);
        
        //if (VERBOSE)
        //    acb_printd(temp1, 10);
        
        acb_addmul(res, weights[i], temp1, prec);
    }
    
    if (VERBOSE)
        acb_printd(res, 10);
    acb_clear(t1), acb_clear(t2), acb_clear(temp1), acb_clear(temp2), acb_clear(temp3);
}

void partial_derivative_a2(acb_t res, arb_t bound, const ulong length, const acb_t a1, const acb_t a2, const acb_t rho, const acb_t* weights,
                            const acb_t* c1, const acb_t* c2,
                            const acb_t* d1, const acb_t* d2, slong prec)
{
    acb_t t1, t2, temp1, temp2, temp3;
    acb_init(t1), acb_init(t2), acb_init(temp1), acb_init(temp2), acb_init(temp3);
    
    acb_zero(res);
    
    for (int i = 0; i < length; i++) {
    
        // soundness - bound * completeness
        // partial derivative = -c2 * normalPDF(-t2) * normalCDF((t1 - rho * t2) / sqrt(1 - rho * rho)) + bound / 4 * (1 + b_1 - rho * sqrt( (1 - b1 * b1) / (1 - b2 * b2) ) * b2)
    
        acb_set(t1, d1[i]);
        acb_addmul(t1, a1, c1[i], prec);

        acb_set(t2, d2[i]);
        acb_addmul(t2, a2, c2[i], prec);
        

        acb_set(temp1, t1); // temp1 = t1
        acb_submul(temp1, rho, t2, prec); // temp1 = t1 - rho * t2;
        acb_one(temp2); // temp2 = 1;
        acb_submul(temp2, rho, rho, prec); // temp2 = 1 - rho * rho;
        acb_sqrt(temp2, temp2, prec); // temp2 = sqrt(1 - rho * rho);
        acb_div(temp1, temp1, temp2, prec); // temp1 = (t1 - rho * t2) / sqrt(1 - rho * rho);
        acb_norm_cdf(temp1, temp1, prec); // temp1 = normalCDF((t1 - rho * t2) / sqrt(1 - rho * rho));
        
        acb_mul(temp2, t2, t2, prec); // temp2 = t2 * t2;
        acb_div_si(temp2, temp2, 2, prec); // temp2 = t2 * t2 / 2;
        acb_neg(temp2, temp2); // temp2 = - t2 * t2 / 2;
        acb_exp(temp2, temp2, prec); // temp2 = exp(- t2 * t2) / 2;
        acb_const_pi(temp3, prec); // temp3 = pi;
        acb_mul_si(temp3, temp3, 2, prec); // temp3 = 2pi;
        acb_sqrt(temp3, temp3, prec); // temp3 = sqrt(2pi);
        acb_div(temp2, temp2, temp3, prec); // temp2 = normalPDF(t2) = normalPDF(-t2);
        
        acb_mul(temp1, temp1, temp2, prec); // temp1 = normalPDF(-t2) * normalCDF((t1 - rho * t2) / sqrt(1 - rho * rho)) ;
        acb_mul(temp1, c2[i], temp1, prec); // temp1 = c2 * normalPDF(-t2) * normalCDF((t1 - rho * t2) / sqrt(1 - rho * rho)) ;
        acb_neg(temp1, temp1); // temp1 = - c2 * normalPDF(-t2) * normalCDF((t1 - rho * t2) / sqrt(1 - rho * rho)) ;
        
        acb_one(temp2); // temp2 = 1;
        acb_submul(temp2, a2, a2, prec); // temp2 = 1 - b2 * b2;
        acb_one(temp3);
        acb_submul(temp3, a1, a1, prec); // temp3 = 1 - b1 * b1;
        acb_div(temp2, temp3, temp2, prec); // temp2 = (1 - b1 * b1) / (1 - b2 * b2);
        acb_sqrt(temp2, temp2, prec); // temp2 = sqrt( (1 - b1 * b1) / (1 - b2 * b2) );
        acb_mul(temp2, temp2, a2, prec); // temp2 =  sqrt( (1 - b1 * b1) / (1 - b2 * b2) ) * b2;
        acb_mul(temp2, temp2, rho, prec); // temp2 =  rho * sqrt( (1 - b1 * b1) / (1 - b2 * b2) ) * b2;
        acb_one(temp3);
        acb_sub(temp2, temp3, temp2, prec); // temp2 = 1 - rho * sqrt( (1 - b1 * b1) / (1 - b2 * b2) ) * b2;
        acb_add(temp2, temp2, a1, prec); // temp2 = 1 + b1 - rho * sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1;
        
        acb_mul_arb(temp2, temp2, bound, prec); // temp2 = bound * (1 + b1 - rho * sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1);
        acb_div_si(temp2, temp2, 4, prec);  // temp2 = bound / 4 * (1 + b1 - rho * sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1);
        
        acb_add(temp1, temp1, temp2, prec);
        
        //if (VERBOSE)
        //    acb_printd(temp1, 10);
        
        acb_addmul(res, weights[i], temp1, prec);
    }
    
    if (VERBOSE)
        acb_printd(res, 10);
    acb_clear(t1), acb_clear(t2), acb_clear(temp1), acb_clear(temp2), acb_clear(temp3);
}


void partial_derivative_rho(acb_t res, arb_t bound, const ulong length, const acb_t a1, const acb_t a2, const acb_t rho, const acb_t* weights,
                            const acb_t* c1, const acb_t* c2,
                            const acb_t* d1, const acb_t* d2, slong prec)
{
    acb_t t1, t2, temp1, temp2, temp3;
    acb_init(t1), acb_init(t2), acb_init(temp1), acb_init(temp2), acb_init(temp3);
    
    acb_zero(res);
    
    for (int i = 0; i < length; i++) {
    
        // soundness - bound * completeness
        // partial derivative * 2pi * sqrt(1 - rho^2) = - exp( - (t_1^2 - 2*rho*t1*t2 + t2^2) / 2(1 - rho^2)) + pi/2 * sqrt(1 - rho^2) * bound * sqrt( (1 - a1^2) * (1 - a2^2)).
    
        acb_set(t1, d1[i]);
        acb_addmul(t1, a1, c1[i], prec);

        acb_set(t2, d2[i]);
        acb_addmul(t2, a2, c2[i], prec);
        
        acb_set_ui(temp1, 2);
        acb_mul(temp1, temp1, rho, prec);
        acb_mul(temp1, temp1, t1, prec);
        acb_mul(temp1, temp1, t2, prec);
        acb_neg(temp1, temp1);
        acb_addmul(temp1, t2, t2, prec);
        acb_addmul(temp1, t1, t1, prec);
        acb_neg(temp1, temp1);
        
        acb_one(temp2);
        acb_submul(temp2, rho, rho, prec);
        acb_mul_ui(temp2, temp2, 2, prec);
        acb_div(temp1, temp1, temp2, prec);
        acb_exp(temp1, temp1, prec);
        acb_neg(temp1, temp1); // - exp( - (t_1^2 - 2*rho*t1*t2 + t2^2) / 2(1 - rho^2))
        
        acb_one(temp2);
        acb_one(temp3);
        acb_submul(temp3, rho, rho, prec);
        acb_mul(temp2, temp2, temp3, prec);
        
        acb_one(temp3);
        acb_submul(temp3, a1, a1, prec);
        acb_mul(temp2, temp2, temp3, prec);
        
        acb_one(temp3);
        acb_submul(temp3, a2, a2, prec);
        acb_mul(temp2, temp2, temp3, prec);
        
        acb_sqrt(temp2, temp2, prec);
        
        acb_const_pi(temp3, prec);
        acb_mul(temp2, temp2, temp3, prec);
        acb_div_ui(temp2, temp2, 2, prec);
        acb_mul_arb(temp2, temp2, bound, prec); // pi/2 * sqrt(1 - rho^2) * bound * sqrt( (1 - a1^2) * (1 - a2^2))
        
        acb_add(temp1, temp1, temp2, prec);
        
        //if (VERBOSE)
        //    acb_printd(temp1, 10);
        
        acb_addmul(res, weights[i], temp1, prec);
    }
    
    if (VERBOSE)
        acb_printd(res, 10);
    acb_clear(t1), acb_clear(t2), acb_clear(temp1), acb_clear(temp2), acb_clear(temp3);
}

int linear_rounding_certify_relative_pairwise(arb_t bound, const ulong length, 
                            double a1lo, double a1hi, double a2lo, double a2hi,
                            double rholo, double rhohi, const acb_t* weights,
                            const acb_t* c1, const acb_t* c2,
                            const acb_t* d1, const acb_t* d2, slong prec) 
{
    
    
    if (VERBOSE)
        flint_printf("AT %08d: %.9f %.9f %.9f %.9f %.9f %.9f: ", STEPS++, a1lo, a1hi, a2lo, a2hi, rholo, rhohi);
    
    acb_t a1, a2, rho;
    acb_init(a1), acb_init(a2), acb_init(rho);


    acb_union_d(a1, a1lo, a1hi, prec);
    acb_union_d(a2, a2lo, a2hi, prec);
    acb_union_d(rho, rholo, rhohi, prec);
    
    int all_valid = (tri_check_rho(a1, a2, rho, prec) == 2);
    
    acb_t comp;
    acb_init(comp);
    value_from_rho(comp, a1, a2, rho, prec);
    
    arb_t c_limit;
    arb_init(c_limit);
    arb_set_d(c_limit, COMP_LIMIT);
    
    all_valid = all_valid && (arb_ge(acb_realref(comp), c_limit));
    
    arb_clear(c_limit);
    acb_clear(comp);
    
    if (all_valid) {
        if ((rholo != rhohi) && (rhohi - rholo <= DERIVATIVE_THRESHOLD_RHO)) {
        
            if (rholo != -1 && rhohi != 1) {
            

                if (VERBOSE)
                    flint_printf("   Testing der_rho: ");
                acb_t der;
                acb_init(der);
                partial_derivative_rho(der, bound, length, a1, a2, rho, weights, c1, c2, d1, d2, prec); // partial derivative of (s - r*c) w.r.t. rho
                if (arb_is_nonnegative(acb_realref(der))) {
                    acb_clear(der);
                    if (VERBOSE)
                        flint_printf("\n");
                    return linear_rounding_certify_relative_pairwise(bound, length, a1lo, a1hi, a2lo, a2hi, rholo, rholo,
                                        weights, c1, c2, d1, d2, prec);
                } else if (arb_is_nonpositive(acb_realref(der))) {
                    acb_clear(der);
                    if (VERBOSE)
                        flint_printf("\n");
                    return linear_rounding_certify_relative_pairwise(bound, length, a1lo, a1hi, a2lo, a2hi, rhohi, rhohi,
                                        weights, c1, c2, d1, d2, prec);
                } else {
                    if (VERBOSE)
                        flint_printf(" ");
                    acb_clear(der);
                } 
                
            }
        }
    
        if ((a1lo != a1hi) && (a1hi - a1lo <= DERIVATIVE_THRESHOLD_A)) {
            if (a1lo != -1 && a1hi != 1 && rholo != -1 && rhohi != 1) {
                if (VERBOSE)
                    flint_printf("   Testing der_b1: ");
                acb_t der;
                acb_init(der);
                partial_derivative_a1(der, bound, length, a1, a2, rho, weights, c1, c2, d1, d2, prec);
                
                if (arb_is_nonnegative(acb_realref(der))) {
                    acb_clear(der);
                    if (VERBOSE)
                        flint_printf("\n");
                    return linear_rounding_certify_relative_pairwise(bound, length, a1lo, a1lo, a2lo, a2hi, rholo, rhohi,
                                        weights, c1, c2, d1, d2, prec);
                } else if (arb_is_nonpositive(acb_realref(der))) {
                    acb_clear(der);
                    if (VERBOSE)
                        flint_printf("\n");
                    return linear_rounding_certify_relative_pairwise(bound, length, a1hi, a1hi, a2lo, a2hi, rholo, rhohi,
                                        weights, c1, c2, d1, d2, prec);
                } else {
                    if (VERBOSE)
                        flint_printf(" ");
                    acb_clear(der);
                } 
            }
        }
            
        if ((a2lo != a2hi) && (a2hi - a2lo <= DERIVATIVE_THRESHOLD_A)) {
            if (a2lo != -1 && a2hi != 1 && rholo != -1 && rhohi != 1) {
                if (VERBOSE)
                    flint_printf("   Testing der_b2: ");
                acb_t der;
                acb_init(der);
                partial_derivative_a2(der, bound, length, a1, a2, rho, weights, c1, c2, d1, d2, prec);
                
                if (arb_is_nonnegative(acb_realref(der))) {
                    acb_clear(der);
                    if (VERBOSE)
                        flint_printf("\n");
                    return linear_rounding_certify_relative_pairwise(bound, length, a1lo, a1hi, a2lo, a2lo, rholo, rhohi,
                                        weights, c1, c2, d1, d2, prec);
                } else if (arb_is_nonpositive(acb_realref(der))) {
                    acb_clear(der);
                    if (VERBOSE)
                        flint_printf("\n");
                    return linear_rounding_certify_relative_pairwise(bound, length, a1lo, a1hi, a2hi, a2hi, rholo, rhohi,
                                        weights, c1, c2, d1, d2, prec);
                } else {
                    if (VERBOSE)
                        flint_printf(" ");
                    acb_clear(der);
                } 
            }
        }
    
    }
    
   
    
    int good = rounding_linear_check_relative_pairwise(bound, length, a1, a2, rho, weights, c1, c2, d1, d2, prec);
    
    // flint_printf("Not_nan: %d\n", not_nan);
    
    acb_clear(a1), acb_clear(a2), acb_clear(rho);

    if (good == 1) {
        if (VERBOSE) flint_printf(" GOOD\n");
        return 1;
    }

    if (good == 2) {
        if (VERBOSE) flint_printf(" INVALID\n");
        return 1;
    }
        
    if (good == -1) {
        if (VERBOSE) flint_printf(" BAD\n");
        return 0;
    }

    double a1len = a1hi - a1lo;

    double a2len = a2hi - a2lo;

    double rholen = rhohi - rholo;
    
    

    if (a1len >= a2len && a1len >= rholen) {
        if (VERBOSE) flint_printf(" SPLIT A1\n");
        double a1mid = (a1lo + a1hi) / 2.0;
        good = linear_rounding_certify_relative_pairwise(bound, length, a1lo, a1mid, a2lo, a2hi, rholo, rhohi,
                                       weights, c1, c2, d1, d2, prec) &&
            linear_rounding_certify_relative_pairwise(bound, length, a1mid, a1hi, a2lo, a2hi, rholo, rhohi,
                                    weights, c1, c2, d1, d2, prec);  
                                    
        return good;
        
    } else if (a2len >= rholen) {
        if (VERBOSE) flint_printf(" SPLIT A2\n");
        double a2mid = (a2lo + a2hi) / 2.0;
        good = linear_rounding_certify_relative_pairwise(bound, length, a1lo, a1hi, a2lo, a2mid, rholo, rhohi,
                                       weights, c1, c2, d1, d2, prec) &&
            linear_rounding_certify_relative_pairwise(bound, length, a1lo, a1hi, a2mid, a2hi, rholo, rhohi,
                                    weights, c1, c2, d1, d2, prec);      
        return good;
    } else {
        if (VERBOSE) flint_printf(" SPLIT RHO\n");
        double rhomid = (rholo + rhohi) / 2.0;
        good = linear_rounding_certify_relative_pairwise(bound, length, a1lo, a1hi, a2lo, a2hi, rholo, rhomid,
                                       weights, c1, c2, d1, d2, prec) &&
            linear_rounding_certify_relative_pairwise(bound, length, a1lo, a1hi, a2lo, a2hi, rhomid, rhohi,
                                    weights, c1, c2, d1, d2, prec);      
        return good;
    }
}

void lower_bound_dicut_engine(int ii, int jj, double dbound) {
    slong prec = 32;

    int length, brkpts;
    
    // read input from file

    flint_scanf("%d %d", &brkpts, &length);
    flint_printf("%d %d\n", brkpts, length);

    arb_t bound;
    arb_init(bound);
    arb_set_d(bound, dbound);

    
    double brkpts_list[brkpts];
    
    for (int i = 0; i < brkpts; i++) {
        flint_scanf("%lf", &brkpts_list[i]);
        flint_printf("%d: %lf\n", i, brkpts_list[i]);
    }
    double dw[length], thresh[length][brkpts];

    for (int i = 0; i < length; i++) {
        flint_scanf("%lf", &dw[i]);
        flint_printf("%lf:", dw[i]);
        for (int j = 0; j < brkpts; j++) {
            flint_scanf("%lf", &thresh[i][j]);
            flint_printf(" %lf", thresh[i][j]);
        }
        flint_printf("\n");
    }

    acb_t weights[length], c1[length], c2[length], d1[length], d2[length];
    
    for (int i = 0; i < length; i++) {
        acb_init(weights[i]);
        acb_set_d(weights[i], dw[i]);
        acb_init(c1[i]);
        acb_init(c2[i]);
        acb_init(d1[i]);
        acb_init(d2[i]);
    }

    assert (ii >= 0 && ii < brkpts - 1);
    assert (jj >= 0 && jj < brkpts - 1);

    for (int i = ii; i <= ii; i++) {
        for (int j = jj; j <= jj; j++) {
            for (int k = 0; k < length; k++) {
                double z1 = (thresh[k][i+1] - thresh[k][i]) / (brkpts_list[i+1] - brkpts_list[i]);
                double z2 = (thresh[k][j+1] - thresh[k][j]) / (brkpts_list[j+1] - brkpts_list[j]);
                acb_set_d(c1[k], z1);
                acb_set_d(c2[k], z2);
                acb_set_d(d1[k], thresh[k][i] - z1 * brkpts_list[i]);
                acb_set_d(d2[k], thresh[k][j] - z2 * brkpts_list[j]);
            }

            flint_printf("CHECKING: %d x %d\n", i, j);

            double alo1 = brkpts_list[i];
            double ahi1 = brkpts_list[i+1];
            double alo2 = brkpts_list[j];
            double ahi2 = brkpts_list[j+1];

            flint_printf("(%lf, %lf) x (%lf, %lf)\n", alo1, ahi1, alo2, ahi2);

            //double b12lo = -1, b12hi = +1;

            //int good = linear_rounding_certify(bound, length, alo1, ahi1, alo2, ahi2, b12lo, b12hi,
            //                           weights, c1, c2, d1, d2, prec, NULL, NULL, 0);
            double rholo = -1, rhohi = +1;

            int good = linear_rounding_certify_relative_pairwise(bound, length, alo1, ahi1, alo2, ahi2, rholo, rhohi,
                                       weights, c1, c2, d1, d2, prec);

            flint_printf("RESULT %d x %d: %d in %lf s\n", i, j, good, 1.0 * clock() / CLOCKS_PER_SEC);
            fflush(stdout);
            assert(good);
            if (!good)
                goto a;

        }
    }
    
    a:
    // cleanup

    for (int i = 0; i < length; i++) {
        acb_clear(weights[i]);
        acb_clear(c1[i]);
        acb_clear(c2[i]);
        acb_clear(d1[i]);
        acb_clear(d2[i]);
    }

    arb_clear(bound);
    
}
