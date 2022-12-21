/*
  Copyright (c) 2022 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#ifndef LOWER_H
#define LOWER_H

#include <assert.h>
#include <math.h>
#include <time.h>
#include "arb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"
#include "common.h"

#define DERIVATIVE_THRESHOLD_RHO 0.005 // arbitrarily chosen
#define DERIVATIVE_THRESHOLD_A 0.002 // arbitrarily chosen

int tri_check(const acb_t a1, const acb_t a2, const acb_t b12, slong prec);
int tri_check_rho(const acb_t a1, const acb_t a2, const acb_t rho, slong prec);


void partial_derivative_a1(acb_t res, arb_t bound, const ulong length, const acb_t a1, const acb_t a2, const acb_t rho, const acb_t* weights,
                            const acb_t* c1, const acb_t* c2,
                            const acb_t* d1, const acb_t* d2, slong prec);

void partial_derivative_a2(acb_t res, arb_t bound, const ulong length, const acb_t a1, const acb_t a2, const acb_t rho, const acb_t* weights,
                            const acb_t* c1, const acb_t* c2,
                            const acb_t* d1, const acb_t* d2, slong prec);

void partial_derivative_rho(acb_t res, arb_t bound, const ulong length, const acb_t a1, const acb_t a2, const acb_t rho, const acb_t* weights,
                            const acb_t* c1, const acb_t* c2,
                            const acb_t* d1, const acb_t* d2, slong prec);
                            
                            
int rounding_linear_check_relative_pairwise(const arb_t bound, const ulong length,
                          const acb_t a1, const acb_t a2, const acb_t rho,
                          const acb_t* weights, const acb_t* c1, const acb_t* c2,
                          const acb_t* d1, const acb_t* d2, slong prec);

int linear_rounding_certify_relative_pairwise(arb_t bound, const ulong length, 
                            double a1lo, double a1hi, double a2lo, double a2hi,
                            double rholo, double rhohi, const acb_t* weights,
                            const acb_t* c1, const acb_t* c2,
                            const acb_t* d1, const acb_t* d2, slong prec);

void lower_bound_dicut_engine(int ii, int jj, double dbound);

#endif
