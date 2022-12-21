/*
  Copyright (c) 2022 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#ifndef COMMON_H
#define COMMON_H

#include <assert.h>
#include <math.h>
#include <time.h>
#include "arb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"

#define RHO_LIMIT 0.999 
#define VERBOSE 1
#define COMP_LIMIT 0.000001

// simple functions

int acb_is_nan(const acb_t z);

int arb_is_nan(const arb_t z);

void arb_union_d(arb_t z, double a, double b, slong prec);

void acb_union_d(acb_t z, double a, double b, slong prec);

void acb_norm_cdf(acb_t res, const acb_t z, slong prec);

void acb_norm_2d_cdf_zero_mean(acb_t res, const acb_t rho, slong prec);

int acb_crop_d(acb_t res, const acb_t z, double lo, double hi, slong prec);

void b12_from_rho(acb_t res, const acb_t a1, const acb_t a2, const acb_t rho, slong prec);
// configs

typedef struct config
{
    ulong i1; // index into a table of biases
    ulong i2; 
    acb_t b12; // an explicit pairwise bias
} config;

typedef config* config_t;

config_t config_init(const ulong i1, const ulong i2, const acb_t b12);

config_t config_init_d(const ulong i1, const ulong i2, const double b12);

int config_clear(config_t c);

// threshold computations

int thresh_eval_helper(acb_ptr res, const acb_t z, void * param, slong order, slong prec);

void thresh_eval_relative(acb_ptr res, const acb_t rho, const acb_t t1r, const acb_t t2r, slong prec);

// completeness

int value(acb_ptr res, const acb_t a1, const acb_t a2, const acb_t b12, slong prec);

void value_from_rho(acb_t res, const acb_t a1, const acb_t a2, const acb_t rho, slong prec);

#endif
