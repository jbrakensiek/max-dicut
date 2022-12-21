/*
  Copyright (c) 2022 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <assert.h>
#include <math.h>
#include <time.h>
#include "arb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"
#include "common.h"
#include "lower.h"
#define NUM_THREADS 1

int main(int argc, char* argv[])
{
    // for multi-threaded runs
    flint_set_num_threads(NUM_THREADS);

    printf("START: %d %s %s %s\n", argc, argv[1], argv[2], argv[3]);

    lower_bound_dicut_engine(atoi(argv[1]), atoi(argv[2]), atof(argv[3]));

    flint_cleanup_master();

    return 0;
}
