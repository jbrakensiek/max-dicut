/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <vector>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include "config.hpp"
#include "maxdicut.hpp"
#define NUM_THREADS 1
#define COMP_LIMIT 0.000099

class Linear_scheme {
public:
    int length;
    Arb total_weight;
    std::vector<Arb> weights;
    std::vector<Arb> c1, d1, c2, d2;   // unnormalized thresholds: c1 * b1 + d1, c2 * b2 + d2
};

Arb prob_by_linear_scheme_from_rel(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s) {
    Arb result = 0;
    Arb b12 = Config::b12_from_rel_rho(b1, b2, rho);
    MaxDiCut conf(b1, b2, b12);
    for (int i = 0; i < s.length; i++) {
        Arb t1 = s.c1[i] * b1 + s.d1[i];
        Arb t2 = s.c2[i] * b2 + s.d2[i];
        result = result + s.weights[i] * conf.prob_unnormalized_from_rho(t1, t2, rho);
    }
    return result;
}

Arb obj(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound) {
    // sound - bound * comp
    // (\sum_i w_i * sound_i) - bound * total_w * comp
    Arb b12 = Config::b12_from_rel_rho(b1, b2, rho);
    MaxDiCut conf(b1, b2, b12);
    Arb result = prob_by_linear_scheme_from_rel(b1, b2, rho, s) - bound * s.total_weight * conf.value();
    return result;
}

Arb obj_d_b1(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound) {
    Arb result = 0;
    for (int i = 0; i < s.length; i++) {
        Arb t1 = s.c1[i] * b1 + s.d1[i];
        Arb t2 = s.c2[i] * b2 + s.d2[i];
        // partial derivative = c1 * normalPDF(t1) * normalCDF((-t2 + rho * t1) / sqrt(1 - rho * rho)) - bound / 4 * (1 - b_2 + rho * sqrt( (1 - b2 * b2) / (1 - b1 * b1) ) * b1)
        Arb d_sound = s.c1[i] * Arb::norm_pdf(t1) * Arb::norm_cdf((-t2 + rho * t1) / Arb::safe_sqrt(1 - rho * rho));
        result = result + s.weights[i] * ( d_sound - bound * MaxDiCut::value_from_rel_d_b1(b1, b2, rho) );
    }
    return result;
}

Arb obj_d_b2(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound) {
    Arb result = 0;
    for (int i = 0; i < s.length; i++) {
        Arb t1 = s.c1[i] * b1 + s.d1[i];
        Arb t2 = s.c2[i] * b2 + s.d2[i];
        // partial derivative = -c2 * normalPDF(-t2) * normalCDF((t1 - rho * t2) / sqrt(1 - rho * rho)) + bound / 4 * (1 + b_1 - rho * sqrt( (1 - b1 * b1) / (1 - b2 * b2) ) * b2)
        Arb d_sound = -s.c2[i] * Arb::norm_pdf(t2) * Arb::norm_cdf((t1 - rho * t2) / Arb::safe_sqrt(1 - rho * rho));
        result = result + s.weights[i] * ( d_sound - bound * MaxDiCut::value_from_rel_d_b2(b1, b2, rho) );
    }
    return result;
}

Arb obj_d_rho(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound) {
    Arb result = 0;
    for (int i = 0; i < s.length; i++) {
        Arb t1 = s.c1[i] * b1 + s.d1[i];
        Arb t2 = s.c2[i] * b2 + s.d2[i];
        // partial derivative * 2pi * sqrt(1 - rho^2) = - exp( - (t_1^2 - 2*rho*t1*t2 + t2^2) / 2(1 - rho^2)) + pi/2 * sqrt(1 - rho^2) * bound * sqrt( (1 - a1^2) * (1 - a2^2)).
        Arb d_sound = -Arb::exp( - (t1 * t1 - 2 * rho * t1 * t2 + t2 * t2) / 2 / (1 - rho * rho)) / 2 / Arb::pi() / Arb::sqrt(1 - rho*rho);
        result = result + s.weights[i] * ( d_sound - bound * MaxDiCut::value_from_rel_d_rho(b1, b2, rho) );
    }
    return result;
}


int check_linear_scheme_from_rel(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound) {
    Arb result = 0;
    Arb b12 = Config::b12_from_rel_rho(b1, b2, rho);
    MaxDiCut conf(b1, b2, b12);

    result = prob_by_linear_scheme_from_rel(b1, b2, rho, s) - bound * s.total_weight * conf.value();
    if (result >= 0)
        return 1; // success
    else if (result < 0)
        return -1; // failure
    else
        return 0; // uncertain
}

// -- modified from 2-SAT code
int check(const Arb &b1, const Arb &b2, const Arb &rho, const Linear_scheme &s, const Arb &bound) {

    Arb b12 = Config::b12_from_rel_rho(b1, b2, rho);

    if (MaxDiCut::is_invalid_tri_check_from_rho(b1, b2, rho)) {
        // invalid region, so "good" by default
        return 1;
    }

    if (obj(b1, b2, rho, s, bound) >= 0) {
        // vol_est = vol_est + vol(b1, b2, rho);
        return 1;
    } else if (obj(b1, b2, rho, s, bound) < 0) {
        // current implementation may introduce false negative
        // since obj only computes a lower bound
        return 0;
    }

    if ((b12 < 1 - Arb::abs(b1 - b2)) && (b12 > -1 + Arb::abs(b1 + b2))) {
        // derivative check only inside valid regions
        Arb d_b1 = obj_d_b1(b1, b2, rho, s, bound);
        Arb d_b2 = obj_d_b2(b1, b2, rho, s, bound);
        Arb d_rho = obj_d_rho(b1, b2, rho, s, bound);

#ifdef DEBUG
        flint_printf("PARTIALS\n");
        d_b1.pretty_println();
        d_b2.pretty_println();
        d_rho.pretty_println();
#endif

        if (!b1.is_nan() && (d_b1 > 0 || d_b1 < 0)) {
            //excl_est = excl_est + vol(b1, b2, rho);
            return 1;
        }
        else if (!b2.is_nan() && (d_b2 > 0 || d_b2 < 0)) {
            //excl_est = excl_est + vol(b1, b2, rho);
            return 1;
        }
        else if (!rho.is_nan() && (d_rho > 0 || d_rho < 0)) {
            //excl_est = excl_est + vol(b1, b2, rho);
            return 1;
        }
    }

    if (MaxDiCut::value_from_rel(b1, b2, rho) <= COMP_LIMIT) {
        // completeness is too small
        return 1;
    }

    // otherwise we need to split

    Arb rb1 = b1.rad();
    Arb rb2 = b2.rad();
    Arb rrho = rho.rad();

    if (rb1 >= rb2 && rb1 >= rrho) {
        return check(b1.left_half(), b2, rho, s, bound) && check(b1.right_half(), b2, rho, s, bound);
    }
    else if (rb2 >= rrho) {
        return check(b1, b2.left_half(), rho, s, bound) &&
                        check(b1, b2.right_half(), rho, s, bound);
    }
    else {
        return check(b1, b2, rho.left_half(), s, bound) &&
                        check(b1, b2, rho.right_half(), s, bound);
    }
    return 1;
}


int main(int argc, char* argv[]) {
    // for multi-threaded runs
    flint_set_num_threads(NUM_THREADS);

    printf("START: %d %s %s %s\n", argc, argv[1], argv[2], argv[3]);

    int ii = atoi(argv[1]), jj = atoi(argv[2]);
    Arb bound(atof(argv[3]));

    //reading input:

    int length, brkpts;
    flint_scanf("%d %d", &brkpts, &length);
    flint_printf("%d %d\n", brkpts, length);


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

    Linear_scheme s;
    s.length = length;
    s.total_weight = 0;

    assert ((ii >= 0) && (ii < brkpts - 1));
    assert ((jj >= 0) && (jj < brkpts - 1));

    for (int i = ii; i <= ii; i++) {
        for (int j = jj; j <= jj; j++) {
            for (int k = 0; k < length; k++) {
                Arb z1 = (thresh[k][i+1] - thresh[k][i]) / (brkpts_list[i+1] - brkpts_list[i]);
                Arb z2 = (thresh[k][j+1] - thresh[k][j]) / (brkpts_list[j+1] - brkpts_list[j]);
                s.c1.push_back(z1);
                s.c2.push_back(z2);
                s.d1.push_back(thresh[k][i] - z1 * brkpts_list[i]);
                s.d2.push_back(thresh[k][j] - z2 * brkpts_list[j]);
                s.weights.push_back(dw[k]);
                s.total_weight = s.total_weight + dw[k];
            }

            flint_printf("CHECKING: %d x %d\n", i, j);

            Arb b1 = Arb::join(brkpts_list[i], brkpts_list[i+1]);
            Arb b2 = Arb::join(brkpts_list[j], brkpts_list[j+1]);
            Arb rho = Arb::join(-1, 1);

            b1.println();
            b2.println();

            int good = check(b1, b2, rho, s, bound);

            flint_printf("RESULT %d x %d: %d in %lf s\n", i, j, good, 1.0 * clock() / CLOCKS_PER_SEC);
            fflush(stdout);
            assert(good);
        }
    }

    flint_cleanup_master();

    return 0;
}
