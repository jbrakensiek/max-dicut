/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include "maxdicut.hpp"
#include "bivariate_normal.hpp"
#include <cassert>

// (1 + x_1 - x_2 - x_1x_2)/4

Arb MaxDiCut::value() const {
    return (1 + (this -> b1) - (this -> b2) - (this -> b12)) / 4;
}

// t_1, t_2 \in [0, 1]
Arb MaxDiCut::prob_normalized(const Arb& t1, const Arb& t2) const {
    Arb relrho = this -> rho_safe();
    return biv_norm_cdf_norm_thresh(t1, 1 - t2, -relrho);
}

// t_1, t_2 \in \mathbb{R}
Arb MaxDiCut::prob_unnormalized(const Arb& t1, const Arb& t2) const {
    Arb relrho = this -> rho_safe();
    return biv_norm_cdf(t1, -t2, -relrho);
}

// t_1, t_2 \in \mathbb{R}
Arb MaxDiCut::prob_unnormalized_from_rho(const Arb& t1, const Arb& t2, const Arb& rho) const {
    return biv_norm_cdf(t1, -t2, -rho);
}

Arb MaxDiCut::value_from_rel(const Arb& b1, const Arb& b2, const Arb& rho) {
    Arb b12 = b12_from_rel_rho(b1, b2, rho);
    return (1 + b1 - b2 - b12) / 4;
}

Arb MaxDiCut::value_from_rel_d_b1(const Arb& b1, const Arb& b2, const Arb& rho) {
    Arb x = (1-b2)/4;
    Arb y = rho * b1 * Arb::safe_sqrt(1-b2.sqr()) / (4 * Arb::safe_sqrt(1 - b1.sqr()));
    return x + y;
}

Arb MaxDiCut::value_from_rel_d_b2(const Arb& b1, const Arb& b2, const Arb& rho) {
    Arb x = (-1-b2)/4;
    Arb y = rho * b1 * Arb::safe_sqrt(1-b2.sqr()) / (4 * Arb::safe_sqrt(1 - b1.sqr()));
    return x + y;
}

Arb MaxDiCut::value_from_rel_d_rho(const Arb& b1, const Arb& b2, const Arb& rho) {
    return -Arb::safe_sqrt((1-b1.sqr())*(1-b2.sqr())) / 4;
}

// returns 1 if at least one of the triangle inequalities is invalid
// assumes that b1 and b2 are bounded away from 1 and -1
int MaxDiCut::is_invalid_tri_check_from_rho(const Arb& b1, const Arb& b2, const Arb& rho) { 
    int tri_1 = rho * Arb::safe_sqrt( (1 - b1) * (1 - b2)) + Arb::safe_sqrt( (1 + b1) * (1 + b2)) < 0;
    int tri_2 = rho * Arb::safe_sqrt( (1 + b1) * (1 + b2)) + Arb::safe_sqrt( (1 - b1) * (1 - b2)) < 0;
    int tri_3 = -rho * Arb::safe_sqrt( (1 + b1) * (1 - b2)) + Arb::safe_sqrt( (1 - b1) * (1 + b2)) < 0;
    int tri_4 = -rho * Arb::safe_sqrt( (1 - b1) * (1 + b2)) + Arb::safe_sqrt( (1 + b1) * (1 - b2)) < 0;
    return (tri_1 || tri_2 || tri_3 || tri_4);
}
