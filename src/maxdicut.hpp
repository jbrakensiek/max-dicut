/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#ifndef MAXDICUT_HPP
#define MAXDICUT_HPP


#define RHO_RAD_LIMIT 0.001

#include "arb_wrapper.hpp"
#include "config.hpp"

class MaxDiCut : public Config {
public:
    // https://stackoverflow.com/questions/347358/inheriting-constructors
    MaxDiCut(const Arb b1, const Arb b2, const Arb b12) : Config(b1,b2,b12) { }
    MaxDiCut(const Config &c) : Config(c) { }

    Arb value() const;

    // t_1 and t_2 are normalized thresholds \in [0, 1]
    Arb prob_normalized(const Arb& t1, const Arb& t2) const;

    // t_1 and t_2 are unnormalized thresholds \in \mathbb{R}
    Arb prob_unnormalized(const Arb& t1, const Arb& t2) const;
    Arb prob_unnormalized_from_rho(const Arb& t1, const Arb& t2, const Arb& rho) const;

    static Arb value_from_rel(const Arb& b1, const Arb& b2, const Arb& rho);
    static Arb value_from_rel_d_b1(const Arb& b1, const Arb& b2, const Arb& rho);
    static Arb value_from_rel_d_b2(const Arb& b1, const Arb& b2, const Arb& rho);
    static Arb value_from_rel_d_rho(const Arb& b1, const Arb& b2, const Arb& rho);

    // returns 1 if at least one of the triangle inequalities is invalid
    // assumes that b1 and b2 are bounded away from 1 and -1
    static int is_invalid_tri_check_from_rho(const Arb& b1, const Arb& b2, const Arb& rho); 

};

#endif
