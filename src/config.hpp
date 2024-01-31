/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "arb_wrapper.hpp"
#include "bivariate_normal.hpp"

class Config {
public:
    // constructors
    Config(const Arb b1, const Arb b2, const Arb b12);
    Config(const Config &c);

    // pseudo-constructor
    // Here's rel_b12 is in the range [0,1] cleanly
    // parameterizing between the lower and upper limits
    // of b12 allowed by the triangle inequality
    static Config from_relative(const Arb& b1, const Arb& b2, const Arb& rel_b12);

    // destructor not necessary?

    // other utlities
    void print() const;
    void println() const;

    // TODO: a "split" utility to take a config and split it into cases
    
    Arb lower_b12() const;
    static Arb lower_b12(const Arb& b1, const Arb& b2);

    Arb upper_b12() const;
    static Arb upper_b12(const Arb& b1, const Arb& b2);

    // this can return three values
    // +1 if the triangle inequality is *always* satisfied
    // -1 if the triangle inequality is *never* satisfied
    // 0 if ambiguous (this can include NaN)
    int triangle_check() const;

    // relative correlation
    Arb rho() const;
    static Arb rho(const Arb& b1, const Arb& b2, const Arb& b12);

    // version of above that ALWAYS outputs a subset of [-1, 1],
    // even if rho() returns NaN
    Arb rho_safe() const;
    static Arb rho_safe(const Arb& b1, const Arb& b2, const Arb& b12);

    // TODO: partial derivatives of rho

    static Arb b12_from_rel_rho(const Arb& b1, const Arb& b2, const Arb& rho);

    static int tri_check_rel_rho(const Arb& b1, const Arb& b2, const Arb& rho);

    // internal elements
    Arb b1, b2, b12;
};

#endif
