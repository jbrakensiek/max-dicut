/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#ifndef MAX2SAT_HPP
#define MAX2SAT_HPP

#include "arb_wrapper.hpp"
#include "config.hpp"

class Max2Sat : public Config {
public:
    // https://stackoverflow.com/questions/347358/inheriting-constructors
    Max2Sat(const Arb b1, const Arb b2, const Arb b12) : Config(b1,b2,b12) { }
    Max2Sat(const Config &c) : Config(c) { }

    Arb value() const;

    Arb prob(const Arb& beta) const;
    Arb prob(const Arb& t1, const Arb& t2) const;

    // prob(beta) / value
    Arb ratio(const Arb& beta) const;

    // this prob(beta) - alpha * value -- avoids divide by 0..
    Arb obj(const Arb& alpha, const Arb& beta) const;

    static Arb value_from_rel(const Arb& b1, const Arb& b2, const Arb& rho);
    static Arb value_from_rel_d_b1(const Arb& b1, const Arb& b2, const Arb& rho);
    static Arb value_from_rel_d_b2(const Arb& b1, const Arb& b2, const Arb& rho);
    static Arb value_from_rel_d_rho(const Arb& b1, const Arb& b2, const Arb& rho);

    // for type 4/5
    static Arb prob_from_rel(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta);
    static Arb prob_from_rel_d_b1(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta);
    static Arb prob_from_rel_d_b2(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta);
    static Arb prob_from_rel_d_rho(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta);

    // for type 3
    static Arb type3_prob_from_rel(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta);
    static Arb type3_prob_from_rel_d_b1(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta);
    static Arb type3_prob_from_rel_d_b2(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta);
    static Arb type3_prob_from_rel_d_rho(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta);

};

#endif
