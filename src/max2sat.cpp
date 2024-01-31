/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include "max2sat.hpp"
#include "bivariate_normal.hpp"
#include <cassert>

Arb Max2Sat::value() const {
    return (3 - (this -> b1) - (this -> b2) - (this -> b12)) / 4;
}

Arb Max2Sat::prob(const Arb& beta) const {
    Arb t1 = (1 + beta*(this->b1))/2;
    Arb t2 = (1 + beta*(this->b2))/2;
    return this -> prob(t1, t2);
}

Arb Max2Sat::prob(const Arb& t1, const Arb& t2) const {
    Arb relrho = this -> rho_safe();
    return 1 - biv_norm_cdf_norm_thresh(t1, t2, relrho); 
}

Arb Max2Sat::ratio(const Arb& beta) const {
    return this -> prob(beta) / this -> value();
}

Arb Max2Sat::obj(const Arb& alpha, const Arb& beta) const {
    return this -> prob(beta) - alpha * this -> value();
}

Arb Max2Sat::value_from_rel(const Arb& b1, const Arb& b2, const Arb& rho) {
    Arb b12 = b12_from_rel_rho(b1, b2, rho);
    return (3 - b1 - b2 - b12) / 4;
}

Arb Max2Sat::value_from_rel_d_b1(const Arb& b1, const Arb& b2, const Arb& rho) {
    Arb x = -(b2+1)/4;
    Arb y = rho * b1 * Arb::safe_sqrt(1-b2.sqr()) / (4 * Arb::safe_sqrt(1 - b1.sqr()));
    return x + y;
}

Arb Max2Sat::value_from_rel_d_b2(const Arb& b1, const Arb& b2, const Arb& rho) {
    Arb x = -(b1+1)/4;
    Arb y = rho * b2 * Arb::safe_sqrt(1-b1.sqr()) / (4 * Arb::safe_sqrt(1 - b2.sqr()));
    return x + y;
}

Arb Max2Sat::value_from_rel_d_rho(const Arb& b1, const Arb& b2, const Arb& rho) {
    return -Arb::safe_sqrt((1-b1.sqr())*(1-b2.sqr())) / 4;
}

Arb Max2Sat::prob_from_rel(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta) {
    Arb t1 = (1 + beta*b1)/2;
    Arb t2 = (1 + beta*b2)/2;
    return 1 - biv_norm_cdf_norm_thresh(t1, t2, rho); 
}

Arb Max2Sat::prob_from_rel_d_b1(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta) {
    Arb t1 = Arb::norm_cdf_inv((1 + beta*b1)/2);
    Arb t2 = Arb::norm_cdf_inv((1 + beta*b2)/2);
    Arb c = (t2 - rho * t1) / Arb::safe_sqrt(1 - rho*rho);
    return (-beta / 2) * c.norm_cdf();
}

Arb Max2Sat::prob_from_rel_d_b2(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta) {
    Arb t1 = Arb::norm_cdf_inv((1 + beta*b1)/2);
    Arb t2 = Arb::norm_cdf_inv((1 + beta*b2)/2);
    Arb c = (t1 - rho * t2) / Arb::safe_sqrt(1 - rho*rho);
    return (-beta / 2) * c.norm_cdf();
}

Arb Max2Sat::prob_from_rel_d_rho(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta) {
    Arb t1 = Arb::norm_cdf_inv((1 + beta*b1)/2);
    Arb t2 = Arb::norm_cdf_inv((1 + beta*b2)/2);
    Arb x = -2 * Arb::pi() * Arb::safe_sqrt(1 - rho*rho);
    Arb y = -(t1.sqr()+t2.sqr() - 2*rho*t1*t2)/(2*(1-rho.sqr()));
    return Arb::exp(y) / x;
}

Arb Max2Sat::type3_prob_from_rel(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta) {
    Arb t1 = beta*(1 + b1)/2;
    Arb t2 = beta*(1 + b2)/2;
    return 1 - biv_norm_cdf_norm_thresh(t1, t2, rho); 
}

Arb Max2Sat::type3_prob_from_rel_d_b1(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta) {
    Arb t1 = Arb::norm_cdf_inv(beta*(1 + b1)/2);
    Arb t2 = Arb::norm_cdf_inv(beta*(1 + b2)/2);
    Arb c = (t2 - rho * t1) / Arb::safe_sqrt(1 - rho*rho);
    return (-beta / 2) * c.norm_cdf();
}

Arb Max2Sat::type3_prob_from_rel_d_b2(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta) {
    Arb t1 = Arb::norm_cdf_inv(beta*(1 + b1)/2);
    Arb t2 = Arb::norm_cdf_inv(beta*(1 + b2)/2);
    Arb c = (t1 - rho * t2) / Arb::safe_sqrt(1 - rho*rho);
    return (-beta / 2) * c.norm_cdf();
}

Arb Max2Sat::type3_prob_from_rel_d_rho(const Arb& b1, const Arb& b2, const Arb& rho, const Arb& beta) {
    Arb t1 = Arb::norm_cdf_inv(beta*(1 + b1)/2);
    Arb t2 = Arb::norm_cdf_inv(beta*(1 + b2)/2);
    Arb x = -2 * Arb::pi() * Arb::safe_sqrt(1 - rho*rho);
    Arb y = -(t1.sqr()+t2.sqr() - 2*rho*t1*t2)/(2*(1-rho.sqr()));
    return Arb::exp(y) / x;
}
