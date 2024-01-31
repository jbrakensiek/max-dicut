/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include "config.hpp"
#include <cassert>

Config::Config(const Arb b1, const Arb b2, const Arb b12) {
    this->b1 = b1;
    this->b2 = b2;
    this->b12 = b12;
}

Config::Config(const Config &c) {
    this->b1 = c.b1;
    this->b2 = c.b2;
    this->b12 = c.b12;
}

void Config::print() const {
    flint_printf("( ");
    (this -> b1).print();
    flint_printf(" , ");
    (this -> b2).print();
    flint_printf(" , ");
    (this -> b12).print();
    flint_printf(" )");
}

void Config::println() const {
    this -> print();
    flint_printf("\n");
}

Arb Config::lower_b12(const Arb& b1, const Arb& b2) {
    return Arb::abs(b1+b2) - 1.0;
}

Arb Config::lower_b12() const {
    return lower_b12(this->b1, this->b2);
}

Arb Config::upper_b12(const Arb& b1, const Arb& b2) {
    return 1.0 - Arb::abs(b1-b2);
}

Arb Config::upper_b12() const {
    return upper_b12(this->b1, this->b2);
}

Config Config::from_relative(const Arb& b1, const Arb& b2, const Arb& rel_b12) {
    Arb b12 = (1 - rel_b12) * lower_b12(b1, b2) +
                    rel_b12 * upper_b12(b1, b2);
    return Config(b1,b2,b12);
}

int Config::triangle_check() const {
    if ( this->b12 >= this -> lower_b12() &&
         this->b12 <= this -> upper_b12() ) {
        return 1;
    }
    else if ( this->b12 < this -> lower_b12() ||
              this->b12 > this -> upper_b12() ){
        return -1;
    }
    else {
        return 0;
    }
}

Arb Config::rho(const Arb& b1, const Arb& b2, const Arb& b12) {
    Arb num = b12 - b1*b2;
    Arb denom = Arb::sqrt((1-b1.sqr())*(1-b2.sqr()));
    return num/denom;
}


Arb Config::rho() const {
    return rho(this -> b1, this -> b2, this -> b12);
}

Arb Config::rho_safe(const Arb& b1, const Arb& b2, const Arb& b12) {
    Arb rh = rho(b1, b2, b12);
    if (rh.is_nan()) {
        return Arb(-1,1);
    }
    return rh.intersect(Arb(-1,1));
}

Arb Config::rho_safe() const {
    return rho_safe(this -> b1, this -> b2, this -> b12);
}

Arb Config::b12_from_rel_rho(const Arb& b1, const Arb& b2, const Arb& rho) {
    Arb z = Arb::safe_sqrt((1-b1.sqr())*(1-b2.sqr()));
    if (z.is_nan()) {
        return Arb(-1, 1);
    }
    return b1*b2 + rho * z;
}

int Config::tri_check_rel_rho(const Arb& b1, const Arb& b2, const Arb& rho) {
    Arb b12 = b12_from_rel_rho(b1, b2, rho);
    Config c(b1, b2, b12);
    return c.triangle_check();
}
