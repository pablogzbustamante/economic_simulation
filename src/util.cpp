#include "sim/util.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>

namespace sim {

    double clamp_min(double x, double lo) {
    return (x < lo) ? lo : x;
    }

    double safe_div(double num, double den, double fallback) {
    if (!std::isfinite(num) || !std::isfinite(den)) return fallback;
    if (std::abs(den) < EPS) return fallback;
    return num / den;
    }

    bool is_finite(double x) {
    return std::isfinite(x);
    }

    void die(const std::string& msg) {
    std::cerr << msg << "\n";
    std::exit(1);
    }

    void require_finite(double x, const char* name, int tick) {
    if (!std::isfinite(x)) die(std::string("invariant fail tick=") + std::to_string(tick) + " non-finite " + name);
    }

    void require_nonneg(double x, const char* name, int tick) {
    if (!(x >= -1e-12)) die(std::string("invariant fail tick=") + std::to_string(tick) + " negative " + name + " value=" + std::to_string(x));
    }

}