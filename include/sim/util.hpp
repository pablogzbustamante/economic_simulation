#pragma once
#include <string>
#include "sim/constants.hpp"

namespace sim {

    double clamp_min(double x, double lo);
    double safe_div(double num, double den, double fallback = 0.0);
    bool is_finite(double x);
    void die(const std::string& msg);

    void require_finite(double x, const char* name, int tick);
    void require_nonneg(double x, const char* name, int tick);

}