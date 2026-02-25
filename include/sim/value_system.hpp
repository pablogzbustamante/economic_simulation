#pragma once

#include "sim/types.hpp"

namespace sim {

    Vec solve_labor_values(const Mat& A, const Vec& l, double tol = 1e-10, int max_iter = 100000);

}