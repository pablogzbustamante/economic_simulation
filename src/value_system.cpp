#include "sim/value_system.hpp"
#include "sim/util.hpp"
#include <cmath>
#include <algorithm>

namespace sim {

Vec solve_labor_values(const Mat& A, const Vec& l, double tol, int max_iter) {
  const int S = (int)l.size();
  if ((int)A.size() != S) die("solve_labor_values: A rows != S");
  for (int i = 0; i < S; ++i) if ((int)A[i].size() != S) die("solve_labor_values: A cols != S");

  Vec v((std::size_t)S, 0.0);
  Vec v_next((std::size_t)S, 0.0);

  for (int it = 0; it < max_iter; ++it) {
    for (int j = 0; j < S; ++j) {
      double acc = l[(std::size_t)j];
      for (int i = 0; i < S; ++i) acc += A[(std::size_t)i][(std::size_t)j] * v[(std::size_t)i];
      v_next[(std::size_t)j] = acc;
    }

    double md = 0.0;
    for (int j = 0; j < S; ++j) md = std::max(md, std::abs(v_next[(std::size_t)j] - v[(std::size_t)j]));

    v.swap(v_next);

    if (!is_finite(md)) die("solve_labor_values: NaN/Inf during iteration");
    if (md < tol) return v;
  }

  die("solve_labor_values: did not converge");
  return v;
}

}