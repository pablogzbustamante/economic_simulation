#pragma once
#include <fstream>
#include <string>
#include "sim/types.hpp"
#include "sim/state.hpp"

namespace sim {

struct CsvWriter {
  std::ofstream out;
  i32 S;

  CsvWriter(const std::string& path, i32 S_);
  void write_header();
  void write_row(i32 tick,
                 f64 wage,
                 f64 unemployment,
                 f64 output_total,
                 f64 employment_hours,
                 f64 W,
                 f64 S_val,
                 f64 e,
                 f64 C_over_V,
                 f64 r,
                 const Vec& prices);

  void write_row_from_state(i32 tick, const SimState& st);
};

}