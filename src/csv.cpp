#include "sim/csv.hpp"
#include "sim/util.hpp"
#include <filesystem>
#include <iomanip>

namespace sim {

CsvWriter::CsvWriter(const std::string& path, i32 S_) : S(S_) {
  std::filesystem::path p(path);
  if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());
  out.open(path, std::ios::out | std::ios::trunc);
  if (!out.is_open()) die("cannot open output csv: " + path);
  out.setf(std::ios::fixed);
  out << std::setprecision(10);
}

void CsvWriter::write_header() {
  out << "tick,wage,unemployment,output_total,employment_hours,W,S,e,C_over_V,r";
  for (i32 j = 0; j < S; ++j) out << ",price_" << j;
  out << "\n";
}

void CsvWriter::write_row(i32 tick, f64 wage, f64 unemployment,f64 output_total, f64 employment_hours, f64 W, f64 S_val, f64 e, f64 C_over_V, f64 r, const Vec& prices) {
  if ((i32)prices.size() != S) die("prices size != S");
  out << tick << "," << wage << "," << unemployment << "," << output_total << ","
      << employment_hours << "," << W << "," << S_val << "," << e << "," << C_over_V << "," << r;
  for (i32 j = 0; j < S; ++j) out << "," << prices[j];
  out << "\n";
}

void CsvWriter::write_row_from_state(i32 tick, const SimState& st) {
  Vec prices;
  prices.reserve((std::size_t)S);
  for (i32 j = 0; j < S; ++j) prices.push_back(st.sectors[(std::size_t)j].price);

  write_row(tick, st.agg.wage, st.agg.unemployment, st.agg.total_output, st.agg.total_employment_hours, st.agg.W, st.agg.S, st.agg.e, st.agg.C_over_V, st.agg.r, prices);
}

}