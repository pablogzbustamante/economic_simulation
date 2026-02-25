#pragma once
#include "sim/types.hpp"
#include <vector>

namespace sim {

struct FirmState {
  i32 sector;
  f64 K;
  f64 inv;
  f64 L;
  f64 q_plan;
  f64 q_eff;
  f64 sales;
  f64 revenue;
  f64 cost_intermediate;
  f64 wage_bill;
  f64 profit;
};

struct SectorState {
  f64 price;
  f64 demand_final;
  f64 demand_total;
  f64 supply;
  f64 sales;
  f64 value_unit;
};

struct AggState {
  f64 wage;
  f64 unemployment;
  f64 total_output;
  f64 total_employment_hours;
  f64 W;
  f64 S;
  f64 e;
  f64 C_over_V;
  f64 r;
};

struct SimState {
  i32 S;
  i32 N;
  std::vector<FirmState> firms;
  std::vector<SectorState> sectors;
  AggState agg;
};

SimState make_state(i32 S, i32 N);

}