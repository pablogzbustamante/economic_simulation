#include "sim/invariants.hpp"
#include "sim/util.hpp"

namespace sim {

void check_invariants(const SimState& st, int tick) {
  require_finite(st.agg.wage, "agg.wage", tick);
  require_finite(st.agg.unemployment, "agg.unemployment", tick);
  require_finite(st.agg.total_output, "agg.total_output", tick);
  require_finite(st.agg.total_employment_hours, "agg.total_employment_hours", tick);
  require_finite(st.agg.W, "agg.W", tick);
  require_finite(st.agg.S, "agg.S", tick);
  require_finite(st.agg.e, "agg.e", tick);
  require_finite(st.agg.C_over_V, "agg.C_over_V", tick);
  require_finite(st.agg.r, "agg.r", tick);

  require_nonneg(st.agg.wage, "agg.wage", tick);
  if (!(st.agg.unemployment >= -1e-12 && st.agg.unemployment <= 1.0 + 1e-12))
    die(std::string("invariant fail tick=") + std::to_string(tick) + " unemployment out of [0,1] value=" + std::to_string(st.agg.unemployment));

  for (int j = 0; j < st.S; ++j) {
    const auto& s = st.sectors[(std::size_t)j];
    require_finite(s.price, "sector.price", tick);
    require_nonneg(s.price, "sector.price", tick);
    require_finite(s.demand_final, "sector.demand_final", tick);
    require_finite(s.demand_total, "sector.demand_total", tick);
    require_finite(s.supply, "sector.supply", tick);
    require_finite(s.sales, "sector.sales", tick);
    require_finite(s.value_unit, "sector.value_unit", tick);
    require_nonneg(s.demand_final, "sector.demand_final", tick);
    require_nonneg(s.demand_total, "sector.demand_total", tick);
    require_nonneg(s.supply, "sector.supply", tick);
    require_nonneg(s.sales, "sector.sales", tick);
  }

  for (int i = 0; i < st.N; ++i) {
    const auto& f = st.firms[(std::size_t)i];
    require_finite(f.K, "firm.K", tick);
    require_finite(f.inv, "firm.inv", tick);
    require_finite(f.L, "firm.L", tick);
    require_finite(f.q_plan, "firm.q_plan", tick);
    require_finite(f.q_eff, "firm.q_eff", tick);
    require_finite(f.sales, "firm.sales", tick);
    require_finite(f.revenue, "firm.revenue", tick);
    require_finite(f.cost_intermediate, "firm.cost_intermediate", tick);
    require_finite(f.wage_bill, "firm.wage_bill", tick);
    require_finite(f.profit, "firm.profit", tick);

    require_nonneg(f.K, "firm.K", tick);
    require_nonneg(f.inv, "firm.inv", tick);
    require_nonneg(f.L, "firm.L", tick);
    require_nonneg(f.q_plan, "firm.q_plan", tick);
    require_nonneg(f.q_eff, "firm.q_eff", tick);
    require_nonneg(f.sales, "firm.sales", tick);
    require_nonneg(f.revenue, "firm.revenue", tick);
    require_nonneg(f.cost_intermediate, "firm.cost_intermediate", tick);
    require_nonneg(f.wage_bill, "firm.wage_bill", tick);
  }
}

}