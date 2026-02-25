#include "sim/engine.hpp"
#include "sim/util.hpp"
#include "sim/value_system.hpp"
#include "sim/invariants.hpp"
#include "sim/csv.hpp"
#include "sim/thread_pool.hpp"
#include "sim/benchmark.hpp"
#include <algorithm>
#include <vector>
#include <cmath>

namespace sim {

DoubleState::DoubleState(i32 S, i32 N) : cur(make_state(S, N)), next(make_state(S, N)) {}

void DoubleState::swap_buffers() {
  std::swap(cur, next);
}

static void init_firms(const SimConfig&, SimState& st) {
  for (i32 i = 0; i < st.N; ++i) {
    auto& f = st.firms[(std::size_t)i];
    f.sector = i % st.S;
    f.K = 100.0;
    f.inv = 0.0;
    f.L = 0.0;
    f.q_plan = 0.0;
    f.q_eff = 0.0;
    f.sales = 0.0;
    f.revenue = 0.0;
    f.cost_intermediate = 0.0;
    f.wage_bill = 0.0;
    f.profit = 0.0;
  }
}

static void init_sectors(const SimConfig& cfg, SimState& st) {
  for (i32 j = 0; j < st.S; ++j) {
    auto& s = st.sectors[(std::size_t)j];
    s.price = cfg.init_price[(std::size_t)j];
    s.demand_final = 0.0;
    s.demand_total = 0.0;
    s.supply = 0.0;
    s.sales = 0.0;
    s.value_unit = 0.0;
  }
}

static void init_agg(const SimConfig& cfg, SimState& st) {
  st.agg.wage = cfg.init_wage;
  st.agg.unemployment = cfg.u_ref;
  st.agg.total_output = 0.0;
  st.agg.total_employment_hours = cfg.labor_supply * (1.0 - cfg.u_ref);
  st.agg.W = st.agg.wage * st.agg.total_employment_hours;
  st.agg.S = 0.0;
  st.agg.e = 0.0;
  st.agg.C_over_V = 0.0;
  st.agg.r = 0.0;
}

DoubleState init_sim(const SimConfig& cfg) {
  DoubleState ds(cfg.S, cfg.N_firms);
  init_firms(cfg, ds.cur);
  init_sectors(cfg, ds.cur);
  init_agg(cfg, ds.cur);

  auto v = solve_labor_values(cfg.A, cfg.l);
  for (i32 j = 0; j < ds.cur.S; ++j) ds.cur.sectors[(std::size_t)j].value_unit = v[(std::size_t)j];

  ds.next = ds.cur;
  return ds;
}

static void compute_firms_in_sector(const SimState& cur, std::vector<i32>& counts) {
  std::fill(counts.begin(), counts.end(), 0);
  for (const auto& f : cur.firms) counts[(std::size_t)f.sector] += 1;
  for (auto& c : counts) if (c <= 0) c = 1;
}

static void phase_1_plan_production(const SimConfig& cfg, const SimState& cur, SimState& next) {
  next = cur;

  std::vector<i32> firms_in_sector((std::size_t)cur.S, 0);
  compute_firms_in_sector(cur, firms_in_sector);

  Vec D_final((std::size_t)cur.S, 0.0);
  const f64 income = cur.agg.W;
  const f64 consumption = (1.0 - cfg.save_rate) * income;

  for (i32 j = 0; j < cur.S; ++j) {
    const f64 p = clamp_min(cur.sectors[(std::size_t)j].price, EPS);
    const f64 spend = consumption * cfg.household_c_share[(std::size_t)j];
    D_final[(std::size_t)j] = safe_div(spend, p, 0.0);
    next.sectors[(std::size_t)j].demand_final = D_final[(std::size_t)j];
  }

  for (i32 i = 0; i < cur.N; ++i) {
    auto& nf = next.firms[(std::size_t)i];
    const i32 j = nf.sector;

    const f64 cap = safe_div(cur.firms[(std::size_t)i].K, cfg.kappa_cap, 0.0);
    const f64 Dpf = safe_div(D_final[(std::size_t)j], (f64)firms_in_sector[(std::size_t)j], 0.0);
    const f64 inv_target = cfg.inv_target_theta * Dpf;

    const f64 desired = Dpf + (inv_target - cur.firms[(std::size_t)i].inv);
    nf.q_plan = std::clamp(desired, 0.0, cap);
  }
}

static void phase_2_aggregate_plans_and_demands(const SimConfig& cfg, const SimState& cur, SimState& next) {
  const i32 S = cur.S;

  Vec Q_plan((std::size_t)S, 0.0);
  for (const auto& f : next.firms) Q_plan[(std::size_t)f.sector] += f.q_plan;

  Vec D_inter((std::size_t)S, 0.0);
  for (i32 i = 0; i < S; ++i) {
    f64 acc = 0.0;
    for (i32 j = 0; j < S; ++j) acc += cfg.A[(std::size_t)i][(std::size_t)j] * Q_plan[(std::size_t)j];
    D_inter[(std::size_t)i] = acc;
  }

  for (i32 i = 0; i < S; ++i) {
    const f64 df = next.sectors[(std::size_t)i].demand_final;
    next.sectors[(std::size_t)i].demand_total = df + D_inter[(std::size_t)i];
  }
}

static void phase_3_labor_and_output(const SimConfig& cfg, const SimState& cur, SimState& next) {
  f64 L_req_total = 0.0;
  for (i32 i = 0; i < cur.N; ++i) {
    const i32 j = next.firms[(std::size_t)i].sector;
    const f64 L_req = cfg.l[(std::size_t)j] * next.firms[(std::size_t)i].q_plan;
    L_req_total += L_req;
  }

  f64 phi = 1.0;
  if (L_req_total > 0.0) phi = std::min(1.0, safe_div(cfg.labor_supply, L_req_total, 1.0));

  f64 L_emp_total = 0.0;
  f64 Q_total = 0.0;

  for (i32 i = 0; i < cur.N; ++i) {
    auto& nf = next.firms[(std::size_t)i];
    const i32 j = nf.sector;

    const f64 L_req = cfg.l[(std::size_t)j] * nf.q_plan;
    nf.L = phi * L_req;
    nf.q_eff = phi * nf.q_plan;

    L_emp_total += nf.L;
    Q_total += nf.q_eff;
  }

  next.agg.total_employment_hours = L_emp_total;
  next.agg.total_output = Q_total;
  next.agg.unemployment = 1.0 - safe_div(L_emp_total, cfg.labor_supply, 0.0);
  next.agg.unemployment = std::clamp(next.agg.unemployment, 0.0, 1.0);
}

static void phase_4_sales_inventory_prices(const SimConfig& cfg, const SimState& cur, SimState& next) {
  const i32 S = cur.S;

  Vec supply((std::size_t)S, 0.0);
  Vec avail((std::size_t)S, 0.0);

  for (i32 i = 0; i < cur.N; ++i) {
    const i32 j = next.firms[(std::size_t)i].sector;
    supply[(std::size_t)j] += next.firms[(std::size_t)i].q_eff;
    avail[(std::size_t)j] += cur.firms[(std::size_t)i].inv + next.firms[(std::size_t)i].q_eff;
  }

  for (i32 j = 0; j < S; ++j) {
    next.sectors[(std::size_t)j].supply = supply[(std::size_t)j];

    const f64 D = next.sectors[(std::size_t)j].demand_total;
    const f64 Asec = avail[(std::size_t)j];
    const f64 sales_sec = std::min(Asec, D);

    next.sectors[(std::size_t)j].sales = sales_sec;

    const f64 denom = std::max(Asec, EPS);
    for (i32 i = 0; i < cur.N; ++i) {
      auto& nf = next.firms[(std::size_t)i];
      if (nf.sector != j) continue;

      const f64 a_i = cur.firms[(std::size_t)i].inv + nf.q_eff;
      const f64 s_i = sales_sec * safe_div(a_i, denom, 0.0);

      nf.sales = s_i;
      nf.inv = std::max(0.0, a_i - s_i);
    }

    const f64 p_cur = clamp_min(cur.sectors[(std::size_t)j].price, EPS);
    const f64 sup = supply[(std::size_t)j];
    const f64 excess = safe_div((D - sup), std::max(sup, EPS), 0.0);
    const f64 p_next = p_cur * (1.0 + cfg.price_adj * excess);
    next.sectors[(std::size_t)j].price = clamp_min(p_next, EPS);
  }
}

static void phase_5_accounting_investment(const SimConfig& cfg, const SimState& cur, SimState& next) {
  const i32 S = cur.S;
  Vec prices((std::size_t)S, 0.0);
  for (i32 j = 0; j < S; ++j) prices[(std::size_t)j] = clamp_min(cur.sectors[(std::size_t)j].price, EPS);

  for (i32 i = 0; i < cur.N; ++i) {
    auto& nf = next.firms[(std::size_t)i];
    const auto& cf = cur.firms[(std::size_t)i];
    const i32 j = nf.sector;

    nf.revenue = prices[(std::size_t)j] * nf.sales;

    f64 ic = 0.0;
    for (i32 k = 0; k < S; ++k) ic += prices[(std::size_t)k] * cfg.A[(std::size_t)k][(std::size_t)j] * nf.q_eff;
    nf.cost_intermediate = ic;

    nf.wage_bill = cur.agg.wage * nf.L;

    nf.profit = nf.revenue - nf.cost_intermediate - nf.wage_bill;

    const f64 I = cfg.inv_rate * std::max(nf.profit, 0.0);
    const f64 K_next = (1.0 - cfg.delta) * cf.K + I;
    nf.K = std::max(0.0, K_next);
  }
}

static void phase_6_metrics_and_wage_update(const SimConfig& cfg, const SimState& cur, SimState& next) {
  f64 W = 0.0;
  f64 C_inter = 0.0;
  f64 K_dep = 0.0;

  for (i32 i = 0; i < cur.N; ++i) {
    W += next.firms[(std::size_t)i].wage_bill;
    C_inter += next.firms[(std::size_t)i].cost_intermediate;
    K_dep += cfg.delta * cur.firms[(std::size_t)i].K;
  }

  const f64 V = W;
  const f64 C = C_inter + K_dep;
  const f64 L_emp = next.agg.total_employment_hours;

  f64 S_val = L_emp - W;
  if (!is_finite(S_val)) S_val = 0.0;

  next.agg.W = W;
  next.agg.S = S_val;
  next.agg.e = safe_div(S_val, std::max(W, EPS), 0.0);
  next.agg.C_over_V = safe_div(C, std::max(V, EPS), 0.0);
  next.agg.r = safe_div(S_val, std::max(C + V, EPS), 0.0);

  const f64 u = next.agg.unemployment;
  const f64 w_next = cur.agg.wage * (1.0 + cfg.wage_adj * (cfg.u_ref - u));
  next.agg.wage = clamp_min(w_next, EPS);
}

std::uint64_t run_single_thread(const SimConfig& cfg, DoubleState& ds, const std::string& out_csv_path) {
  CsvWriter writer(out_csv_path, cfg.S);
  writer.write_header();

  for (i32 t = 0; t < cfg.ticks; ++t) {
    phase_1_plan_production(cfg, ds.cur, ds.next);
    phase_2_aggregate_plans_and_demands(cfg, ds.cur, ds.next);
    phase_3_labor_and_output(cfg, ds.cur, ds.next);
    phase_4_sales_inventory_prices(cfg, ds.cur, ds.next);
    phase_5_accounting_investment(cfg, ds.cur, ds.next);
    phase_6_metrics_and_wage_update(cfg, ds.cur, ds.next);

    check_invariants(ds.next, t);
    writer.write_row_from_state(t, ds.next);

    ds.swap_buffers();
  }

  return hash_state_fingerprint(ds.cur);
}

static void phase_1_plan_production_mt(const SimConfig& cfg, const SimState& cur, SimState& next, ThreadPool& pool) {
  next = cur;

  std::vector<i32> firms_in_sector((std::size_t)cur.S, 0);
  for (const auto& f : cur.firms) firms_in_sector[(std::size_t)f.sector] += 1;
  for (auto& c : firms_in_sector) if (c <= 0) c = 1;

  Vec D_final((std::size_t)cur.S, 0.0);
  const f64 income = cur.agg.W;
  const f64 consumption = (1.0 - cfg.save_rate) * income;

  for (i32 j = 0; j < cur.S; ++j) {
    const f64 p = clamp_min(cur.sectors[(std::size_t)j].price, EPS);
    const f64 spend = consumption * cfg.household_c_share[(std::size_t)j];
    D_final[(std::size_t)j] = safe_div(spend, p, 0.0);
    next.sectors[(std::size_t)j].demand_final = D_final[(std::size_t)j];
  }

  pool.parallel_for(0, cur.N, [&](int b, int e, int) {
    for (int i = b; i < e; ++i) {
      auto& nf = next.firms[(std::size_t)i];
      const i32 j = nf.sector;

      const f64 cap = safe_div(cur.firms[(std::size_t)i].K, cfg.kappa_cap, 0.0);
      const f64 Dpf = safe_div(D_final[(std::size_t)j], (f64)firms_in_sector[(std::size_t)j], 0.0);
      const f64 inv_target = cfg.inv_target_theta * Dpf;

      const f64 desired = Dpf + (inv_target - cur.firms[(std::size_t)i].inv);
      nf.q_plan = std::clamp(desired, 0.0, cap);
    }
  });
}

static void phase_2_aggregate_plans_and_demands_mt(const SimConfig& cfg, const SimState& cur, SimState& next, ThreadPool& pool) {
  const int S = cur.S;
  const int T = pool.size();

  std::vector<Vec> Qp((std::size_t)T, Vec((std::size_t)S, 0.0));

  pool.parallel_for(0, cur.N, [&](int b, int e, int tid) {
    auto& local = Qp[(std::size_t)tid];
    for (int i = b; i < e; ++i) {
      const auto& f = next.firms[(std::size_t)i];
      local[(std::size_t)f.sector] += f.q_plan;
    }
  });

  Vec Q_plan((std::size_t)S, 0.0);
  for (int tid = 0; tid < T; ++tid) {
    for (int j = 0; j < S; ++j) Q_plan[(std::size_t)j] += Qp[(std::size_t)tid][(std::size_t)j];
  }

  Vec D_inter((std::size_t)S, 0.0);
  for (int i = 0; i < S; ++i) {
    f64 acc = 0.0;
    for (int j = 0; j < S; ++j) acc += cfg.A[(std::size_t)i][(std::size_t)j] * Q_plan[(std::size_t)j];
    D_inter[(std::size_t)i] = acc;
  }

  for (int i = 0; i < S; ++i) {
    const f64 df = next.sectors[(std::size_t)i].demand_final;
    next.sectors[(std::size_t)i].demand_total = df + D_inter[(std::size_t)i];
  }
}

static void phase_3_labor_and_output_mt(const SimConfig& cfg, const SimState& cur, SimState& next, ThreadPool& pool) {
  const int T = pool.size();
  std::vector<f64> Lreq((std::size_t)T, 0.0);

  pool.parallel_for(0, cur.N, [&](int b, int e, int tid) {
    f64 acc = 0.0;
    for (int i = b; i < e; ++i) {
      const i32 j = next.firms[(std::size_t)i].sector;
      acc += cfg.l[(std::size_t)j] * next.firms[(std::size_t)i].q_plan;
    }
    Lreq[(std::size_t)tid] = acc;
  });

  f64 L_req_total = 0.0;
  for (int tid = 0; tid < T; ++tid) L_req_total += Lreq[(std::size_t)tid];

  f64 phi = 1.0;
  if (L_req_total > 0.0) phi = std::min(1.0, safe_div(cfg.labor_supply, L_req_total, 1.0));

  std::vector<f64> Lemp((std::size_t)T, 0.0);
  std::vector<f64> Qsum((std::size_t)T, 0.0);

  pool.parallel_for(0, cur.N, [&](int b, int e, int tid) {
    f64 le = 0.0;
    f64 qs = 0.0;
    for (int i = b; i < e; ++i) {
      auto& nf = next.firms[(std::size_t)i];
      const i32 j = nf.sector;

      const f64 L_req = cfg.l[(std::size_t)j] * nf.q_plan;
      nf.L = phi * L_req;
      nf.q_eff = phi * nf.q_plan;

      le += nf.L;
      qs += nf.q_eff;
    }
    Lemp[(std::size_t)tid] = le;
    Qsum[(std::size_t)tid] = qs;
  });

  f64 L_emp_total = 0.0;
  f64 Q_total = 0.0;
  for (int tid = 0; tid < T; ++tid) {
    L_emp_total += Lemp[(std::size_t)tid];
    Q_total += Qsum[(std::size_t)tid];
  }

  next.agg.total_employment_hours = L_emp_total;
  next.agg.total_output = Q_total;
  next.agg.unemployment = 1.0 - safe_div(L_emp_total, cfg.labor_supply, 0.0);
  next.agg.unemployment = std::clamp(next.agg.unemployment, 0.0, 1.0);
}

static void phase_4_sales_inventory_prices_mt(const SimConfig& cfg, const SimState& cur, SimState& next, ThreadPool& pool) {
  const int S = cur.S;
  const int T = pool.size();

  std::vector<Vec> supply_loc((std::size_t)T, Vec((std::size_t)S, 0.0));
  std::vector<Vec> avail_loc((std::size_t)T, Vec((std::size_t)S, 0.0));

  pool.parallel_for(0, cur.N, [&](int b, int e, int tid) {
    auto& sl = supply_loc[(std::size_t)tid];
    auto& al = avail_loc[(std::size_t)tid];
    for (int i = b; i < e; ++i) {
      const i32 j = next.firms[(std::size_t)i].sector;
      sl[(std::size_t)j] += next.firms[(std::size_t)i].q_eff;
      al[(std::size_t)j] += cur.firms[(std::size_t)i].inv + next.firms[(std::size_t)i].q_eff;
    }
  });

  Vec supply((std::size_t)S, 0.0);
  Vec avail((std::size_t)S, 0.0);
  for (int tid = 0; tid < T; ++tid) {
    for (int j = 0; j < S; ++j) {
      supply[(std::size_t)j] += supply_loc[(std::size_t)tid][(std::size_t)j];
      avail[(std::size_t)j] += avail_loc[(std::size_t)tid][(std::size_t)j];
    }
  }

  Vec sales_sec((std::size_t)S, 0.0);
  Vec denom((std::size_t)S, 0.0);
  for (int j = 0; j < S; ++j) {
    next.sectors[(std::size_t)j].supply = supply[(std::size_t)j];

    const f64 D = next.sectors[(std::size_t)j].demand_total;
    const f64 Asec = avail[(std::size_t)j];
    const f64 ss = std::min(Asec, D);

    next.sectors[(std::size_t)j].sales = ss;
    sales_sec[(std::size_t)j] = ss;
    denom[(std::size_t)j] = std::max(Asec, EPS);

    const f64 p_cur = clamp_min(cur.sectors[(std::size_t)j].price, EPS);
    const f64 sup = supply[(std::size_t)j];
    const f64 excess = safe_div((D - sup), std::max(sup, EPS), 0.0);
    const f64 p_next = p_cur * (1.0 + cfg.price_adj * excess);
    next.sectors[(std::size_t)j].price = clamp_min(p_next, EPS);
  }

  pool.parallel_for(0, cur.N, [&](int b, int e, int) {
    for (int i = b; i < e; ++i) {
      auto& nf = next.firms[(std::size_t)i];
      const i32 j = nf.sector;

      const f64 a_i = cur.firms[(std::size_t)i].inv + nf.q_eff;
      const f64 s_i = sales_sec[(std::size_t)j] * safe_div(a_i, denom[(std::size_t)j], 0.0);

      nf.sales = s_i;
      nf.inv = std::max(0.0, a_i - s_i);
    }
  });
}

static void phase_5_accounting_investment_mt(const SimConfig& cfg, const SimState& cur, SimState& next, ThreadPool& pool) {
  const int S = cur.S;
  Vec prices((std::size_t)S, 0.0);
  for (int j = 0; j < S; ++j) prices[(std::size_t)j] = clamp_min(cur.sectors[(std::size_t)j].price, EPS);

  pool.parallel_for(0, cur.N, [&](int b, int e, int) {
    for (int i = b; i < e; ++i) {
      auto& nf = next.firms[(std::size_t)i];
      const auto& cf = cur.firms[(std::size_t)i];
      const i32 j = nf.sector;

      nf.revenue = prices[(std::size_t)j] * nf.sales;

      f64 ic = 0.0;
      for (int k = 0; k < S; ++k) ic += prices[(std::size_t)k] * cfg.A[(std::size_t)k][(std::size_t)j] * nf.q_eff;
      nf.cost_intermediate = ic;

      nf.wage_bill = cur.agg.wage * nf.L;

      nf.profit = nf.revenue - nf.cost_intermediate - nf.wage_bill;

      const f64 I = cfg.inv_rate * std::max(nf.profit, 0.0);
      const f64 K_next = (1.0 - cfg.delta) * cf.K + I;
      nf.K = std::max(0.0, K_next);
    }
  });
}

static void phase_6_metrics_and_wage_update_mt(const SimConfig& cfg, const SimState& cur, SimState& next, ThreadPool& pool) {
  const int T = pool.size();
  std::vector<f64> Wloc((std::size_t)T, 0.0);
  std::vector<f64> Cloc((std::size_t)T, 0.0);
  std::vector<f64> Kdloc((std::size_t)T, 0.0);

  pool.parallel_for(0, cur.N, [&](int b, int e, int tid) {
    f64 W = 0.0;
    f64 C = 0.0;
    f64 Kd = 0.0;
    for (int i = b; i < e; ++i) {
      W += next.firms[(std::size_t)i].wage_bill;
      C += next.firms[(std::size_t)i].cost_intermediate;
      Kd += cfg.delta * cur.firms[(std::size_t)i].K;
    }
    Wloc[(std::size_t)tid] = W;
    Cloc[(std::size_t)tid] = C;
    Kdloc[(std::size_t)tid] = Kd;
  });

  f64 W = 0.0;
  f64 C_inter = 0.0;
  f64 K_dep = 0.0;
  for (int tid = 0; tid < T; ++tid) {
    W += Wloc[(std::size_t)tid];
    C_inter += Cloc[(std::size_t)tid];
    K_dep += Kdloc[(std::size_t)tid];
  }

  const f64 V = W;
  const f64 C = C_inter + K_dep;
  const f64 L_emp = next.agg.total_employment_hours;

  f64 S_val = L_emp - W;
  if (!is_finite(S_val)) S_val = 0.0;

  next.agg.W = W;
  next.agg.S = S_val;
  next.agg.e = safe_div(S_val, std::max(W, EPS), 0.0);
  next.agg.C_over_V = safe_div(C, std::max(V, EPS), 0.0);
  next.agg.r = safe_div(S_val, std::max(C + V, EPS), 0.0);

  const f64 u = next.agg.unemployment;
  const f64 w_next = cur.agg.wage * (1.0 + cfg.wage_adj * (cfg.u_ref - u));
  next.agg.wage = clamp_min(w_next, EPS);
}

std::uint64_t run_multi_thread(const SimConfig& cfg, DoubleState& ds, const std::string& out_csv_path) {
  ThreadPool pool(cfg.threads);

  CsvWriter writer(out_csv_path, cfg.S);
  writer.write_header();

  for (i32 t = 0; t < cfg.ticks; ++t) {
    phase_1_plan_production_mt(cfg, ds.cur, ds.next, pool);
    phase_2_aggregate_plans_and_demands_mt(cfg, ds.cur, ds.next, pool);
    phase_3_labor_and_output_mt(cfg, ds.cur, ds.next, pool);
    phase_4_sales_inventory_prices_mt(cfg, ds.cur, ds.next, pool);
    phase_5_accounting_investment_mt(cfg, ds.cur, ds.next, pool);
    phase_6_metrics_and_wage_update_mt(cfg, ds.cur, ds.next, pool);

    check_invariants(ds.next, t);
    writer.write_row_from_state(t, ds.next);

    ds.swap_buffers();
  }

  return hash_state_fingerprint(ds.cur);
}

std::uint64_t run(const SimConfig& cfg, DoubleState& ds, const std::string& out_csv_path) {
  if (cfg.threads <= 1) return run_single_thread(cfg, ds, out_csv_path);
  return run_multi_thread(cfg, ds, out_csv_path);
}

}