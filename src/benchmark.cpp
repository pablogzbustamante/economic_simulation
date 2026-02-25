#include "sim/benchmark.hpp"
#include <cstdint>
#include <cstring>

namespace sim {

static std::uint64_t fnv1a_u64(const void* data, std::size_t n) {
  const std::uint8_t* p = (const std::uint8_t*)data;
  std::uint64_t h = 1469598103934665603ull;
  for (std::size_t i = 0; i < n; ++i) {
    h ^= (std::uint64_t)p[i];
    h *= 1099511628211ull;
  }
  return h;
}

static void hash_double(std::uint64_t& h, double x) {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &x, sizeof(bits));
  h ^= fnv1a_u64(&bits, sizeof(bits));
  h *= 1099511628211ull;
}

std::uint64_t hash_state_fingerprint(const SimState& st) {
  std::uint64_t h = 1469598103934665603ull;

  hash_double(h, st.agg.wage);
  hash_double(h, st.agg.unemployment);
  hash_double(h, st.agg.total_output);
  hash_double(h, st.agg.total_employment_hours);
  hash_double(h, st.agg.W);
  hash_double(h, st.agg.S);
  hash_double(h, st.agg.e);
  hash_double(h, st.agg.C_over_V);
  hash_double(h, st.agg.r);

  for (int j = 0; j < st.S; ++j) {
    const auto& s = st.sectors[(std::size_t)j];
    hash_double(h, s.price);
    hash_double(h, s.demand_final);
    hash_double(h, s.demand_total);
    hash_double(h, s.supply);
    hash_double(h, s.sales);
    hash_double(h, s.value_unit);
  }

  const int stride = st.N >= 64 ? (st.N / 64) : 1;
  for (int i = 0; i < st.N; i += stride) {
    const auto& f = st.firms[(std::size_t)i];
    hash_double(h, f.K);
    hash_double(h, f.inv);
    hash_double(h, f.L);
    hash_double(h, f.q_plan);
    hash_double(h, f.q_eff);
    hash_double(h, f.sales);
    hash_double(h, f.revenue);
    hash_double(h, f.cost_intermediate);
    hash_double(h, f.wage_bill);
    hash_double(h, f.profit);
  }

  return h;
}

}