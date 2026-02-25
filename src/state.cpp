#include "sim/state.hpp"

namespace sim {

SimState make_state(i32 S, i32 N) {
  SimState st;
  st.S = S;
  st.N = N;
  st.firms.resize((std::size_t)N);
  st.sectors.resize((std::size_t)S);
  st.agg = AggState{};
  return st;
}

}