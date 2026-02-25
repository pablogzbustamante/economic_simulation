#pragma once
#include <string>
#include <cstdint>
#include "sim/config.hpp"
#include "sim/state.hpp"

namespace sim {

struct DoubleState {
  SimState cur;
  SimState next;

  DoubleState(i32 S, i32 N);
  void swap_buffers();
};

DoubleState init_sim(const SimConfig& cfg);

    std::uint64_t run_single_thread(const SimConfig& cfg, DoubleState& ds, const std::string& out_csv_path);
    std::uint64_t run_multi_thread(const SimConfig& cfg, DoubleState& ds, const std::string& out_csv_path);
    std::uint64_t run(const SimConfig& cfg, DoubleState& ds, const std::string& out_csv_path);

}