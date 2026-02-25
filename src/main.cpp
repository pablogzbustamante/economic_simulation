#include <iostream>
#include <chrono>
#include "sim/config.hpp"
#include "sim/engine.hpp"

int main() {
  auto cfg = sim::load_config("config/base.json");
  sim::validate_config(cfg);

  auto ds = sim::init_sim(cfg);

  const auto t0 = std::chrono::steady_clock::now();
  auto h = sim::run(cfg, ds, "output/out.csv");
  const auto t1 = std::chrono::steady_clock::now();

  const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  std::cout << "run_ok ticks=" << cfg.ticks << " threads=" << cfg.threads << " ms=" << ms << " hash=" << h << "\n";
  return 0;
}