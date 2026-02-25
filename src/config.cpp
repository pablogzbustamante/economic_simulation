#include "sim/config.hpp"
#include "sim/util.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <cmath>

namespace sim {

static nlohmann::json read_json_file(const std::string& path) {
  std::ifstream in(path);
  if (!in.is_open()) die("cannot open config: " + path);
  nlohmann::json j;
  in >> j;
  return j;
}

static Vec read_vec(const nlohmann::json& j, const char* key) {
  if (!j.contains(key) || !j.at(key).is_array()) die(std::string("missing/invalid array: ") + key);
  Vec v;
  for (const auto& x : j.at(key)) v.push_back(x.get<f64>());
  return v;
}

static Mat read_mat(const nlohmann::json& j, const char* key) {
  if (!j.contains(key) || !j.at(key).is_array()) die(std::string("missing/invalid matrix: ") + key);
  Mat m;
  for (const auto& row : j.at(key)) {
    if (!row.is_array()) die(std::string("invalid matrix row: ") + key);
    Vec r;
    for (const auto& x : row) r.push_back(x.get<f64>());
    m.push_back(std::move(r));
  }
  return m;
}

SimConfig load_config(const std::string& path) {
  auto j = read_json_file(path);

  SimConfig cfg{};
  cfg.seed = j.at("seed").get<u64>();
  cfg.ticks = j.at("ticks").get<i32>();
  cfg.threads = j.at("threads").get<i32>();
  cfg.N_firms = j.at("N_firms").get<i32>();
  cfg.S = j.at("S").get<i32>();

  cfg.A = read_mat(j, "A");
  cfg.l = read_vec(j, "l");
  cfg.delta = j.at("delta").get<f64>();
  cfg.kappa_cap = j.at("kappa_cap").get<f64>();

  cfg.init_price = read_vec(j, "init_price");
  cfg.price_adj = j.at("price_adj").get<f64>();
  cfg.household_c_share = read_vec(j, "household_c_share");
  cfg.save_rate = j.at("save_rate").get<f64>();

  cfg.labor_supply = j.at("labor_supply").get<f64>();
  cfg.init_wage = j.at("init_wage").get<f64>();
  cfg.wage_adj = j.at("wage_adj").get<f64>();
  cfg.u_ref = j.at("u_ref").get<f64>();

  cfg.inv_rate = j.at("inv_rate").get<f64>();
  cfg.markup = read_vec(j, "markup");

  cfg.inv_target_theta = j.value("inv_target_theta", 0.2);

  return cfg;
}

void validate_config(const SimConfig& cfg) {
  if (cfg.S < 2 || cfg.S > 5) die("S must be in [2,5]");
  if (cfg.ticks <= 0) die("ticks must be > 0");
  if (cfg.threads <= 0) die("threads must be > 0");
  if (cfg.N_firms <= 0) die("N_firms must be > 0");

  if ((i32)cfg.A.size() != cfg.S) die("A rows != S");
  for (const auto& row : cfg.A) if ((i32)row.size() != cfg.S) die("A cols != S");

  if ((i32)cfg.l.size() != cfg.S) die("l size != S");
  if ((i32)cfg.init_price.size() != cfg.S) die("init_price size != S");
  if ((i32)cfg.household_c_share.size() != cfg.S) die("household_c_share size != S");
  if ((i32)cfg.markup.size() != cfg.S) die("markup size != S");

  for (i32 j = 0; j < cfg.S; ++j) {
    if (!(cfg.l[j] >= 0.0)) die("l has negative");
    if (!(cfg.init_price[j] > 0.0)) die("init_price must be > 0");
    if (!is_finite(cfg.init_price[j])) die("init_price has NaN/Inf");
  }

  if (!(cfg.delta >= 0.0 && cfg.delta <= 1.0)) die("delta must be in [0,1]");
  if (!(cfg.kappa_cap > 0.0)) die("kappa_cap must be > 0");

  if (!(cfg.price_adj >= 0.0)) die("price_adj must be >= 0");
  if (!(cfg.save_rate >= 0.0 && cfg.save_rate <= 1.0)) die("save_rate must be in [0,1]");

  if (!(cfg.labor_supply > 0.0)) die("labor_supply must be > 0");
  if (!(cfg.init_wage > 0.0)) die("init_wage must be > 0");
  if (!(cfg.wage_adj >= 0.0)) die("wage_adj must be >= 0");
  if (!(cfg.u_ref >= 0.0 && cfg.u_ref <= 1.0)) die("u_ref must be in [0,1]");

  if (!(cfg.inv_rate >= 0.0 && cfg.inv_rate <= 1.0)) die("inv_rate must be in [0,1]");

  f64 share_sum = 0.0;
  for (auto s : cfg.household_c_share) share_sum += s;
  if (std::abs(share_sum - 1.0) > SHARE_TOL) die("household_c_share must sum to 1");
}

}