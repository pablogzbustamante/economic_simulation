#pragma once
#include <string>
#include "sim/types.hpp"

namespace sim {

    struct SimConfig{

        u64 seed;
        i32 ticks;
        i32 threads;
        i32 N_firms;
        i32 S;
        
        Mat A; 
        Vec l;
        f64 delta;
        f64 kappa_cap;
        
        Vec init_price;
        f64 price_adj;
        Vec household_c_share;
        f64 save_rate;
        
        f64 labor_supply;
        f64 init_wage;
        f64 wage_adj;
        f64 u_ref;
        
        f64 inv_rate;
        Vec markup;
        
        f64 inv_target_theta;
    };
        
    SimConfig load_config(const std::string& path);
    void validate_config(const SimConfig& cfg);
}
