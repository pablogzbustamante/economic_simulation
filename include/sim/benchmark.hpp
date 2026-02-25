#pragma once 
#include <cstdint>
#include "sim/state.hpp"

namespace sim {
    std::uint64_t hash_state_fingerprint(const SimState& st);
}

