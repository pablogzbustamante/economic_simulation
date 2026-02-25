# Economic Simulation (Marxist Multi-Sector Model)

Discrete-time multi-sector economic simulation in C++20.

This project implements a tick-based simulation of a simplified
Marxist production economy with:

- Input–Output technology matrix (A)
- Direct labor coefficients (l)
- Fixed capital with depreciation
- Wage adjustment (Goodwin-style rule)
- Price adjustment by excess demand
- Firm-level accounting
- Double-buffered state engine
- Deterministic seed support

Current Status:
Bootstrap complete.
Config parser + validation implemented.
Core engine under development.

---

## Build Instructions

### Requirements
- CMake >= 3.20
- C++20 compiler (MSVC / Clang / GCC)

### Build (Windows)

cmake -S . -B build
cmake --build build --config Release

Run:

.\build\Release\economic_simulation.exe

### Build (Linux/macOS)

cmake -S . -B build
cmake --build build -j

Run:

./build/economic_simulation

Author: Pablo González Bustamante