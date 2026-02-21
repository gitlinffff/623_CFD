/**
 * Euler finite-volume solver - unified interface.
 *
 * Usage: #include "solver.hpp"
 *
 * Core API:
 *   - GriMesh, read_gri()         mesh loading
 *   - ProblemParams, initialize_uniform()  problem init
 *   - solve_steady()               steady-state solve (SSP-RK3 + local time stepping)
 *   - fluxROE / fluxHLLC / fluxRusanov  numerical flux
 */
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "readgri.hpp"
#include "problem.hpp"
#include "reconstruct.hpp"
#include "advance.hpp"
#include "flux.hpp"
#include "write_vtu.hpp"

#endif
