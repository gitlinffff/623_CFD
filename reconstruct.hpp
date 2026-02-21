#ifndef RECONSTRUCT_HPP
#define RECONSTRUCT_HPP

#include "readgri.hpp"

/** Reconstruction function: recon(mesh, U, UL_int, UR_int, UL_bnd, UR_bnd, gamma) */
typedef void (*ReconFn)(const GriMesh& mesh, const double* U,
                       double* UL_int, double* UR_int,
                       double* UL_bnd, double* UR_bnd, double gamma);

/**
 * No reconstruction (constant / Godunov): UL = U[elemL], UR = U[elemR].
 * For boundary faces: UL = U[elem], UR = U[elem] (caller overwrites UR with ghost before flux).
 */
void reconstruct_const(const GriMesh& mesh, const double* U,
                      double* UL_int, double* UR_int,
                      double* UL_bnd, double* UR_bnd, double gamma);

/**
 * Second-order linear reconstruction without limiter.
 * U_face = U_cell + gradU · (x_face - x_cell). Falls back to cell average if invalid.
 */
void reconstruct_nolimiter(const GriMesh& mesh, const double* U,
                           double* UL_int, double* UR_int,
                           double* UL_bnd, double* UR_bnd, double gamma);

/** Clip U to valid physical bounds (rho>0, p>0, finite). Call after time update to prevent blow-up. */
void clip_cons_state(double U[4], double gamma);

#endif
