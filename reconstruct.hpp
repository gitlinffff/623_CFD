#ifndef RECONSTRUCT_HPP
#define RECONSTRUCT_HPP

#include "readgri.hpp"

/** Reconstruction function: recon(mesh, U, UL_int, UR_int, UL_bnd, UR_bnd, gammad) */
typedef void (*ReconFn)(const GriMesh& mesh, const double* U,
                       double* UL_int, double* UR_int,
                       double* UL_bnd, double* UR_bnd, double gammad);

/**
 * No reconstruction (constant / Godunov): UL = U[elemL], UR = U[elemR].
 * For boundary faces: UL = U[elem], UR = U[elem] (caller overwrites UR with ghost before flux).
 */
void reconstruct_const(const GriMesh& mesh, const double* U,
                      double* UL_int, double* UR_int,
                      double* UL_bnd, double* UR_bnd, double gammad);

/**
 * Second-order linear reconstruction without limiter.
 * U_face = U_cell + gradU · (x_face - x_cell). Falls back to cell average if invalid.
 */
void reconstruct_nolimiter(const GriMesh& mesh, const double* U,
                           double* UL_int, double* UR_int,
                           double* UL_bnd, double* UR_bnd, double gammad);

/**
 * Second-order linear reconstruction with Barth-Jespersen limiter.
 * Same as nolimiter but gradients are limited so that extrapolated values lie within
 * [U_min, U_max] over the cell and its face-neighbors (monotonicity).
 */
void reconstruct_BJ(const GriMesh& mesh, const double* U,
                    double* UL_int, double* UR_int,
                    double* UL_bnd, double* UR_bnd, double gammad);

/**
 * Second-order linear reconstruction with maximum-principle (component-wise) region limiter.
 * Same as nolimiter but each conserved variable is limited so that its face values lie within
 * [U_min_k, U_max_k] over the cell and face-neighbors. Phi is per variable (less restrictive
 * than Barth-Jespersen's single phi per cell), preserving the maximum principle per component.
 */
void reconstruct_MP(const GriMesh& mesh, const double* U,
                    double* UL_int, double* UR_int,
                    double* UL_bnd, double* UR_bnd, double gammad);

/**
 * Second-order linear reconstruction with LCD (Limited Central Difference) limiter.
 * Uses one unlimited gradient L̃; for each edge k computes α_k so that edge-midpoint
 * values lie in [min(u_k, u₀), max(u_k, u₀)]; then L_LCD = (min_k α_k) * L̃.
 * Gradient direction preserved, only magnitude reduced. Robust and cheap, but diffusive.
 */
void reconstruct_LCD(const GriMesh& mesh, const double* U,
                     double* UL_int, double* UR_int,
                     double* UL_bnd, double* UR_bnd, double gammad);

/** Clip U to valid physical bounds (rho>0, p>0, finite). Call after time update to prevent blow-up. */
void clip_cons_state(double U[4], double gammad);

#endif
