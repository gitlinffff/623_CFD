#ifndef RECONSTRUCT_HPP
#define RECONSTRUCT_HPP

#include "readgri.hpp"

/** Reconstruction function: recon(mesh, U, UL_int, UR_int, UL_bnd, UR_bnd) */
typedef void (*ReconFn)(const GriMesh& mesh, const double* U,
                       double* UL_int, double* UR_int,
                       double* UL_bnd, double* UR_bnd);

/**
 * No reconstruction (constant / Godunov): UL = U[elemL], UR = U[elemR].
 * For boundary faces: UL = U[elem], UR = U[elem] (caller overwrites UR with ghost before flux).
 *
 * U: cell-averaged state, size Ne*4 [rho, rho*u, rho*v, E] per element
 * UL_int, UR_int: output for interior faces, size num_interior_faces * 4 each
 * UL_bnd, UR_bnd: output for boundary faces, size num_boundary_faces * 4 each
 */
void reconstruct_const(const GriMesh& mesh, const double* U,
                      double* UL_int, double* UR_int,
                      double* UL_bnd, double* UR_bnd);

#endif
