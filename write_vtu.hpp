#ifndef WRITE_VTU_HPP
#define WRITE_VTU_HPP

#include "readgri.hpp"

/**
 * Export cell-centered flow field to Paraview-readable VTU (UnstructuredGrid).
 * Converts conserved U to primitive (rho, u, v, p) and writes to file.
 *
 * filepath: e.g. "data/solution.vtu" or "data/t=0.123.vtu"
 * Returns true on success.
 */
bool write_vtu(const GriMesh& mesh, const double* U, double gamma,
              const char* filepath);

/**
 * Read cell-centered flow field from a VTU file written by write_vtu().
 * Parses primitive variables (rho, u, v, p) and converts back to conserved U.
 * U must be pre-allocated with size mesh.Ne * 4.
 * Returns true on success; on failure U is unchanged.
 */
bool read_vtu(const GriMesh& mesh, const char* filepath, double gamma, double* U);

#endif
