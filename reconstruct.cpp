#include "reconstruct.hpp"
#include <cstring>

void reconstruct_const(const GriMesh& mesh, const double* U,
                      double* UL_int, double* UR_int,
                      double* UL_bnd, double* UR_bnd) {
    for (int i = 0; i < mesh.num_interior_faces; ++i) {
        int elemL = mesh.I2E[4 * i + 0];
        int elemR = mesh.I2E[4 * i + 2];
        std::memcpy(UL_int + i * 4, U + elemL * 4, 4 * sizeof(double));
        std::memcpy(UR_int + i * 4, U + elemR * 4, 4 * sizeof(double));
    }
    for (int i = 0; i < mesh.num_boundary_faces; ++i) {
        int elem = mesh.B2E[3 * i + 0];
        std::memcpy(UL_bnd + i * 4, U + elem * 4, 4 * sizeof(double));
        std::memcpy(UR_bnd + i * 4, U + elem * 4, 4 * sizeof(double));
    }
}
