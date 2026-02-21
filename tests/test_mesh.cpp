#include "readgri.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>

static void verify_mesh(const char* fname) {
    std::printf("Verifying mesh: %s\n", fname);
    GriMesh mesh;
    if (!read_gri(fname, mesh)) {
        std::printf("Failed to read %s\n", fname);
        std::printf("------------------------------\n");
        return;
    }

    int Ne = mesh.Ne;
    std::printf("Number of elements: %d\n", Ne);

    std::vector<double> element_errors(Ne * 2, 0.0);

    for (int i = 0; i < mesh.num_interior_faces; ++i) {
        int elemL = mesh.I2E[4 * i + 0];
        int elemR = mesh.I2E[4 * i + 2];
        double nx = mesh.In[2 * i + 0];
        double ny = mesh.In[2 * i + 1];
        double length = mesh.In_len[i];
        element_errors[elemL * 2 + 0] += nx * length;
        element_errors[elemL * 2 + 1] += ny * length;
        element_errors[elemR * 2 + 0] -= nx * length;
        element_errors[elemR * 2 + 1] -= ny * length;
    }

    for (int i = 0; i < mesh.num_boundary_faces; ++i) {
        int elem = mesh.B2E[3 * i + 0];
        double nx = mesh.Bn[2 * i + 0];
        double ny = mesh.Bn[2 * i + 1];
        double length = mesh.Bn_len[i];
        element_errors[elem * 2 + 0] += nx * length;
        element_errors[elem * 2 + 1] += ny * length;
    }

    double max_error = 0.0;
    int error_elements = 0;
    const double tol = 1e-12;
    for (int e = 0; e < Ne; ++e) {
        double ex = element_errors[e * 2 + 0];
        double ey = element_errors[e * 2 + 1];
        double mag = std::sqrt(ex * ex + ey * ey);
        if (mag > max_error) max_error = mag;
        if (mag > tol) error_elements++;
    }

    std::printf("Maximum error magnitude: %.6e\n", max_error);
    std::printf("Number of elements with errors: %d\n", error_elements);

    if (max_error < tol) {
        std::printf("Verification PASSED: Errors are within machine precision.\n");
    } else {
        std::printf("Verification FAILED: Errors exceed machine precision.\n");
    }
    std::printf("------------------------------\n");
}

int main() {
    verify_mesh("mesh/initial_mesh_3.gri");
    return 0;
}
