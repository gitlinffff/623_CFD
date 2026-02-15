/**
 * Euler finite-volume solver - main driver.
 * Build: mkdir build && cd build && cmake .. && make
 */
#include "solver.hpp"
#include <iostream>
#include <vector>

int main() {
    GriMesh mesh;
    if (!read_gri("mesh/initial_mesh_3.gri", mesh)) {
        std::cerr << "Failed to read mesh.\n";
        return 1;
    }

    const double gamma = 1.4;
    ProblemParams params;
    std::vector<double> U(mesh.Ne * 4);
    initialize_uniform(U.data(), mesh.Ne, 0.1, params);

    std::cout << "Mesh: " << mesh.Ne << " elements, "
              << mesh.num_interior_faces << " interior, "
              << mesh.num_boundary_faces << " boundary faces.\n";

    solve_steady(mesh, U.data(), gamma, params, fluxROE, reconstruct_const, 0.3, 50);

    std::cout << "Done. Cell 0: rho=" << U[0] << " rhoU=" << U[1]
              << " rhoV=" << U[2] << " E=" << U[3] << "\n";

    if (write_vtu(mesh, U.data(), gamma, "data/solution.vtu"))
        std::cout << "Exported: data/solution.vtu\n";
    else
        std::cerr << "Failed to export VTU (ensure data/ exists).\n";

    return 0;
}
