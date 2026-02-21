/**
 * Euler finite-volume solver - steady-state only.
 * Build: mkdir build && cd build && cmake .. && make
 * Usage: ./main
 * Output: data/results/solution.vtu
 */
#include "solver.hpp"
#include <iostream>
#include <vector>
#include <cstdio>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

static void ensure_dir(const char* path) {
#ifdef _WIN32
    (void)_mkdir(path);
#else
    (void)mkdir(path, 0755);
#endif
}

int main() {
    GriMesh mesh;
    if (!read_gri("mesh/global_refine_2.gri", mesh)) {
        std::cerr << "Failed to read mesh.\n";
        return 1;
    }

    const double gamma = 1.4;
    FluxFn flux_fn = fluxROE;
    ReconFn recon_fn = reconstruct_nolimiter;

    ProblemParams params;
    std::vector<double> U(mesh.Ne * 4);

    std::cout << "Mesh: " << mesh.Ne << " elements, "
              << mesh.num_interior_faces << " interior, "
              << mesh.num_boundary_faces << " boundary faces.\n";

    ensure_dir("data");
    ensure_dir("data/results");

    std::cout << "Restart from VTU? (y/n): ";
    char choice;
    std::cin >> choice;
    if (choice == 'y' || choice == 'Y') {
        std::string restart_path;
        std::cout << "VTU file path: ";
        std::cin >> restart_path;
        if (!read_vtu(mesh, restart_path.c_str(), gamma, U.data())) {
            std::cerr << "Error: failed to read or mesh mismatch: " << restart_path << "\n";
            return 1;
        }
        std::cout << "Loaded: " << restart_path << "\n";
    } else {
        initialize_uniform(U.data(), mesh.Ne, 0.1, params);
    }

    solve_steady(mesh, U.data(), gamma, params, flux_fn, recon_fn, 0.3, 50, 1000000);

    const char* out_path = "data/results/refine2_2nd.vtu";
    if (write_vtu(mesh, U.data(), gamma, out_path))
        std::cout << "Output: " << out_path << "\n";
    else
        std::cerr << "Failed to write VTU.\n";

    return 0;
}
