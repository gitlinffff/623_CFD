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

int rst_unsteady() {
    /* restart unsteady run*/
    const double gammad = 1.4;
    FluxFn flux_fn = fluxROE;
    ReconFn recon_fn = reconstruct_nolimiter;

    const double t_end = 400;         /* run until periodic; adjust as needed */
    const double vtu_interval = 0.2;  
    const double CFL = 0.3;
    const char* gri_file = "/home/linfel/umich_course/623_CFD/mesh/ver2/coarse_mesh.gri";
    const char* rst_file = "/home/linfel/umich_course/623_CFD/data/steady_results/coarse_2nd.vtu";
    const char* out_dir = "/home/linfel/umich_course/623_CFD/data/unsteady/2nd_nolimit_coarse_solutions";

    GriMesh mesh;
    if (!read_gri(gri_file, mesh)) {
        std::cerr << "Failed to read mesh.\n";
        return 1;
    }
    ensure_dir(out_dir);

    ProblemParams params;
    std::vector<double> U(mesh.Ne * 4);

    std::cout << "Mesh: " << mesh.Ne << " elements, "
              << mesh.num_interior_faces << " interior, "
              << mesh.num_boundary_faces << " boundary faces.\n";

    if (!read_vtu(mesh, rst_file, gammad, U.data())) {
        std::cerr << "Error: failed to read or mesh mismatch: " << rst_file << "\n";
        return 1;
    }
    std::cout << "Loaded: " << rst_file << "\n";

    solve_unsteady(mesh, U.data(), gammad, params, flux_fn, recon_fn, CFL, t_end, vtu_interval, 50, out_dir);

    return 0;
}

int main() {
    rst_unsteady();
}
