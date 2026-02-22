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

    GriMesh mesh;
    if (!read_gri("mesh/global_refine_2.gri", mesh)) {
        std::cerr << "Failed to read mesh.\n";
        return 1;
    }

    const double gammad = 1.4;
    FluxFn flux_fn = fluxROE;
    ReconFn recon_fn = reconstruct_nolimiter;

    const char* rst_file = "/home/linfel/umich_course/623_CFD/data/unsteady/order1_unsteady_globalrefine2_solutions/steady_solution.vtu";
    const char* out_dir = "/home/linfel/umich_course/623_CFD/data/unsteady/order1_unsteady_globalrefine2_solutions";
    const double t_end = 200;         /* run until periodic; adjust as needed */
    const double vtu_interval = 0.2;	
    const double CFL = 0.3;

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

    //solve_steady(mesh, U.data(), gammad, params, flux_fn, recon_fn, 0.3, 50, 1000000);
    solve_unsteady(mesh, U.data(), gammad, params, flux_fn, recon_fn, CFL, t_end, vtu_interval, 50, out_dir);

    return 0;
}

int run_unsteady() {
    GriMesh mesh;
    if (!read_gri("mesh/global_refine_2.gri", mesh)) {
        std::cerr << "Failed to read mesh.\n";
        return 1;
    }

    const double gammad = 1.4;
    FluxFn flux_fn = fluxROE;
    ReconFn recon_fn = reconstruct_nolimiter;

    const char* out_dir = "/home/linfel/umich_course/623_CFD/data/unsteady/order1_unsteady_globalrefine2_solutions";
    const double t_end = 200;         /* run until periodic; adjust as needed */
    const double vtu_interval = 0.2;	
    const double CFL = 0.3;

    ProblemParams params;
    std::vector<double> U(mesh.Ne * 4);

    std::cout << "Mesh: " << mesh.Ne << " elements, "
              << mesh.num_interior_faces << " interior, "
              << mesh.num_boundary_faces << " boundary faces.\n";

    ensure_dir("data");
    ensure_dir("data/results");

    initialize_uniform(U.data(), mesh.Ne, 0.1, params);

    solve_steady(mesh, U.data(), gammad, params, flux_fn, recon_fn, CFL, 50, 1000000);

    std::string filename = std::string(out_dir) + "/steady_solution.vtu";
    if (write_vtu(mesh, U.data(), gammad, filename.c_str()))
        std::cout << "Output: " << filename << "\n";
    else
        std::cerr << "Failed to write VTU.\n";

    solve_unsteady(mesh, U.data(), gammad, params, flux_fn, recon_fn, CFL, t_end, vtu_interval, 50, out_dir);
    
		return 0;
}

int main() {
    rst_unsteady();	
}
