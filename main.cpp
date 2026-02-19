/**
 * Euler finite-volume solver - main driver.
 * Build: mkdir build && cd build && cmake .. && make
 *
 * Usage: ./main [steady|unsteady]
 *   steady   - steady-state simulation (default), output: data/solution.vtu, data/residual_history.dat
 *   unsteady - unsteady with stator wake inflow, output: data/unsteady/current/
 */
#include "solver.hpp"
#include <iostream>
#include <vector>
#include <cstring>
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

static void run_steady(const GriMesh& mesh, double* U, double gamma,
                       ProblemParams& params) {
    ensure_dir("data");
    std::cout << "--- Steady-state ---\n";
    solve_steady(mesh, U, gamma, params, fluxROE, reconstruct_const, 0.3, 50);

    std::cout << "Done. Cell 0: rho=" << U[0] << " rhoU=" << U[1]
              << " rhoV=" << U[2] << " E=" << U[3] << "\n";

    if (write_vtu(mesh, U, gamma, "data/solution.vtu"))
        std::cout << "Exported: data/solution.vtu\n";
    else
        std::cerr << "Failed to export VTU (ensure data/ exists).\n";
}

static void run_unsteady(const GriMesh& mesh, double* U, double gamma,
                        ProblemParams& params) {
    const char* out_dir = "data/unsteady/current";
    const double t_end = 200;         /* run until periodic; adjust as needed */
    const double vtu_interval = 0.2;
    const double CFL = 0.5;

    ensure_dir("data");
    ensure_dir("data/unsteady");
    ensure_dir(out_dir);

    std::cout << "--- Steady-state (initial condition) ---\n";
    solve_steady(mesh, U, gamma, params, fluxROE, reconstruct_const, 0.3, 50);

    std::cout << "--- Unsteady (t_end=" << t_end << ", Vrot=" << params.Vrot << ") ---\n";
    solve_unsteady(mesh, U, gamma, params, fluxROE, reconstruct_const,
                   CFL, t_end, vtu_interval, 50, out_dir);
}

int main(int argc, char* argv[]) {
    GriMesh mesh;
    if (!read_gri("mesh/global_refine_1.gri", mesh)) {
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

    bool unsteady = true;
    if (argc >= 2 && std::strcmp(argv[1], "unsteady") == 0)
        unsteady = true;

    if (unsteady)
        run_unsteady(mesh, U.data(), gamma, params);
    else
        run_steady(mesh, U.data(), gamma, params);

    return 0;
}
