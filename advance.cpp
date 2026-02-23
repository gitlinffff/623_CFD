/**
 * Time advance and residual for finite-volume Euler solver.
 *
 * Periodic boundaries: readgri merges periodic pairs (BGroup1-7, 3-5) into interior faces (I2E).
 * B2E contains only BGroup2,4,6,8 (wall, outflow, inflow). No special handling for periodic BC.
 *
 * Local time stepping: dt_i = (2*A_i*CFL) / (sum_e |s|_{i,e} * L_{i,e}) per cell.
 */
#include "advance.hpp"
#include "flux.hpp"
#include "physics.hpp"
#include "reconstruct.hpp"
#include "bc.hpp"
#include "write_vtu.hpp"
#include <cmath>
#include <cstdio>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

void calcRes(const GriMesh& mesh, const double* U, double* R, double gammad,
             const ProblemParams& params, FluxFn flux_fn, ReconFn recon_fn,
             double* dt_per_cell, double CFL, double t) {
    std::memset(R, 0, mesh.Ne * 4 * sizeof(double));
    std::vector<double> sum_s(mesh.Ne, 0.0);

    std::vector<double> UL_int(mesh.num_interior_faces * 4);
    std::vector<double> UR_int(mesh.num_interior_faces * 4);
    std::vector<double> UL_bnd(mesh.num_boundary_faces * 4);
    std::vector<double> UR_bnd(mesh.num_boundary_faces * 4);

    recon_fn(mesh, U, UL_int.data(), UR_int.data(), UL_bnd.data(), UR_bnd.data(), gammad);

    double Fhat[4];
    double smag;
    const double nin[2] = {std::cos(params.alpha), std::sin(params.alpha)};
    const double R_gas = 1.0 / gammad;

    for (int i = 0; i < mesh.num_interior_faces; ++i) {
        flux_fn(&UL_int[i * 4], &UR_int[i * 4], &mesh.In[2 * i], gammad, Fhat, smag);
        double L = mesh.In_len[i];
        int elemL = mesh.I2E[4 * i + 0];
        int elemR = mesh.I2E[4 * i + 2];
        double AL = mesh.Area[elemL];
        double AR = mesh.Area[elemR];

        R[elemL * 4 + 0] -= Fhat[0] * L / AL;
        R[elemL * 4 + 1] -= Fhat[1] * L / AL;
        R[elemL * 4 + 2] -= Fhat[2] * L / AL;
        R[elemL * 4 + 3] -= Fhat[3] * L / AL;
        R[elemR * 4 + 0] += Fhat[0] * L / AR;
        R[elemR * 4 + 1] += Fhat[1] * L / AR;
        R[elemR * 4 + 2] += Fhat[2] * L / AR;
        R[elemR * 4 + 3] += Fhat[3] * L / AR;

        if (dt_per_cell) {
            sum_s[elemL] += smag * L;
            sum_s[elemR] += smag * L;
        }
    }

    for (int i = 0; i < mesh.num_boundary_faces; ++i) {
        const double* UL = &UL_bnd[i * 4];
        const double* n = &mesh.Bn[2 * i];
        double L = mesh.Bn_len[i];
        int elem = mesh.B2E[3 * i + 0];
        double A = mesh.Area[elem];

        BCType bt = get_bc_type(mesh, i);
        if (bt == BC_WALL)
            WallFlux(UL, n, gammad, Fhat, smag);
        else if (bt == BC_INFLOW) {
            double rho0_in = params.rho0;
            if (t >= 0.0) {
                /* Unsteady: rho0(eta) = rho0*[1 - fwake*exp(-eta^2/(2*delta^2))] */
                int elem = mesh.B2E[3 * i + 0];
                int face = mesh.B2E[3 * i + 1];
                int v0 = (face + 1) % 3, v1 = (face + 2) % 3;
                int vid0 = mesh.E[elem * 3 + v0], vid1 = mesh.E[elem * 3 + v1];
                double y_rot = 0.5 * (mesh.V[vid0 * 2 + 1] + mesh.V[vid1 * 2 + 1]);
                double ystator = y_rot + params.Vrot * t;
                double frac = ystator / params.delta_y - std::floor(ystator / params.delta_y);
                double eta = frac - 0.5;
                double fac = 1.0 - params.fwake * std::exp(-eta * eta / (2.0 * params.delta_wake * params.delta_wake));
                rho0_in = params.rho0 * fac;
            }
            try {
                InflowFlux(UL, n, nin, rho0_in, params.a0, gammad, R_gas, flux_fn, Fhat, smag);
            } catch (const std::runtime_error&) {
                flux_fn(UL, UL, n, gammad, Fhat, smag);
            }
        } else if (bt == BC_OUTFLOW)
            OutflowFlux(UL, n, params.pout, gammad, flux_fn, Fhat, smag);
        else
            flux_fn(UL, UL, n, gammad, Fhat, smag);

        R[elem * 4 + 0] -= Fhat[0] * L / A;
        R[elem * 4 + 1] -= Fhat[1] * L / A;
        R[elem * 4 + 2] -= Fhat[2] * L / A;
        R[elem * 4 + 3] -= Fhat[3] * L / A;

        if (dt_per_cell)
            sum_s[elem] += smag * L;
    }

    if (dt_per_cell) {
        for (int i = 0; i < mesh.Ne; ++i) {
            double denom = sum_s[i];
            if (denom < 1e-14) denom = 1e-14;
            dt_per_cell[i] = (2.0 * mesh.Area[i] * CFL) / denom;
        }
    }
}

double SSPRK3(const GriMesh& mesh, double* U, double gammad,
              const ProblemParams& params, FluxFn flux_fn, ReconFn recon_fn,
              double CFL, double t) {
    std::vector<double> U1(mesh.Ne * 4);
    std::vector<double> U2(mesh.Ne * 4);
    std::vector<double> R(mesh.Ne * 4);
    double dt = compute_dt(mesh, U, gammad, CFL);  /* global dt, same for all cells */

    calcRes(mesh, U, R.data(), gammad, params, flux_fn, recon_fn, nullptr, CFL, t);
    for (int i = 0; i < mesh.Ne; ++i) {
        for (int k = 0; k < 4; ++k)
            U1[i * 4 + k] = U[i * 4 + k] + dt * R[i * 4 + k];
        clip_cons_state(&U1[i * 4], gammad);
    }

    calcRes(mesh, U1.data(), R.data(), gammad, params, flux_fn, recon_fn, nullptr, CFL, t + dt);
    for (int i = 0; i < mesh.Ne; ++i) {
        for (int k = 0; k < 4; ++k)
            U2[i * 4 + k] = 0.75 * U[i * 4 + k] + 0.25 * (U1[i * 4 + k] + dt * R[i * 4 + k]);
        clip_cons_state(&U2[i * 4], gammad);
    }

    calcRes(mesh, U2.data(), R.data(), gammad, params, flux_fn, recon_fn, nullptr, CFL, t + dt);
    for (int i = 0; i < mesh.Ne; ++i) {
        for (int k = 0; k < 4; ++k)
            U[i * 4 + k] = (1.0 / 3.0) * U[i * 4 + k] + (2.0 / 3.0) * (U2[i * 4 + k] + dt * R[i * 4 + k]);
        clip_cons_state(&U[i * 4], gammad);
    }
    return dt;
}

void SSPRK3_local(const GriMesh& mesh, double* U, double gammad,
                  const ProblemParams& params, FluxFn flux_fn, ReconFn recon_fn, double CFL) {
    std::vector<double> U1(mesh.Ne * 4);
    std::vector<double> U2(mesh.Ne * 4);
    std::vector<double> R(mesh.Ne * 4);
    std::vector<double> dt_local(mesh.Ne);

    calcRes(mesh, U, R.data(), gammad, params, flux_fn, recon_fn, dt_local.data(), CFL);
    for (int i = 0; i < mesh.Ne; ++i) {
        double dt = dt_local[i];
        for (int k = 0; k < 4; ++k)
            U1[i * 4 + k] = U[i * 4 + k] + dt * R[i * 4 + k];
    }

    calcRes(mesh, U1.data(), R.data(), gammad, params, flux_fn, recon_fn, dt_local.data(), CFL);
    for (int i = 0; i < mesh.Ne; ++i) {
        double dt = dt_local[i];
        for (int k = 0; k < 4; ++k)
            U2[i * 4 + k] = 0.75 * U[i * 4 + k] + 0.25 * (U1[i * 4 + k] + dt * R[i * 4 + k]);
    }

    calcRes(mesh, U2.data(), R.data(), gammad, params, flux_fn, recon_fn, dt_local.data(), CFL);
    for (int i = 0; i < mesh.Ne; ++i) {
        double dt = dt_local[i];
        for (int k = 0; k < 4; ++k)
            U[i * 4 + k] = (1.0 / 3.0) * U[i * 4 + k] + (2.0 / 3.0) * (U2[i * 4 + k] + dt * R[i * 4 + k]);
    }
}

double residual_L1_norm(const GriMesh& mesh, const double* R) {
    double sum = 0.0;
    for (int i = 0; i < mesh.Ne * 4; ++i)
        sum += std::fabs(R[i]);
    return sum;
}

double residual_L2_norm(const GriMesh& mesh, const double* R) {
    double sum = 0.0;
    for (int i = 0; i < mesh.Ne * 4; ++i)
        sum += R[i] * R[i];
    return std::sqrt(sum);
}

double compute_dt(const GriMesh& mesh, const double* U, double gammad, double CFL) {
    double h_min_sq = mesh.Area[0];
    for (int i = 1; i < mesh.Ne; ++i)
        h_min_sq = std::min(h_min_sq, mesh.Area[i]);
    double h_min = std::sqrt(h_min_sq);
    double speed_max = 0.0;
    for (int i = 0; i < mesh.Ne; ++i) {
        double rho, u, v, p, c;
        consToPrim(&U[i * 4], gammad, rho, u, v, p, c);
        speed_max = std::max(speed_max, std::sqrt(u * u + v * v) + c);
    }
    if (speed_max < 1e-14) speed_max = 1e-14;
    return CFL * h_min / speed_max;
}

void solve_steady(const GriMesh& mesh, double* U, double gammad, const ProblemParams& params,
                  FluxFn flux_fn, ReconFn recon_fn, double CFL, int residual_stride, int max_iter) {
    std::vector<double> R(mesh.Ne * 4);
    calcRes(mesh, U, R.data(), gammad, params, flux_fn, recon_fn);
    double R0 = residual_L1_norm(mesh, R.data());
    std::cout << "Initial L1 residual: " << R0 << "\n";

    int step = 0;
    double t = 0.0;

    while (step < max_iter) {
        SSPRK3_local(mesh, U, gammad, params, flux_fn, recon_fn, CFL);
        t += compute_dt(mesh, U, gammad, CFL);
        step++;

        if (residual_stride > 0 && step % residual_stride == 0) {
            calcRes(mesh, U, R.data(), gammad, params, flux_fn, recon_fn);
            double R1 = residual_L1_norm(mesh, R.data());
            std::cout << "Step " << step << "  t=" << t << "  L1=" << R1;
            if (R0 > 1e-30)
                std::cout << "  ratio=" << (R1 / R0);
            std::cout << "\n";
            if (R1 < R0 * 1e-5) {
                std::cout << "Converged (L1 < 1e-5 * R0).\n";
                break;
            }
        }
    }
}

void solve_steady_2nd(const GriMesh& mesh, double* U, double gammad, const ProblemParams& params,
                      FluxFn flux_fn, ReconFn recon_2nd, double CFL, int residual_stride,
                      int max_iter, int max_iter_1st) {
    std::vector<double> R(mesh.Ne * 4);
    calcRes(mesh, U, R.data(), gammad, params, flux_fn, reconstruct_const);
    double R0 = residual_L1_norm(mesh, R.data());
    std::cout << "solve_steady_2nd: Initial L1 (1st order) = " << R0 << "\n";
    std::cout << "  Phase 1: 1st order for " << max_iter_1st << " steps, then 2nd order until convergence.\n";

    int step = 0;
    double t = 0.0;
    ReconFn recon_fn = reconstruct_const;

    while (step < max_iter) {
        SSPRK3_local(mesh, U, gammad, params, flux_fn, recon_fn, CFL);
        t += compute_dt(mesh, U, gammad, CFL);
        step++;

        if (step == max_iter_1st && recon_fn == reconstruct_const) {
            recon_fn = recon_2nd;
            calcRes(mesh, U, R.data(), gammad, params, flux_fn, recon_fn);
            double R1 = residual_L1_norm(mesh, R.data());
            std::cout << "  Phase 2 (2nd order): Step " << step << "  t=" << t << "  L1=" << R1;
            if (R0 > 1e-30) std::cout << "  ratio=" << (R1 / R0);
            std::cout << "\n";
        }

        if (residual_stride > 0 && step % residual_stride == 0) {
            calcRes(mesh, U, R.data(), gammad, params, flux_fn, recon_fn);
            double R1 = residual_L1_norm(mesh, R.data());
            std::cout << "Step " << step << "  t=" << t << "  L1=" << R1;
            if (R0 > 1e-30)
                std::cout << "  ratio=" << (R1 / R0);
            std::cout << (recon_fn == reconstruct_const ? "  [1st]" : "  [2nd]") << "\n";
            if (R1 < R0 * 1e-5) {
                std::cout << "Converged (L1 < 1e-5 * R0).\n";
                break;
            }
        }
    }
}

void solve_unsteady(const GriMesh& mesh, double* U, double gammad,
                    const ProblemParams& params, FluxFn flux_fn, ReconFn recon_fn,
                    double CFL, double t_end, double vtu_interval, int residual_stride,
                    const char* out_dir) {
    std::vector<double> R(mesh.Ne * 4);
    std::vector<double> vtu_times;
    char path[512];

    std::ofstream hist(std::string(out_dir) + "/residual_history.dat");
    if (hist.is_open())
        hist << "# step  t  L1  L2\n";

    double t = 0.0;
    int step = 0;
    double next_vtu_t = vtu_interval;

    /* Output initial */
    std::snprintf(path, sizeof(path), "%s/solution_t_%06.2f.vtu", out_dir, t);
    if (write_vtu(mesh, U, gammad, path))
        vtu_times.push_back(t);

    calcRes(mesh, U, R.data(), gammad, params, flux_fn, recon_fn, nullptr, CFL, t);
    double R0 = residual_L1_norm(mesh, R.data());
    std::cout << "Unsteady: t=" << t << "  L1=" << R0 << "\n";
    if (hist.is_open())
        hist << "0  " << t << "  " << R0 << "  " << residual_L2_norm(mesh, R.data()) << "\n";

    while (t < t_end) {
        double dt = SSPRK3(mesh, U, gammad, params, flux_fn, recon_fn, CFL, t);
        t += dt;
        step++;

        if (residual_stride > 0 && step % residual_stride == 0) {
            calcRes(mesh, U, R.data(), gammad, params, flux_fn, recon_fn, nullptr, CFL, t);
            double R1 = residual_L1_norm(mesh, R.data());
            double R1_L2 = residual_L2_norm(mesh, R.data());
            std::cout << "  Step " << step << "  t=" << t << "  L1=" << R1 << "  L2=" << R1_L2 << "\n";
            if (hist.is_open())
                hist << step << "  " << t << "  " << R1 << "  " << R1_L2 << "\n";
        }

        if (t >= next_vtu_t - 1e-12 || t >= t_end - 1e-12) {
            double out_t = (t >= t_end - 1e-12) ? t_end : next_vtu_t;
            calcRes(mesh, U, R.data(), gammad, params, flux_fn, recon_fn, nullptr, CFL, t);
            if (hist.is_open())
                hist << step << "  " << t << "  " << residual_L1_norm(mesh, R.data()) << "  " << residual_L2_norm(mesh, R.data()) << "\n";
            std::snprintf(path, sizeof(path), "%s/solution_t_%06.2f.vtu", out_dir, out_t);
            if (write_vtu(mesh, U, gammad, path))
                vtu_times.push_back(out_t);
            next_vtu_t += vtu_interval;
        }
    }

    /* Write PVD */
    std::ofstream pvd(std::string(out_dir) + "/solution.pvd");
    if (pvd.is_open()) {
        pvd << "<?xml version=\"1.0\"?>\n";
        pvd << "<VTKFile type=\"Collection\" version=\"1.0\">\n";
        pvd << "  <Collection>\n";
        for (size_t i = 0; i < vtu_times.size(); ++i) {
            char fname[64];
            std::snprintf(fname, sizeof(fname), "solution_t_%06.2f.vtu", vtu_times[i]);
            pvd << "    <DataSet timestep=\"" << vtu_times[i] << "\" part=\"0\" file=\"" << fname << "\"/>\n";
        }
        pvd << "  </Collection>\n";
        pvd << "</VTKFile>\n";
    }
    std::cout << "Unsteady done. Output: " << out_dir << "/ (" << vtu_times.size() << " VTU, solution.pvd)\n";
}