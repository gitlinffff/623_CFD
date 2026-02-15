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
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

void calcRes(const GriMesh& mesh, const double* U, double* R, double gamma,
             const ProblemParams& params, FluxFn flux_fn, ReconFn recon_fn,
             double* dt_per_cell, double CFL) {
    std::memset(R, 0, mesh.Ne * 4 * sizeof(double));
    std::vector<double> sum_s(mesh.Ne, 0.0);

    std::vector<double> UL_int(mesh.num_interior_faces * 4);
    std::vector<double> UR_int(mesh.num_interior_faces * 4);
    std::vector<double> UL_bnd(mesh.num_boundary_faces * 4);
    std::vector<double> UR_bnd(mesh.num_boundary_faces * 4);

    recon_fn(mesh, U, UL_int.data(), UR_int.data(), UL_bnd.data(), UR_bnd.data());

    double Fhat[4];
    double smag;
    const double nin[2] = {std::cos(params.alpha), std::sin(params.alpha)};
    const double R_gas = 1.0 / gamma;

    for (int i = 0; i < mesh.num_interior_faces; ++i) {
        flux_fn(&UL_int[i * 4], &UR_int[i * 4], &mesh.In[2 * i], gamma, Fhat, smag);
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
            WallFlux(UL, n, gamma, Fhat, smag);
        else if (bt == BC_INFLOW) {
            try {
                InflowFlux(UL, n, nin, params.rho0, params.a0, gamma, R_gas, flux_fn, Fhat, smag);
            } catch (const std::runtime_error&) {
                flux_fn(UL, UL, n, gamma, Fhat, smag);
            }
        } else if (bt == BC_OUTFLOW)
            OutflowFlux(UL, n, params.pout, gamma, flux_fn, Fhat, smag);
        else
            flux_fn(UL, UL, n, gamma, Fhat, smag);

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

void SSPRK3(const GriMesh& mesh, double* U, double gamma,
            const ProblemParams& params, FluxFn flux_fn, ReconFn recon_fn, double CFL) {
    std::vector<double> U1(mesh.Ne * 4);
    std::vector<double> U2(mesh.Ne * 4);
    std::vector<double> R(mesh.Ne * 4);
    std::vector<double> dt_local(mesh.Ne);

    calcRes(mesh, U, R.data(), gamma, params, flux_fn, recon_fn, dt_local.data(), CFL);
    for (int i = 0; i < mesh.Ne; ++i) {
        double dt = dt_local[i];
        for (int k = 0; k < 4; ++k)
            U1[i * 4 + k] = U[i * 4 + k] + dt * R[i * 4 + k];
    }

    calcRes(mesh, U1.data(), R.data(), gamma, params, flux_fn, recon_fn, dt_local.data(), CFL);
    for (int i = 0; i < mesh.Ne; ++i) {
        double dt = dt_local[i];
        for (int k = 0; k < 4; ++k)
            U2[i * 4 + k] = 0.75 * U[i * 4 + k] + 0.25 * (U1[i * 4 + k] + dt * R[i * 4 + k]);
    }

    calcRes(mesh, U2.data(), R.data(), gamma, params, flux_fn, recon_fn, dt_local.data(), CFL);
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

double compute_dt(const GriMesh& mesh, const double* U, double gamma, double CFL) {
    double h_min_sq = mesh.Area[0];
    for (int i = 1; i < mesh.Ne; ++i)
        h_min_sq = std::min(h_min_sq, mesh.Area[i]);
    double h_min = std::sqrt(h_min_sq);
    double speed_max = 0.0;
    for (int i = 0; i < mesh.Ne; ++i) {
        double rho, u, v, p, c;
        consToPrim(&U[i * 4], gamma, rho, u, v, p, c);
        speed_max = std::max(speed_max, std::sqrt(u * u + v * v) + c);
    }
    if (speed_max < 1e-14) speed_max = 1e-14;
    return CFL * h_min / speed_max;
}

void solve_steady(const GriMesh& mesh, double* U, double gamma, const ProblemParams& params,
                  FluxFn flux_fn, ReconFn recon_fn, double CFL, int residual_stride, int max_iter) {
    std::vector<double> R(mesh.Ne * 4);
    calcRes(mesh, U, R.data(), gamma, params, flux_fn, recon_fn);
    double R0 = residual_L1_norm(mesh, R.data());
    std::cout << "Initial L1 residual: " << R0 << "\n";

    int step = 0;
    double t = 0.0;

    while (step < max_iter) {
        SSPRK3(mesh, U, gamma, params, flux_fn, recon_fn, CFL);
        t += compute_dt(mesh, U, gamma, CFL);
        step++;

        if (residual_stride > 0 && step % residual_stride == 0) {
            calcRes(mesh, U, R.data(), gamma, params, flux_fn, recon_fn);
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
