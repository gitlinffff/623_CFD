/**
 * Verify boundary flux functions: WallFlux, InflowFlux, OutflowFlux.
 * Run from project root: ./build/tests/test_bc_flux
 */
#include "flux.hpp"
#include "physics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

static const double gamma = 1.4;
static const double tol = 1e-10;
static const double tol_machine = 1e-14;  /* machine precision for exactness tests */

static void U_from_prim(double rho, double u, double v, double p, double U[4]) {
    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
}

/* ---------------------------------------------------------------------------
 * WallFlux tests
 * ---------------------------------------------------------------------------
 * Wall: F·n = [0, p*nx, p*ny, 0] (zero mass/energy flux, momentum = pressure)
 */
static bool test_wall_zero_mass_flux() {
    double U[4];
    U_from_prim(1.0, 0.5, 0.3, 1.0, U);
    double n[2] = {1.0, 0.0};
    double Fhat[4];
    double smag;
    WallFlux(U, n, gamma, Fhat, smag);
    bool ok = std::fabs(Fhat[0]) < tol;
    std::printf("  WallFlux: Fhat[0] (mass) = %.2e  %s\n", Fhat[0], ok ? "PASS" : "FAIL");
    return ok;
}

static bool test_wall_zero_energy_flux() {
    double U[4];
    U_from_prim(1.0, 0.5, 0.3, 1.0, U);
    double n[2] = {1.0, 0.0};
    double Fhat[4];
    double smag;
    WallFlux(U, n, gamma, Fhat, smag);
    bool ok = std::fabs(Fhat[3]) < tol;
    std::printf("  WallFlux: Fhat[3] (energy) = %.2e  %s\n", Fhat[3], ok ? "PASS" : "FAIL");
    return ok;
}

static bool test_wall_momentum_equals_pressure() {
    /* Use flow parallel to wall (un=0) so wall pressure = interior pressure */
    double U[4];
    double rho = 1.0, p = 0.7;
    double n[2] = {1.0, 0.0};  /* wall normal */
    double u = 0.0, v = 0.3;   /* tangent flow, un = u*nx+v*ny = 0 */
    U_from_prim(rho, u, v, p, U);
    double Fhat[4];
    double smag;
    WallFlux(U, n, gamma, Fhat, smag);
    bool ok = (std::fabs(Fhat[1] - p * n[0]) < tol && std::fabs(Fhat[2] - p * n[1]) < tol);
    std::printf("  WallFlux: Fhat[1,2] = p*n (parallel flow)  err=%.2e  %s\n",
                std::max(std::fabs(Fhat[1] - p * n[0]), std::fabs(Fhat[2] - p * n[1])),
                ok ? "PASS" : "FAIL");
    return ok;
}

/* Rotated Wall Exactness: n=(cos θ, sin θ), Fhat must match [0, p_L*nx, p_L*ny, 0] to 1e-14.
 * Use u=v=0 so u·n=0 for all n, ensuring wall pressure = interior pressure. */
static bool test_wall_rotated_exactness() {
    double U[4];
    double rho = 1.0, u = 0.0, v = 0.0, p = 0.85;
    U_from_prim(rho, u, v, p, U);
    double rhoL, uL, vL, pL, cL;
    consToPrim(U, gamma, rhoL, uL, vL, pL, cL);
    double max_err = 0.0;
    for (int k = 0; k < 8; ++k) {
        double theta = 0.3 + 0.7 * k;  /* arbitrary angles */
        double n[2] = {std::cos(theta), std::sin(theta)};
        double Fhat[4];
        double smag;
        WallFlux(U, n, gamma, Fhat, smag);
        max_err = std::max(max_err, std::fabs(Fhat[0]));
        max_err = std::max(max_err, std::fabs(Fhat[1] - pL * n[0]));
        max_err = std::max(max_err, std::fabs(Fhat[2] - pL * n[1]));
        max_err = std::max(max_err, std::fabs(Fhat[3]));
    }
    bool ok = max_err < tol_machine;
    std::printf("  WallFlux: rotated exactness (8 angles)  max_err=%.2e  %s\n",
                max_err, ok ? "PASS" : "FAIL");
    return ok;
}

/* ---------------------------------------------------------------------------
 * InflowFlux tests
 * ---------------------------------------------------------------------------
 * When interior state UL equals the inflow state Ub, flux should match physical flux.
 */
static bool test_inflow_consistency() {
    const double rho0 = 1.0, a0 = 1.0, R = 1.0 / gamma;
    const double alpha = 50.0 * 3.14159265 / 180.0;
    double rho, u, v, p;
    isentropic_prim_from_M(rho0, rho0 * a0 * a0 / gamma, gamma, 0.1, alpha, rho, u, v, p);
    double nin[2] = {std::cos(alpha), std::sin(alpha)};
    double UL[4];
    U_from_prim(rho, u, v, p, UL);
    double n[2] = {-nin[0], -nin[1]};  /* outward normal opposite to inflow */
    double Fhat[4], Fphys[4], smag;
    InflowFlux(UL, n, nin, rho0, a0, gamma, R, fluxROE, Fhat, smag);
    physicalFlux(UL, n, gamma, Fphys);
    double err = 0.0;
    for (int k = 0; k < 4; ++k) err = std::max(err, std::fabs(Fhat[k] - Fphys[k]));
    bool ok = err < 1e-8;  /* looser: inflow may differ slightly when UL≈Ub */
    std::printf("  InflowFlux: |Fhat - F(UL)| when UL=inflow  err=%.2e  %s\n", err, ok ? "PASS" : "FAIL");
    return ok;
}

static bool test_inflow_flux_direction() {
    const double rho0 = 1.0, a0 = 1.0, R = 1.0 / gamma;
    /* Use state very close to inflow (M~0.1) to avoid quadratic discriminant issues */
    double rho, u, v, p;
    isentropic_prim_from_M(rho0, rho0 * a0 * a0 / gamma, gamma, 0.08,
                           50.0 * 3.14159265 / 180.0, rho, u, v, p);
    double UL[4];
    U_from_prim(rho, u, v, p, UL);
    double nin[2] = {1.0, 0.0};
    double n[2] = {-1.0, 0.0};  /* outward normal */
    double Fhat[4], smag;
    try {
        InflowFlux(UL, n, nin, rho0, a0, gamma, R, fluxROE, Fhat, smag);
    } catch (const std::exception& e) {
        std::printf("  InflowFlux: threw %s  FAIL\n", e.what());
        return false;
    }
    /* For inflow, Fhat = F·n with n outward: mass flux < 0 (flow into domain) */
    bool ok = Fhat[0] < 0.0;
    std::printf("  InflowFlux: Fhat[0] (mass) < 0 for inflow  Fhat[0]=%.4e  %s\n", Fhat[0], ok ? "PASS" : "FAIL");
    return ok;
}

/* Inflow Reverse Engineering: total_pressure(Ub)=Pt, total_temperature(Ub)=Tt, v_b || nin */
static bool test_inflow_reverse_engineering() {
    const double rho0 = 1.0, a0 = 1.0, R = 1.0 / gamma;
    const double input_Pt = rho0 * a0 * a0 / gamma;
    const double input_Tt = (a0 * a0) / (gamma * R);
    const double alpha = 50.0 * 3.14159265 / 180.0;
    double nin[2] = {std::cos(alpha), std::sin(alpha)};
    double rho, u, v, p;
    isentropic_prim_from_M(rho0, input_Pt, gamma, 0.1, alpha, rho, u, v, p);
    double UL[4];
    U_from_prim(rho, u, v, p, UL);
    double n[2] = {-nin[0], -nin[1]};
    double Ub[4];
    try {
        InflowFlux_compute_Ub(UL, n, nin, rho0, a0, gamma, R, Ub);
    } catch (const std::exception& e) {
        std::printf("  InflowFlux: reverse eng threw %s  FAIL\n", e.what());
        return false;
    }
    double rhob, ub, vb, pb, cb;
    consToPrim(Ub, gamma, rhob, ub, vb, pb, cb);
    double pt_Ub = total_pressure(rhob, ub, vb, pb, gamma);
    double Tt_Ub = total_temperature(rhob, ub, vb, pb, gamma);
    double err_pt = std::fabs(pt_Ub - input_Pt);
    double err_Tt = std::fabs(Tt_Ub - input_Tt);
    double vmag = std::sqrt(ub*ub + vb*vb);
    double nin_mag = std::sqrt(nin[0]*nin[0] + nin[1]*nin[1]);
    double err_dir = (vmag > 1e-14 && nin_mag > 1e-14) ?
        std::fabs(ub*nin[1] - vb*nin[0]) : 0.0;  /* cross product = 0 if parallel */
    bool ok = (err_pt < tol_machine && err_Tt < 5e-14 && err_dir < tol_machine);
    std::printf("  InflowFlux: reverse eng  |pt-pt_in|=%.2e |Tt-Tt_in|=%.2e |v×nin|=%.2e  %s\n",
                err_pt, err_Tt, err_dir, ok ? "PASS" : "FAIL");
    return ok;
}

/* ---------------------------------------------------------------------------
 * OutflowFlux tests
 * ---------------------------------------------------------------------------
 * When pL = pout, boundary state Ub ≈ UL, flux should match physical flux.
 */
static bool test_outflow_consistency() {
    const double pout = 0.7;
    double UL[4];
    U_from_prim(0.9, 0.3, 0.1, pout, UL);  /* pL = pout */
    double n[2] = {1.0, 0.0};
    double Fhat[4], Fphys[4], smag;
    OutflowFlux(UL, n, pout, gamma, fluxROE, Fhat, smag);
    physicalFlux(UL, n, gamma, Fphys);
    double err = 0.0;
    for (int k = 0; k < 4; ++k) err = std::max(err, std::fabs(Fhat[k] - Fphys[k]));
    bool ok = err < 1e-8;
    std::printf("  OutflowFlux: |Fhat - F(UL)| when pL=pout  err=%.2e  %s\n", err, ok ? "PASS" : "FAIL");
    return ok;
}

static bool test_outflow_flux_direction() {
    const double pout = 0.7;
    double UL[4];
    U_from_prim(0.9, 0.3, 0.1, 0.8, UL);  /* pL > pout, flow out */
    double n[2] = {1.0, 0.0};
    double Fhat[4], smag;
    OutflowFlux(UL, n, pout, gamma, fluxROE, Fhat, smag);
    /* For outflow, mass flux Fhat[0] should be positive (flow out, n outward) */
    bool ok = Fhat[0] > 0.0;
    std::printf("  OutflowFlux: Fhat[0] (mass) > 0 for outflow  Fhat[0]=%.4e  %s\n", Fhat[0], ok ? "PASS" : "FAIL");
    return ok;
}

/* Outflow Riemann Invariant: p_b=pout, S_b=S_L, J_b^+=J_L^+ */
static bool test_outflow_riemann_invariant() {
    const double pout = 0.7;
    double UL[4];
    U_from_prim(0.9, 0.3, 0.1, 0.8, UL);
    double n[2] = {1.0, 0.0};
    double Ub[4];
    OutflowFlux_compute_Ub(UL, n, pout, gamma, Ub);
    double rhoL, uL, vL, pL, cL;
    consToPrim(UL, gamma, rhoL, uL, vL, pL, cL);
    double rhob, ub, vb, pb, cb;
    consToPrim(Ub, gamma, rhob, ub, vb, pb, cb);
    double SL = pL / std::pow(rhoL, gamma);
    double Sb = pb / std::pow(rhob, gamma);
    double unL = uL*n[0] + vL*n[1];
    double unb = ub*n[0] + vb*n[1];
    double J_L = unL + 2*cL/(gamma - 1);
    double J_b = unb + 2*cb/(gamma - 1);
    double err_pb = std::fabs(pb - pout);
    double err_S = std::fabs(Sb - SL);
    double err_J = std::fabs(J_b - J_L);
    bool ok = (err_pb < tol_machine && err_S < tol_machine && err_J < tol_machine);
    std::printf("  OutflowFlux: Riemann inv  |pb-pout|=%.2e |Sb-SL|=%.2e |Jb-JL|=%.2e  %s\n",
                err_pb, err_S, err_J, ok ? "PASS" : "FAIL");
    return ok;
}

int main() {
    std::printf("=== WallFlux ===\n");
    int nw = 0;
    nw += test_wall_zero_mass_flux();
    nw += test_wall_zero_energy_flux();
    nw += test_wall_momentum_equals_pressure();
    nw += test_wall_rotated_exactness();
    std::printf("  Result: %d/4\n\n", nw);

    std::printf("=== InflowFlux ===\n");
    int ni = 0;
    ni += test_inflow_consistency();
    ni += test_inflow_flux_direction();
    ni += test_inflow_reverse_engineering();
    std::printf("  Result: %d/3\n\n", ni);

    std::printf("=== OutflowFlux ===\n");
    int no = 0;
    no += test_outflow_consistency();
    no += test_outflow_flux_direction();
    no += test_outflow_riemann_invariant();
    std::printf("  Result: %d/3\n\n", no);

    int total = nw + ni + no;
    std::printf("=== Summary ===\n");
    std::printf("  WallFlux: %d/4\n", nw);
    std::printf("  InflowFlux: %d/3\n", ni);
    std::printf("  OutflowFlux: %d/3\n", no);
    std::printf("  Total: %d/10  %s\n", total, total == 10 ? "ALL PASS" : "FAIL");

    return (total == 10) ? 0 : 1;
}
