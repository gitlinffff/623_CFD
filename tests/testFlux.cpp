#include "flux.hpp"
#include "physics.hpp"
#include <iostream>
#include <cmath>

static const double gammad = 1.4;
static const double tol = 1e-10;

static void U_from_prim(double rho, double u, double v, double p, double U[4]) {
    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = p / (gammad - 1.0) + 0.5 * rho * (u * u + v * v);
}

static double vec4_diff(const double a[4], const double b[4]) {
    double d = 0;
    for (int k = 0; k < 4; ++k) {
        double t = std::abs(a[k] - b[k]);
        if (t > d) d = t;
    }
    return d;
}

/* Consistency: when UL = UR, Fhat(UL, UR, n) = F(U, n) = physical flux */
static bool test_consistency(const char* scheme, void (*flux)(const double*, const double*, const double*, double, double*, double&)) {
    double U[4];
    U_from_prim(1.0, 0.5, 0.3, 1.0, U);
    double n[2] = { 1.0, 0.0 };

    double Fhat[4], Fphys[4], smag;
    flux(U, U, n, gammad, Fhat, smag);
    physicalFlux(U, n, gammad, Fphys);

    double err = vec4_diff(Fhat, Fphys);
    bool ok = (err < tol);
    std::cout << "  " << scheme << ": max|Fhat - F| = " << err << (ok ? "  PASS" : "  FAIL") << "\n";
    return ok;
}

/* Exchange symmetry: Fhat(UL, UR, n) = -Fhat(UR, UL, -n) */
static bool test_exchange_symmetry(const char* scheme, void (*flux)(const double*, const double*, const double*, double, double*, double&)) {
    double UL[4], UR[4];
    U_from_prim(1.0, 0.0, 0.0, 1.0, UL);
    U_from_prim(0.125, 0.0, 0.0, 0.1, UR);
    double n[2] = { 1.0, 0.0 };
    double n_neg[2] = { -1.0, 0.0 };

    double Fhat_LR[4], Fhat_RL[4], smag;
    flux(UL, UR, n, gammad, Fhat_LR, smag);
    flux(UR, UL, n_neg, gammad, Fhat_RL, smag);
    for (int k = 0; k < 4; ++k) Fhat_RL[k] = -Fhat_RL[k];

    double err = vec4_diff(Fhat_LR, Fhat_RL);
    bool ok = (err < tol);
    std::cout << "  " << scheme << ": max|Fhat(L,R,n) + Fhat(R,L,-n)| = " << err << (ok ? "  PASS" : "  FAIL") << "\n";
    return ok;
}

/* Rotation invariance: Fhat(UL_rot, UR_rot, n_rot) = R(Fhat(UL, UR, n)) */
static bool test_rotation_invariance(const char* scheme, void (*flux)(const double*, const double*, const double*, double, double*, double&)) {
    double UL[4], UR[4];
    U_from_prim(1.0, 0.5, 0.2, 1.0, UL);
    U_from_prim(0.5, 0.3, 0.1, 0.5, UR);
    double n[2] = { 1.0, 0.0 };

    const double theta = 0.5;  /* rad */
    double c = std::cos(theta), s = std::sin(theta);

    double UL_rot[4], UR_rot[4], n_rot[2];
    UL_rot[0] = UL[0];
    UL_rot[1] = UL[1] * c - UL[2] * s;
    UL_rot[2] = UL[1] * s + UL[2] * c;
    UL_rot[3] = UL[3];
    UR_rot[0] = UR[0];
    UR_rot[1] = UR[1] * c - UR[2] * s;
    UR_rot[2] = UR[1] * s + UR[2] * c;
    UR_rot[3] = UR[3];
    n_rot[0] = n[0] * c - n[1] * s;
    n_rot[1] = n[0] * s + n[1] * c;

    double Fhat_orig[4], Fhat_rot[4], smag;
    flux(UL, UR, n, gammad, Fhat_orig, smag);
    flux(UL_rot, UR_rot, n_rot, gammad, Fhat_rot, smag);

    double Fhat_expected[4];
    Fhat_expected[0] = Fhat_orig[0];
    Fhat_expected[1] = Fhat_orig[1] * c - Fhat_orig[2] * s;
    Fhat_expected[2] = Fhat_orig[1] * s + Fhat_orig[2] * c;
    Fhat_expected[3] = Fhat_orig[3];

    double err = vec4_diff(Fhat_rot, Fhat_expected);
    bool ok = (err < tol);
    std::cout << "  " << scheme << ": max|Fhat_rot - R(Fhat)| = " << err << (ok ? "  PASS" : "  FAIL") << "\n";
    return ok;
}

int main() {
    typedef void (*FluxFunc)(const double*, const double*, const double*, double, double*, double&);
    struct { const char* name; FluxFunc f; } schemes[] = {
        {"Rusanov", fluxRusanov},
        {"HLLC", fluxHLLC},
        {"ROE", fluxROE}
    };
    int nschemes = 3;

    std::cout << "=== 1. Consistency: Fhat(U,U,n) = F(U,n) when UL = UR ===\n";
    int ok1 = 0;
    for (int i = 0; i < nschemes; ++i)
        if (test_consistency(schemes[i].name, schemes[i].f)) ok1++;
    std::cout << "  Result: " << ok1 << "/" << nschemes << "\n\n";

    std::cout << "=== 2. Exchange symmetry: Fhat(UL,UR,n) = -Fhat(UR,UL,-n) ===\n";
    int ok2 = 0;
    for (int i = 0; i < nschemes; ++i)
        if (test_exchange_symmetry(schemes[i].name, schemes[i].f)) ok2++;
    std::cout << "  Result: " << ok2 << "/" << nschemes << "\n\n";

    std::cout << "=== 3. Rotation invariance: Fhat(UL',UR',n') = R(Fhat(UL,UR,n)) ===\n";
    int ok3 = 0;
    for (int i = 0; i < nschemes; ++i)
        if (test_rotation_invariance(schemes[i].name, schemes[i].f)) ok3++;
    std::cout << "  Result: " << ok3 << "/" << nschemes << "\n\n";

    int total = ok1 + ok2 + ok3;
    std::cout << "=== Summary ===\n";
    std::cout << "  Consistency: " << ok1 << "/" << nschemes << "\n";
    std::cout << "  Exchange symmetry: " << ok2 << "/" << nschemes << "\n";
    std::cout << "  Rotation invariance: " << ok3 << "/" << nschemes << "\n";
    std::cout << "  Total: " << total << "/" << (nschemes * 3) << (total == nschemes * 3 ? "  ALL PASS" : "  FAIL") << "\n";

    return (total == nschemes * 3) ? 0 : 1;
}
