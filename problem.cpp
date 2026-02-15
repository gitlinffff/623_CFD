#include "problem.hpp"
#include "physics.hpp"
#include <cstring>

namespace {
const double PI = 3.14159265358979323846;
}

ProblemParams::ProblemParams() {
    rho0 = 1.0;
    a0 = 1.0;
    gamma = 1.4;
    alpha = 50.0 * PI / 180.0;  /* 50 deg */
    double p0 = getp0(rho0, a0, gamma);
    pout = 0.7 * p0;
}

double getp0(const ProblemParams& p) {
    return getp0(p.rho0, p.a0, p.gamma);
}

void initialize_uniform(double* U, int Ne, double M, const ProblemParams& params) {
    double rho, u, v, p;
    double p0 = getp0(params.rho0, params.a0, params.gamma);
    isentropic_prim_from_M(params.rho0, p0, params.gamma, M, params.alpha,
                          rho, u, v, p);

    double Ucell[4];
    primToCons(rho, u, v, p, params.gamma, Ucell);

    for (int i = 0; i < Ne; ++i) {
        std::memcpy(U + i * 4, Ucell, 4 * sizeof(double));
    }
}
