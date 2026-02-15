#include "physics.hpp"
#include <cmath>
#include <algorithm>

void consToPrim(const double U[4], double gamma,
                  double& rho, double& u, double& v, double& p, double& c) {
    const double rho_floor = 1e-12;
    const double p_floor = 1e-12;

    rho = std::max(U[0], rho_floor);
    double mom_x = U[1];
    double mom_y = U[2];
    double E = U[3];

    u = mom_x / rho;
    v = mom_y / rho;

    p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
    p = std::max(p, p_floor);

    c = std::sqrt(gamma * p / rho);
}

void primToCons(double rho, double u, double v, double p, double gamma, double U[4]) {
    double E = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = E;
}

double getp0(double rho0, double a0, double gamma) {
    return rho0 * a0 * a0 / gamma;
}

void isentropic_prim_from_M(double rho0, double p0, double gamma, double M, double alpha,
                            double& rho, double& u, double& v, double& p) {
    /* fac = 1 + (g-1)/2 * M^2 */
    double fac = 1.0 + 0.5 * (gamma - 1.0) * M * M;
    p = p0 / std::pow(fac, gamma / (gamma - 1.0));
    rho = rho0 / std::pow(fac, 1.0 / (gamma - 1.0));
    double a = std::sqrt(gamma * p / rho);
    u = M * a * std::cos(alpha);
    v = M * a * std::sin(alpha);
}

double total_temperature(double rho, double u, double v, double p, double gamma){
    /* Tt/T = 1 + (g-1)/2 * M^2, M^2 = (u^2+v^2)/a^2, a^2 = gamma*p/rho */
    const double eps = 1e-14;
    double a2 = gamma * p / (rho + eps);
    a2 = std::max(a2, eps);

    double M2 = (u*u + v*v) / a2;
    double fac = 1.0 + 0.5*(gamma-1.0)*M2;

    /* Tt = T * fac, T = gamma*p/rho (R=1/gamma) */
    return gamma * (p / (rho + eps)) * fac;
}

double total_pressure(double rho, double u, double v, double p, double gamma) {
    /* pt = p * (1 + (g-1)/2 * M^2)^(g/(g-1)), M^2 = (u^2+v^2)/a^2 */
    const double eps = 1e-14;
    double a2 = gamma * p / (rho + eps);
    if (a2 < eps) a2 = eps;
    double M2 = (u * u + v * v) / a2;
    double fac = 1.0 + 0.5 * (gamma - 1.0) * M2;
    return p * std::pow(fac, gamma / (gamma - 1.0));
}

void physicalFlux(const double U[4], const double n[2], double gamma, double F[4]) {
    double nx = n[0];
    double ny = n[1];

    double rho, u, v, p, c;
    consToPrim(U, gamma, rho, u, v, p, c);
    double E = U[3];

    double un = u * nx + v * ny;

    F[0] = rho * un;
    F[1] = rho * u * un + p * nx;
    F[2] = rho * v * un + p * ny;
    F[3] = (E + p) * un;
}
