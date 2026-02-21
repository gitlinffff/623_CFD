#include "reconstruct.hpp"
#include "physics.hpp"
#include <algorithm>
#include <cstring>
#include <cmath>

namespace {

/* Cell centroid for triangle i. */
void cell_centroid(const GriMesh& mesh, int i, double& cx, double& cy) {
    const int* tri = &mesh.E[i * 3];
    cx = (mesh.V[tri[0] * 2 + 0] + mesh.V[tri[1] * 2 + 0] + mesh.V[tri[2] * 2 + 0]) / 3.0;
    cy = (mesh.V[tri[0] * 2 + 1] + mesh.V[tri[1] * 2 + 1] + mesh.V[tri[2] * 2 + 1]) / 3.0;
}

/* Face centroid for interior face: edge between elemL's face faceL.
 * Convention (matches readgri): face j = edge (tri[j], tri[(j+1)%3]). */
void interior_face_centroid(const GriMesh& mesh, int elemL, int faceL, double& fx, double& fy) {
    const int* tri = &mesh.E[elemL * 3];
    int v0 = tri[faceL], v1 = tri[(faceL + 1) % 3];
    fx = 0.5 * (mesh.V[v0 * 2 + 0] + mesh.V[v1 * 2 + 0]);
    fy = 0.5 * (mesh.V[v0 * 2 + 1] + mesh.V[v1 * 2 + 1]);
}

/* Face centroid for boundary face.
 * Convention (matches readgri): face j = edge (tri[j], tri[(j+1)%3]). */
void boundary_face_centroid(const GriMesh& mesh, int elem, int face, double& fx, double& fy) {
    const int* tri = &mesh.E[elem * 3];
    int v0 = tri[face], v1 = tri[(face + 1) % 3];
    fx = 0.5 * (mesh.V[v0 * 2 + 0] + mesh.V[v1 * 2 + 0]);
    fy = 0.5 * (mesh.V[v0 * 2 + 1] + mesh.V[v1 * 2 + 1]);
}

bool is_valid_state(const double U[4], double gamma) {
    if (U[0] <= 0.0) return false;
    double rho, u, v, p, c;
    consToPrim(U, gamma, rho, u, v, p, c);
    return (rho > 0.0 && p > 0.0 && std::isfinite(rho) && std::isfinite(p));
}

}  // namespace

void reconstruct_const(const GriMesh& mesh, const double* U,
                      double* UL_int, double* UR_int,
                      double* UL_bnd, double* UR_bnd, double /*gamma*/) {
    for (int i = 0; i < mesh.num_interior_faces; ++i) {
        int elemL = mesh.I2E[4 * i + 0];
        int elemR = mesh.I2E[4 * i + 2];
        std::memcpy(UL_int + i * 4, U + elemL * 4, 4 * sizeof(double));
        std::memcpy(UR_int + i * 4, U + elemR * 4, 4 * sizeof(double));
    }
    for (int i = 0; i < mesh.num_boundary_faces; ++i) {
        int elem = mesh.B2E[3 * i + 0];
        std::memcpy(UL_bnd + i * 4, U + elem * 4, 4 * sizeof(double));
        std::memcpy(UR_bnd + i * 4, U + elem * 4, 4 * sizeof(double));
    }
}

/* Max gradient magnitude (safeguard against ill-conditioned stencil / QR singularity).
 * If |grad U_k| > U_scale * grad_max_scale, clamp to avoid singular values. */
static const double grad_max_scale = 1e6;
/* Min gradient: below this, treat as zero (avoids spurious gradients in uniform flow). */
static const double grad_min_scale = 1e-14;

void compute_gradients(const GriMesh& mesh, const double* U, double* gradU) {
    /* Green-Gauss: gradU_i = (1/A_i) * sum_f U_face * n_f * L_f
     * U_face = (U_L + U_R)/2 for interior; U_cell for boundary.
     * n_f = outward normal from cell i, L_f = face length.
     */
    std::memset(gradU, 0, mesh.Ne * 8 * sizeof(double));

    for (int f = 0; f < mesh.num_interior_faces; ++f) {
        int elemL = mesh.I2E[4 * f + 0];
        int elemR = mesh.I2E[4 * f + 2];
        double nx = mesh.In[2 * f + 0];
        double ny = mesh.In[2 * f + 1];
        double L = mesh.In_len[f];
        double Uface[4];
        for (int k = 0; k < 4; ++k)
            Uface[k] = 0.5 * (U[elemL * 4 + k] + U[elemR * 4 + k]);

        /* n = In points L->R; outward from L is +n, outward from R is -n */
        for (int k = 0; k < 4; ++k) {
            double uv = Uface[k] * L;
            gradU[elemL * 8 + k * 2 + 0] += uv * nx / mesh.Area[elemL];
            gradU[elemL * 8 + k * 2 + 1] += uv * ny / mesh.Area[elemL];
            gradU[elemR * 8 + k * 2 + 0] -= uv * nx / mesh.Area[elemR];
            gradU[elemR * 8 + k * 2 + 1] -= uv * ny / mesh.Area[elemR];
        }
    }

    for (int f = 0; f < mesh.num_boundary_faces; ++f) {
        int elem = mesh.B2E[3 * f + 0];
        double nx = mesh.Bn[2 * f + 0];
        double ny = mesh.Bn[2 * f + 1];
        double L = mesh.Bn_len[f];
        for (int k = 0; k < 4; ++k) {
            double uv = U[elem * 4 + k] * L;
            gradU[elem * 8 + k * 2 + 0] += uv * nx / mesh.Area[elem];
            gradU[elem * 8 + k * 2 + 1] += uv * ny / mesh.Area[elem];
        }
    }

    /* Safeguard: clamp abnormally large gradients (ill-conditioned stencil). */
    for (int i = 0; i < mesh.Ne; ++i) {
        double h = std::sqrt(mesh.Area[i]);
        if (h < 1e-20) h = 1e-20;
        double U_scale = std::max(U[i * 4 + 0], 1e-12) / h;
        double gmax = grad_max_scale * U_scale;
        double gmin = grad_min_scale * U_scale;
        for (int k = 0; k < 4; ++k) {
            double gx = gradU[i * 8 + k * 2 + 0];
            double gy = gradU[i * 8 + k * 2 + 1];
            if (!std::isfinite(gx) || !std::isfinite(gy)) {
                gradU[i * 8 + k * 2 + 0] = 0;
                gradU[i * 8 + k * 2 + 1] = 0;
            } else {
                double mag = std::sqrt(gx * gx + gy * gy);
                if (mag < gmin) {
                    gradU[i * 8 + k * 2 + 0] = 0;
                    gradU[i * 8 + k * 2 + 1] = 0;
                } else if (mag > gmax) {
                    double fac = gmax / mag;
                    gradU[i * 8 + k * 2 + 0] = gx * fac;
                    gradU[i * 8 + k * 2 + 1] = gy * fac;
                }
            }
        }
    }
}

void reconstruct_nolimiter(const GriMesh& mesh, const double* U,
                                double* UL_int, double* UR_int,
                                double* UL_bnd, double* UR_bnd, double gamma) {
    std::vector<double> gradU(mesh.Ne * 8);
    compute_gradients(mesh, U, gradU.data());

    std::vector<double> cx(mesh.Ne), cy(mesh.Ne);
    for (int i = 0; i < mesh.Ne; ++i)
        cell_centroid(mesh, i, cx[i], cy[i]);

    auto extrapolate = [&](int elem, double fx, double fy, double* Uface) {
        double dx = fx - cx[elem], dy = fy - cy[elem];
        for (int k = 0; k < 4; ++k) {
            Uface[k] = U[elem * 4 + k] + gradU[elem * 8 + k * 2 + 0] * dx
                                           + gradU[elem * 8 + k * 2 + 1] * dy;
        }
        if (!is_valid_state(Uface, gamma))
            std::memcpy(Uface, &U[elem * 4], 4 * sizeof(double));
    };

    for (int i = 0; i < mesh.num_interior_faces; ++i) {
        int elemL = mesh.I2E[4 * i + 0];
        int elemR = mesh.I2E[4 * i + 2];
        double fx_L, fy_L;
        interior_face_centroid(mesh, elemL, mesh.I2E[4 * i + 1], fx_L, fy_L);
        extrapolate(elemL, fx_L, fy_L, &UL_int[i * 4]);
        /* Use elemR's own side face centroid to get the correct local dx for extrapolation.
         * For regular interior faces, this equals (fx_L, fy_L) (same edge).
         * For periodic faces, elemL and elemR are on opposite sides of the domain,
         * so fx_L != fx_R; using fx_L for elemR produces dx ~ domain_width, causing blow-up. */
        double fx_R, fy_R;
        interior_face_centroid(mesh, elemR, mesh.I2E[4 * i + 3], fx_R, fy_R);
        extrapolate(elemR, fx_R, fy_R, &UR_int[i * 4]);
    }

    for (int i = 0; i < mesh.num_boundary_faces; ++i) {
        int elem = mesh.B2E[3 * i + 0];
        int face = mesh.B2E[3 * i + 1];
        double fx, fy;
        boundary_face_centroid(mesh, elem, face, fx, fy);
        extrapolate(elem, fx, fy, &UL_bnd[i * 4]);
        std::memcpy(&UR_bnd[i * 4], &U[elem * 4], 4 * sizeof(double));
    }
}

void reconstruct_BJ(const GriMesh& mesh, const double* U,
                    double* UL_int, double* UR_int,
                    double* UL_bnd, double* UR_bnd, double gamma) {
    std::vector<double> gradU(mesh.Ne * 8);
    compute_gradients(mesh, U, gradU.data());

    /* U_max, U_min over cell and its face-neighbors (for Barth-Jespersen). */
    std::vector<double> U_max(mesh.Ne * 4), U_min(mesh.Ne * 4);
    for (int i = 0; i < mesh.Ne; ++i) {
        for (int k = 0; k < 4; ++k) {
            double u = U[i * 4 + k];
            U_max[i * 4 + k] = u;
            U_min[i * 4 + k] = u;
        }
    }
    for (int f = 0; f < mesh.num_interior_faces; ++f) {
        int elemL = mesh.I2E[4 * f + 0];
        int elemR = mesh.I2E[4 * f + 2];
        for (int k = 0; k < 4; ++k) {
            double uL = U[elemL * 4 + k], uR = U[elemR * 4 + k];
            U_max[elemL * 4 + k] = std::max(U_max[elemL * 4 + k], uR);
            U_min[elemL * 4 + k] = std::min(U_min[elemL * 4 + k], uR);
            U_max[elemR * 4 + k] = std::max(U_max[elemR * 4 + k], uL);
            U_min[elemR * 4 + k] = std::min(U_min[elemR * 4 + k], uL);
        }
    }

    /* Barth-Jespersen: scalar phi per cell so that limited extrapolation stays in [U_min, U_max]. */
    std::vector<double> phi_cell(mesh.Ne, 1.0);
    std::vector<double> cx(mesh.Ne), cy(mesh.Ne);
    for (int i = 0; i < mesh.Ne; ++i)
        cell_centroid(mesh, i, cx[i], cy[i]);

    for (int i = 0; i < mesh.Ne; ++i) {
        double phi = 1.0;
        for (int j = 0; j < 3; ++j) {
            double fx, fy;
            interior_face_centroid(mesh, i, j, fx, fy);
            double dx = fx - cx[i], dy = fy - cy[i];
            for (int k = 0; k < 4; ++k) {
                double u_cell = U[i * 4 + k];
                double delta = gradU[i * 8 + k * 2 + 0] * dx + gradU[i * 8 + k * 2 + 1] * dy;
                double u_face = u_cell + delta;
                double phi_k = 1.0;
                if (delta > 1e-14)
                    phi_k = (U_max[i * 4 + k] - u_cell) / delta;
                else if (delta < -1e-14)
                    phi_k = (U_min[i * 4 + k] - u_cell) / delta;
                phi_k = std::max(0.0, std::min(1.0, phi_k));
                phi = std::min(phi, phi_k);
            }
        }
        phi_cell[i] = phi;
    }

    /* Limited gradient: gradU_limited = phi * gradU */
    std::vector<double> gradU_lim(mesh.Ne * 8);
    for (int i = 0; i < mesh.Ne; ++i) {
        double p = phi_cell[i];
        for (int k = 0; k < 4; ++k) {
            gradU_lim[i * 8 + k * 2 + 0] = p * gradU[i * 8 + k * 2 + 0];
            gradU_lim[i * 8 + k * 2 + 1] = p * gradU[i * 8 + k * 2 + 1];
        }
    }

    /* Same extrapolation as nolimiter, using limited gradient; fallback to cell average if invalid. */
    auto extrapolate = [&](int elem, double fx, double fy, double* Uface) {
        double dx = fx - cx[elem], dy = fy - cy[elem];
        const double* g = &gradU_lim[elem * 8];
        for (int k = 0; k < 4; ++k) {
            Uface[k] = U[elem * 4 + k] + g[k * 2 + 0] * dx + g[k * 2 + 1] * dy;
        }
        if (!is_valid_state(Uface, gamma))
            std::memcpy(Uface, &U[elem * 4], 4 * sizeof(double));
    };

    for (int i = 0; i < mesh.num_interior_faces; ++i) {
        int elemL = mesh.I2E[4 * i + 0];
        int elemR = mesh.I2E[4 * i + 2];
        double fx_L, fy_L;
        interior_face_centroid(mesh, elemL, mesh.I2E[4 * i + 1], fx_L, fy_L);
        extrapolate(elemL, fx_L, fy_L, &UL_int[i * 4]);
        double fx_R, fy_R;
        interior_face_centroid(mesh, elemR, mesh.I2E[4 * i + 3], fx_R, fy_R);
        extrapolate(elemR, fx_R, fy_R, &UR_int[i * 4]);
    }

    for (int i = 0; i < mesh.num_boundary_faces; ++i) {
        int elem = mesh.B2E[3 * i + 0];
        int face = mesh.B2E[3 * i + 1];
        double fx, fy;
        boundary_face_centroid(mesh, elem, face, fx, fy);
        extrapolate(elem, fx, fy, &UL_bnd[i * 4]);
        std::memcpy(&UR_bnd[i * 4], &U[elem * 4], 4 * sizeof(double));
    }
}

void reconstruct_MP(const GriMesh& mesh, const double* U,
                    double* UL_int, double* UR_int,
                    double* UL_bnd, double* UR_bnd, double gamma) {
    std::vector<double> gradU(mesh.Ne * 8);
    compute_gradients(mesh, U, gradU.data());

    /* U_max, U_min over cell and face-neighbors (same region as BJ). */
    std::vector<double> U_max(mesh.Ne * 4), U_min(mesh.Ne * 4);
    for (int i = 0; i < mesh.Ne; ++i) {
        for (int k = 0; k < 4; ++k) {
            double u = U[i * 4 + k];
            U_max[i * 4 + k] = u;
            U_min[i * 4 + k] = u;
        }
    }
    for (int f = 0; f < mesh.num_interior_faces; ++f) {
        int elemL = mesh.I2E[4 * f + 0];
        int elemR = mesh.I2E[4 * f + 2];
        for (int k = 0; k < 4; ++k) {
            double uL = U[elemL * 4 + k], uR = U[elemR * 4 + k];
            U_max[elemL * 4 + k] = std::max(U_max[elemL * 4 + k], uR);
            U_min[elemL * 4 + k] = std::min(U_min[elemL * 4 + k], uR);
            U_max[elemR * 4 + k] = std::max(U_max[elemR * 4 + k], uL);
            U_min[elemR * 4 + k] = std::min(U_min[elemR * 4 + k], uL);
        }
    }

    /* Maximum-principle (component-wise): phi per variable per cell. */
    std::vector<double> phi_var(mesh.Ne * 4, 1.0);
    std::vector<double> cx(mesh.Ne), cy(mesh.Ne);
    for (int i = 0; i < mesh.Ne; ++i)
        cell_centroid(mesh, i, cx[i], cy[i]);

    for (int i = 0; i < mesh.Ne; ++i) {
        for (int k = 0; k < 4; ++k) {
            double phi = 1.0;
            double u_cell = U[i * 4 + k];
            for (int j = 0; j < 3; ++j) {
                double fx, fy;
                interior_face_centroid(mesh, i, j, fx, fy);
                double dx = fx - cx[i], dy = fy - cy[i];
                double delta = gradU[i * 8 + k * 2 + 0] * dx + gradU[i * 8 + k * 2 + 1] * dy;
                double u_face = u_cell + delta;
                double phi_k = 1.0;
                if (delta > 1e-14)
                    phi_k = (U_max[i * 4 + k] - u_cell) / delta;
                else if (delta < -1e-14)
                    phi_k = (U_min[i * 4 + k] - u_cell) / delta;
                phi_k = std::max(0.0, std::min(1.0, phi_k));
                phi = std::min(phi, phi_k);
            }
            phi_var[i * 4 + k] = phi;
        }
    }

    /* Limited gradient: per-variable phi. */
    std::vector<double> gradU_lim(mesh.Ne * 8);
    for (int i = 0; i < mesh.Ne; ++i) {
        for (int k = 0; k < 4; ++k) {
            double p = phi_var[i * 4 + k];
            gradU_lim[i * 8 + k * 2 + 0] = p * gradU[i * 8 + k * 2 + 0];
            gradU_lim[i * 8 + k * 2 + 1] = p * gradU[i * 8 + k * 2 + 1];
        }
    }

    /* Same extrapolation as nolimiter, using limited gradient; fallback to cell average if invalid. */
    auto extrapolate = [&](int elem, double fx, double fy, double* Uface) {
        double dx = fx - cx[elem], dy = fy - cy[elem];
        const double* g = &gradU_lim[elem * 8];
        for (int k = 0; k < 4; ++k) {
            Uface[k] = U[elem * 4 + k] + g[k * 2 + 0] * dx + g[k * 2 + 1] * dy;
        }
        if (!is_valid_state(Uface, gamma))
            std::memcpy(Uface, &U[elem * 4], 4 * sizeof(double));
    };

    for (int i = 0; i < mesh.num_interior_faces; ++i) {
        int elemL = mesh.I2E[4 * i + 0];
        int elemR = mesh.I2E[4 * i + 2];
        double fx_L, fy_L;
        interior_face_centroid(mesh, elemL, mesh.I2E[4 * i + 1], fx_L, fy_L);
        extrapolate(elemL, fx_L, fy_L, &UL_int[i * 4]);
        double fx_R, fy_R;
        interior_face_centroid(mesh, elemR, mesh.I2E[4 * i + 3], fx_R, fy_R);
        extrapolate(elemR, fx_R, fy_R, &UR_int[i * 4]);
    }

    for (int i = 0; i < mesh.num_boundary_faces; ++i) {
        int elem = mesh.B2E[3 * i + 0];
        int face = mesh.B2E[3 * i + 1];
        double fx, fy;
        boundary_face_centroid(mesh, elem, face, fx, fy);
        extrapolate(elem, fx, fy, &UL_bnd[i * 4]);
        std::memcpy(&UR_bnd[i * 4], &U[elem * 4], 4 * sizeof(double));
    }
}

void clip_cons_state(double U[4], double gamma) {
    double rho, u, v, p, c;
    consToPrim(U, gamma, rho, u, v, p, c);
    if (rho <= 0.0 || !std::isfinite(rho)) rho = 1e-12;
    if (p <= 0.0 || !std::isfinite(p)) p = 1e-12;
    primToCons(rho, u, v, p, gamma, U);
}
