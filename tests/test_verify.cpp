/**
 * Verify mesh read and boundary condition mapping.
 * Run from project root: ./build/tests/test_verify
 * Or: make run_test_verify (from build/)
 */
#include "bc.hpp"
#include "readgri.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

static const char* bc_type_str(BCType t) {
    switch (t) {
        case BC_WALL:    return "WALL";
        case BC_INFLOW:  return "INFLOW";
        case BC_OUTFLOW: return "OUTFLOW";
        default:         return "UNKNOWN";
    }
}

int main() {
    GriMesh mesh;
    if (!read_gri("mesh/initial_mesh_3.gri", mesh)) {
        std::printf("FAIL: Could not read mesh.\n");
        return 1;
    }
    std::printf("PASS: Mesh read successfully.\n\n");

    /* Mesh stats */
    std::printf("=== Mesh stats ===\n");
    std::printf("  Vertices: %d\n", mesh.Nn);
    std::printf("  Elements: %d\n", mesh.Ne);
    std::printf("  Interior faces: %d\n", mesh.num_interior_faces);
    std::printf("  Boundary faces: %d\n", mesh.num_boundary_faces);
    std::printf("  Boundary groups (Bname): %zu\n", mesh.Bname.size());
    for (size_t i = 0; i < mesh.Bname.size(); ++i)
        std::printf("    [%zu] %s\n", i, mesh.Bname[i].c_str());
    std::printf("\n");

    /* Topology check: sum(n*L) per element ~ 0 */
    std::vector<double> elem_sum(2 * mesh.Ne, 0.0);
    for (int i = 0; i < mesh.num_interior_faces; ++i) {
        int eL = mesh.I2E[4*i+0], eR = mesh.I2E[4*i+2];
        double nx = mesh.In[2*i], ny = mesh.In[2*i+1], L = mesh.In_len[i];
        elem_sum[eL*2+0] += nx*L; elem_sum[eL*2+1] += ny*L;
        elem_sum[eR*2+0] -= nx*L; elem_sum[eR*2+1] -= ny*L;
    }
    for (int i = 0; i < mesh.num_boundary_faces; ++i) {
        int e = mesh.B2E[3*i+0];
        double nx = mesh.Bn[2*i], ny = mesh.Bn[2*i+1], L = mesh.Bn_len[i];
        elem_sum[e*2+0] += nx*L; elem_sum[e*2+1] += ny*L;
    }
    double max_err = 0.0;
    for (int e = 0; e < mesh.Ne; ++e) {
        double m = std::sqrt(elem_sum[e*2]*elem_sum[e*2] + elem_sum[e*2+1]*elem_sum[e*2+1]);
        if (m > max_err) max_err = m;
    }
    std::printf("=== Topology check (sum n*L per element) ===\n");
    std::printf("  Max |sum n*L|: %.2e  %s\n\n", max_err,
                max_err < 1e-10 ? "PASS" : "FAIL");

    /* BC mapping: per boundary face */
    std::printf("=== Boundary condition mapping ===\n");
    int count_wall = 0, count_inflow = 0, count_outflow = 0, count_unknown = 0;
    for (int i = 0; i < mesh.num_boundary_faces; ++i) {
        int bgroup = mesh.B2E[3*i+2];
        BCType bt = get_bc_type(mesh, i);
        const char* name = mesh.Bname[bgroup-1].c_str();
        if (i < 5 || i >= mesh.num_boundary_faces - 2) {  /* sample */
            std::printf("  bface %3d: bgroup=%d (%s) -> %s\n", i, bgroup, name, bc_type_str(bt));
        } else if (i == 5) {
            std::printf("  ...\n");
        }
        if (bt == BC_WALL)    count_wall++;
        else if (bt == BC_INFLOW)  count_inflow++;
        else if (bt == BC_OUTFLOW) count_outflow++;
        else count_unknown++;
    }
    std::printf("\n  Summary: WALL=%d  INFLOW=%d  OUTFLOW=%d  UNKNOWN=%d\n",
                count_wall, count_inflow, count_outflow, count_unknown);

    if (count_unknown > 0) {
        std::printf("\nWARN: %d faces have UNKNOWN BC (unmapped BGroup).\n", count_unknown);
    }

    /* Print each boundary edge's node coordinates, grouped by BC type */
    std::printf("\n=== Boundary edges by type (node coordinates) ===\n");
    const char* types[] = {"INFLOW", "OUTFLOW", "WALL", "UNKNOWN"};
    BCType order[] = {BC_INFLOW, BC_OUTFLOW, BC_WALL, BC_UNKNOWN};
    for (int t = 0; t < 4; ++t) {
        BCType bt = order[t];
        std::printf("\n--- %s ---\n", types[t]);
        int idx = 0;
        for (int i = 0; i < mesh.num_boundary_faces; ++i) {
            if (get_bc_type(mesh, i) != bt) continue;
            int elem = mesh.B2E[3*i+0];
            int face = mesh.B2E[3*i+1];
            int n0 = mesh.E[elem*3 + face];
            int n1 = mesh.E[elem*3 + (face+1)%3];
            double x0 = mesh.V[n0*2], y0 = mesh.V[n0*2+1];
            double x1 = mesh.V[n1*2], y1 = mesh.V[n1*2+1];
            std::printf("  [%2d] (%12.6f, %12.6f) -- (%12.6f, %12.6f)\n",
                        idx, x0, y0, x1, y1);
            idx++;
        }
    }
    std::printf("\nVerification complete.\n");
    return 0;
}
