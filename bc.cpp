#include "bc.hpp"
#include <string>

BCType get_bc_type(const GriMesh& mesh, int bface_idx) {
    int bgroup = mesh.B2E[3 * bface_idx + 2];
    const std::string& name = mesh.Bname[bgroup - 1];
    if (name == "BGroup2" || name == "BGroup6") return BC_WALL;
    if (name == "BGroup8") return BC_INFLOW;
    if (name == "BGroup4") return BC_OUTFLOW;
    return BC_UNKNOWN;
}
