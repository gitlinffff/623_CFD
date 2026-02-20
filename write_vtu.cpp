/**
 * VTU export for Paraview visualization.
 * Writes ASCII VTU (UnstructuredGrid) with cell-centered rho, u, v, p, and wall_marker.
 * wall_marker: 0 = non-wall; +1,+2,+3 = top wall (BGroup2), 1-indexed edge; -1,-2,-3 = bottom wall (BGroup6).
 */
#include "write_vtu.hpp"
#include "physics.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

bool write_vtu(const GriMesh& mesh, const double* U, double gamma,
              const char* filepath) {
    std::ofstream f(filepath);
    if (!f.is_open())
        return false;

    /* Convert U to primitive per cell */
    std::vector<double> rho(mesh.Ne), u(mesh.Ne), v(mesh.Ne), p(mesh.Ne);
    for (int i = 0; i < mesh.Ne; ++i) {
        double c;
        consToPrim(&U[i * 4], gamma, rho[i], u[i], v[i], p[i], c);
    }

    /* wall_marker: 0 = non-wall; +1,+2,+3 = top wall (BGroup2), edge 1-indexed; -1,-2,-3 = bottom wall (BGroup6) */
    std::vector<int> wall_marker(mesh.Ne, 0);
    for (int i = 0; i < mesh.num_boundary_faces; ++i) {
        int elem = mesh.B2E[3 * i + 0];
        int face = mesh.B2E[3 * i + 1];
        int bgroup = mesh.B2E[3 * i + 2];
        const std::string& name = mesh.Bname[bgroup - 1];
        if (name == "BGroup2" || name == "BGroup6") {
            int edge_1idx = face + 1;  /* 1-indexed local edge */
            wall_marker[elem] = (name == "BGroup2") ? edge_1idx : -edge_1idx;
        }
    }

    /* Points: vertices with z=0 for 2D */
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "  <UnstructuredGrid>\n";
    f << "    <Piece NumberOfPoints=\"" << mesh.Nn << "\" NumberOfCells=\"" << mesh.Ne << "\">\n";

    f << "      <Points>\n";
    f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.Nn; ++i)
        f << "          " << mesh.V[2 * i] << " " << mesh.V[2 * i + 1] << " 0\n";
    f << "        </DataArray>\n";
    f << "      </Points>\n";

    /* Cells: triangles (VTK type 5) */
    f << "      <Cells>\n";
    f << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.Ne; ++i)
        f << "          " << mesh.E[3 * i] << " " << mesh.E[3 * i + 1] << " " << mesh.E[3 * i + 2] << "\n";
    f << "        </DataArray>\n";
    f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.Ne; ++i)
        f << "          " << (3 * (i + 1));
    f << "\n        </DataArray>\n";
    f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.Ne; ++i)
        f << "          5";
    f << "\n        </DataArray>\n";
    f << "      </Cells>\n";

    /* CellData */
    f << "      <CellData>\n";
    f << "        <DataArray type=\"Float64\" Name=\"rho\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.Ne; ++i)
        f << "          " << rho[i] << "\n";
    f << "        </DataArray>\n";
    f << "        <DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.Ne; ++i)
        f << "          " << u[i] << "\n";
    f << "        </DataArray>\n";
    f << "        <DataArray type=\"Float64\" Name=\"v\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.Ne; ++i)
        f << "          " << v[i] << "\n";
    f << "        </DataArray>\n";
    f << "        <DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.Ne; ++i)
        f << "          " << p[i] << "\n";
    f << "        </DataArray>\n";
    f << "        <DataArray type=\"Int32\" Name=\"wall_marker\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.Ne; ++i)
        f << "          " << wall_marker[i] << "\n";
    f << "        </DataArray>\n";
    f << "      </CellData>\n";

    f << "    </Piece>\n";
    f << "  </UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    return f.good();
}

namespace {
bool parse_dataarray(const std::string& content, const char* name, std::vector<double>& out) {
    std::string tag = "Name=\"" + std::string(name) + "\"";
    size_t pos = content.find(tag);
    if (pos == std::string::npos) return false;
    size_t start = content.find('>', pos);
    if (start == std::string::npos) return false;
    start++;
    size_t end = content.find("</DataArray>", start);
    if (end == std::string::npos) return false;
    std::string block = content.substr(start, end - start);
    out.clear();
    std::istringstream ss(block);
    double v;
    while (ss >> v) out.push_back(v);
    return !out.empty();
}
}  // namespace

bool read_vtu(const GriMesh& mesh, const char* filepath, double gamma, double* U) {
    std::ifstream f(filepath);
    if (!f.is_open()) return false;
    std::stringstream buf;
    buf << f.rdbuf();
    std::string content = buf.str();

    std::vector<double> rho, u, v, p;
    if (!parse_dataarray(content, "rho", rho) || !parse_dataarray(content, "u", u) ||
        !parse_dataarray(content, "v", v) || !parse_dataarray(content, "p", p))
        return false;
    if ((int)rho.size() != mesh.Ne || (int)u.size() != mesh.Ne ||
        (int)v.size() != mesh.Ne || (int)p.size() != mesh.Ne)
        return false;

    for (int i = 0; i < mesh.Ne; ++i) {
        double Ui[4];
        primToCons(rho[i], u[i], v[i], p[i], gamma, Ui);
        U[i * 4 + 0] = Ui[0];
        U[i * 4 + 1] = Ui[1];
        U[i * 4 + 2] = Ui[2];
        U[i * 4 + 3] = Ui[3];
    }
    return true;
}
