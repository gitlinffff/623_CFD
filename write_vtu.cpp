/**
 * VTU export for Paraview visualization.
 * Writes ASCII VTU (UnstructuredGrid) with cell-centered rho, u, v, p.
 */
#include "write_vtu.hpp"
#include "physics.hpp"
#include <fstream>
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
    f << "      </CellData>\n";

    f << "    </Piece>\n";
    f << "  </UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    return f.good();
}
