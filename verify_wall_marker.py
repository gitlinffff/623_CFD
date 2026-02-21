#!/usr/bin/env python3
"""
Verify wall_marker in VTU output and optionally plot wall boundaries.
Usage: python verify_wall_marker.py [path/to/file.vtu] [--plot]
Default: data/add_wall_marker.vtu
"""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def parse_vtu_full(filepath):
    """Parse VTU: points, cells (connectivity), and cell data."""
    tree = ET.parse(filepath)
    root = tree.getroot()

    # Points
    points_section = root.find(".//Points")
    points_elem = points_section.find("DataArray") if points_section is not None else None
    if points_elem is None:
        points_elem = root.find(".//Points/DataArray")
    points_text = (points_elem.text or "").strip() if points_elem is not None else ""
    pts = [float(x) for x in points_text.split()]
    points = [(pts[i], pts[i + 1]) for i in range(0, len(pts), 3)]

    # Cells connectivity
    for da in root.iter("DataArray"):
        if da.get("Name") == "connectivity":
            conn_text = (da.text or "").strip()
            conn = [int(x) for x in conn_text.split()]
            break
    else:
        conn = []

    # CellData
    celldata = root.find(".//CellData")
    arrays = {}
    if celldata is not None:
        for da in celldata.findall("DataArray"):
            name = da.get("Name")
            if name:
                text = (da.text or "").strip()
                arrays[name] = [x.strip() for x in text.split() if x.strip()]

    piece = root.find(".//Piece")
    n_cells = int(piece.get("NumberOfCells", 0)) if piece is not None else len(conn) // 3
    return points, conn, arrays, n_cells


def parse_vtu_celldata(filepath):
    """Parse VTU and extract CellData arrays."""
    tree = ET.parse(filepath)
    root = tree.getroot()
    ns = {"vtk": "http://www.vtk.org/schema/VTK"}
    # VTU may not use namespace
    celldata = root.find(".//CellData")
    if celldata is None:
        celldata = root.find(".//vtk:CellData", ns)
    if celldata is None:
        return None, {}

    arrays = {}
    for da in celldata.findall("DataArray"):
        name = da.get("Name")
        if name is None:
            da = celldata.find("vtk:DataArray", ns)
            if da is not None:
                name = da.get("Name")
        if name:
            text = (da.text or "").strip()
            values = [x.strip() for x in text.split() if x.strip()]
            arrays[name] = values

    # Get number of cells from Piece
    piece = root.find(".//Piece")
    if piece is None:
        piece = root.find(".//vtk:Piece", ns)
    n_cells = int(piece.get("NumberOfCells", 0)) if piece is not None else 0
    return n_cells, arrays


def verify_wall_marker(filepath):
    """Verify wall_marker array in VTU."""
    path = Path(filepath)
    if not path.exists():
        print(f"Error: File not found: {path}")
        return False

    n_cells, arrays = parse_vtu_celldata(path)
    if not arrays:
        print(f"Error: No CellData found in {path}")
        return False

    if "wall_marker" not in arrays:
        print(f"Error: wall_marker not found in VTU. Available: {list(arrays.keys())}")
        return False

    wm = arrays["wall_marker"]
    if len(wm) != n_cells:
        print(f"Warning: wall_marker length {len(wm)} != NumberOfCells {n_cells}")

    # Parse as int
    valid = {-3, -2, -1, 0, 1, 2, 3}
    counts = {v: 0 for v in valid}
    invalid = []

    for i, s in enumerate(wm):
        try:
            v = int(float(s))
            if v in valid:
                counts[v] += 1
            else:
                invalid.append((i, v))
        except ValueError:
            invalid.append((i, s))

    # Report
    print(f"=== wall_marker verification: {path} ===\n")
    print(f"Total cells: {n_cells}")
    print(f"wall_marker length: {len(wm)}\n")

    print("Value distribution:")
    for v in sorted(valid):
        c = counts[v]
        label = "non-wall" if v == 0 else ("top wall" if v > 0 else "bottom wall")
        print(f"  {v:3d}: {c:5d}  ({label})")

    n_wall = sum(counts[v] for v in valid if v != 0)
    print(f"\nWall cells (non-zero): {n_wall}")

    if invalid:
        print(f"\nInvalid values: {len(invalid)}")
        for i, v in invalid[:10]:
            print(f"  cell {i}: {v}")
        if len(invalid) > 10:
            print(f"  ... and {len(invalid) - 10} more")
        return False

    # Sanity: wall cells should have edge 1, 2, or 3
    for v in [1, 2, 3, -1, -2, -3]:
        if counts[v] > 0:
            print(f"\n  Edge |{v}| (1-indexed): {counts[v]} cells")

    print("\n=== PASS ===")
    return True


def plot_wall_boundaries(filepath, outpath=None):
    """Plot wall boundary edges from VTU. Top wall=blue, bottom wall=red."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Error: matplotlib required for --plot. Install with: pip install matplotlib")
        return False

    path = Path(filepath)
    if not path.exists():
        print(f"Error: File not found: {path}")
        return False

    points, conn, arrays, n_cells = parse_vtu_full(path)
    if "wall_marker" not in arrays:
        print("Error: wall_marker not found in VTU")
        return False

    wm = [int(float(x)) for x in arrays["wall_marker"]]
    if len(conn) != n_cells * 3:
        print("Warning: connectivity length mismatch")

    # For each wall cell: face f (0-indexed) = abs(wm)-1.
    # In readgri, face f connects tri[f] and tri[(f+1)%3], i.e. conn[3*i+f] and conn[3*i+(f+1)%3]
    top_edges = []  # list of ((x0,y0),(x1,y1))
    bottom_edges = []

    for i in range(min(n_cells, len(wm))):
        val = wm[i]
        if val == 0:
            continue
        f = abs(val) - 1  # 0-indexed face
        v0 = conn[3 * i + f]
        v1 = conn[3 * i + (f + 1) % 3]
        p0 = points[v0]
        p1 = points[v1]
        seg = (p0, p1)
        if val > 0:
            top_edges.append(seg)
        else:
            bottom_edges.append(seg)

    fig, ax = plt.subplots(figsize=(10, 6))
    for (p0, p1) in top_edges:
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "b-", linewidth=1.5)
    for (p0, p1) in bottom_edges:
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "r-", linewidth=1.5)

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Wall boundaries (blue=top/BGroup2, red=bottom/BGroup6)")
    ax.grid(True, alpha=0.3)

    if outpath is None:
        outpath = path.with_suffix(".wall_plot.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Saved: {outpath}")
    plt.close()
    return True


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--plot"]
    do_plot = "--plot" in sys.argv
    default = Path(__file__).parent / "data" / "add_wall_marker.vtu"
    filepath = args[0] if args else str(default)

    ok = verify_wall_marker(filepath)
    if ok and do_plot:
        plot_wall_boundaries(filepath)
    sys.exit(0 if ok else 1)
