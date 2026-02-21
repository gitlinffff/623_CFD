#!/usr/bin/env python3
"""
Create a Paraview PVD collection file from VTU files in a directory.
Usage: python write_pvd.py [directory] [output.pvd]
  directory: path containing t=*.vtu files (default: data)
  output.pvd: output PVD filename (default: solution.pvd)
"""
import os
import re
import sys


def write_pvd(directory="data", pvd_name="solution.pvd"):
    """
    Scan directory for t=*.vtu files and create a PVD collection.
    """
    if not os.path.isdir(directory):
        print(f"Error: directory not found: {directory}")
        return False

    pattern = re.compile(r"t=([0-9eE\.\+\-]+)\.vtu")
    entries = []
    for fname in os.listdir(directory):
        match = pattern.match(fname)
        if match:
            t = float(match.group(1))
            entries.append((t, fname))

    if not entries:
        print(f"Error: no t=*.vtu files found in {directory}")
        return False

    entries.sort(key=lambda x: x[0])
    pvd_path = os.path.join(directory, pvd_name)

    with open(pvd_path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        f.write("  <Collection>\n")
        for t, fname in entries:
            f.write(f'    <DataSet timestep="{t}" file="{fname}"/>\n')
        f.write("  </Collection>\n")
        f.write("</VTKFile>\n")

    print(f"[PVD] Written: {pvd_path}")
    return True


if __name__ == "__main__":
    directory = sys.argv[1] if len(sys.argv) > 1 else "data"
    pvd_name = sys.argv[2] if len(sys.argv) > 2 else "solution.pvd"
    write_pvd(directory, pvd_name)
