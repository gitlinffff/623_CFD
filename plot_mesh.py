#!/usr/bin/env python3
"""
Simple plot of .gri mesh (displays on screen).
Usage: python plot_mesh.py [mesh.gri]
"""
import sys
import matplotlib.pyplot as plt
from pathlib import Path


def read_gri(path):
    with open(path) as f:
        lines = [L.strip() for L in f if L.strip()]
    it = iter(lines)
    Nn, Ne, dim = map(int, next(it).split())
    V = []
    for _ in range(Nn):
        x, y = map(float, next(it).split())
        V.append((x, y))
    # Skip boundary groups: "Nb 2 Name" then Nb edge lines each
    nb = int(next(it))
    for _ in range(nb):
        parts = next(it).split()
        n_edges = int(parts[0])
        for _ in range(n_edges):
            next(it)
    # Elements: "Ne 1 TriLagrange" then Ne lines (periodic groups come after in .gri)
    parts = next(it).split()
    ne = int(parts[0])
    E = []
    for _ in range(ne):
        a, b, c = map(int, next(it).split())
        E.append((a - 1, b - 1, c - 1))  # 1-based -> 0-based
    return V, E


def main():
    fpath = sys.argv[1] if len(sys.argv) > 1 else "mesh/vortex_100x100.gri"

    V, E = read_gri(fpath)
    x = [v[0] for v in V]
    y = [v[1] for v in V]
    tri = [[e[0], e[1], e[2]] for e in E]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.triplot(x, y, tri, "b-", linewidth=0.3, alpha=0.7)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Mesh: {Path(fpath).name} ({len(E)} triangles)")
    ax.grid(True, alpha=0.3)

    plt.show()


if __name__ == "__main__":
    main()
