import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.spatial import Delaunay
from domain import build_domain

def equilateral_tri_mesh(extent,a):
    """
    extent: list of 4 floats [x_min, x_max, y_min, y_max]
    a: edge length of an equilateral triangle
    """
    a = float(a)

    # Calculate height of the equilateral triangle
    h = a * np.sqrt(3.) / 2.
    
    # Calculate number of nodes
    nx = int(np.ceil((extent[1] - extent[0]) / a)) + 2  # number of nodes in x direction
    ny = int(np.ceil((extent[3] - extent[2]) / h)) + 2  # number of nodes in y direction

    # Generate basic grid
    x = np.arange(nx) * a
    y = np.arange(ny) * h
    X, Y = np.meshgrid(x, y)

    # Apply global offsets
    X += extent[0]
    Y += extent[2]

    # avoid boundary touching by shifting grid slightly
    X -= a/4.; Y -= 0.4 * h

    # Apply offsets to every second row to create equilateral triangles
    X[1::2, :] += a/2.

    return X, Y

def generate_unstructured_mesh(a, output_dir):
    """
    a: edge length of an equilateral triangle
    output_dir: directory to save output files
    """
    os.makedirs(output_dir, exist_ok=True)

    # Build initial equilateral triangle mesh
    boundary, extent = build_domain()
    X, Y = equilateral_tri_mesh(extent, a)

    # Check which nodes are inside the boundary
    outer_path = Path(boundary)
    grid_points = np.vstack((X.flatten(), Y.flatten())).T
    inside_mask = outer_path.contains_points(grid_points)
    inside_nodes = grid_points[inside_mask]

    # Plot the nodes and boundary
    plt.figure(figsize=(12, 8))
    plt.plot(inside_nodes[:,0], inside_nodes[:,1], 'ko', ms=1)
    plt.plot(boundary[:, 0], boundary[:, 1], 'o', ms=1, color='tab:orange',
            label='Boundary')
    plt.axis('equal')
    plt.savefig(os.path.join(output_dir, 'nodes.png'), dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()

    # Combine boundary and inside nodes
    all_nodes = np.vstack((boundary, inside_nodes))

    # Perform Delaunay triangulation
    tri = Delaunay(all_nodes)
    simp = tri.simplices # triangle simplices defined by indices of all_nodes
    
    # remove triangles whose centroids are outside the boundary
    centroids = np.mean(all_nodes[simp], axis=1)
    inside_centroid_mask = outer_path.contains_points(centroids)
    simp = simp[inside_centroid_mask]

    # Save nodes coordinates
    np.savetxt(os.path.join(output_dir, "boundary_nodes.txt"), boundary, fmt='%.4f', delimiter=',')
    np.savetxt(os.path.join(output_dir, "inside_nodes.txt"), inside_nodes, fmt='%.4f', delimiter=',')
    np.savetxt(os.path.join(output_dir, "all_nodes.txt"), all_nodes, fmt='%.4f', delimiter=',')

    # Save elements with 1-based indexing
    np.savetxt(os.path.join(output_dir, "E.txt"), simp+1, fmt='%d', delimiter=',')

    # Plot the triangulation
    plt.figure(figsize=(12, 8))
    plt.triplot(all_nodes[:, 0], all_nodes[:, 1], simp, color='gray', lw=0.5)
    #plt.triplot(all_nodes[:, 0], all_nodes[:, 1], tri.simplices, color='gray', lw=0.5)
    plt.savefig(os.path.join(output_dir, 'tri_mesh.png'), dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    plt.close()

    # return all nodes, boundary nodes, inside nodes, and triangle elements
    return all_nodes, boundary, inside_nodes, simp

if __name__ == "__main__":
    a = 17./13. # edge length of an equilateral triangle
    output_dir = '../../output/project1/'
    generate_unstructured_mesh(a, output_dir)
