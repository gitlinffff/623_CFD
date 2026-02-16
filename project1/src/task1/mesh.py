import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.spatial import Delaunay, cKDTree
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

def make_ccw(nodes, simplices):
    """
    Ensure 3 nodes of all triangles are counter-clockwise labeled.
    nodes: array of (x, y) coordinates
    simplices: array of triangle node indices (0-based)
    """
    # Get coordinates of the three vertices for all triangles
    A = nodes[simplices[:, 0]]
    B = nodes[simplices[:, 1]]
    C = nodes[simplices[:, 2]]
    
    # Compute signed area
    # Area = 0.5 * [(x2-x1)(y3-y1) - (x3-x1)(y2-y1)]
    val = (B[:, 0] - A[:, 0]) * (C[:, 1] - A[:, 1]) - \
          (C[:, 0] - A[:, 0]) * (B[:, 1] - A[:, 1])
    
    # Find triangles that are clockwise (val < 0)
    flip_mask = val < 0
    num_flipped = np.sum(flip_mask)
    print(num_flipped, "triangles flipped to ensure counter-clockwise ordering.")

    # Swap the second and third nodes for those triangles
    simplices[flip_mask, 1], simplices[flip_mask, 2] = \
        simplices[flip_mask, 2], simplices[flip_mask, 1]
    
    return simplices


def generate_unstructured_mesh(Bnodes, extent, a, output_dir):
    """
    Bnodes: boundary nodes as np.array of shape (nBNodes, ndim)
    extent: list of 4 floats [x_min, x_max, y_min, y_max]
    a: edge length of an equilateral triangle
    output_dir: directory to save output files
    """
    os.makedirs(output_dir, exist_ok=True)

    # Build initial equilateral triangle mesh
    X, Y = equilateral_tri_mesh(extent, a)

    # Check which nodes are inside the boundary
    outer_path = Path(Bnodes)
    grid_points = np.vstack((X.flatten(), Y.flatten())).T
    inside_mask = outer_path.contains_points(grid_points)
    inside_nodes = grid_points[inside_mask]

    threshold = 0.6 * a 
    # Build a KDTree for the boundary nodes for fast distance lookup
    boundary_tree = cKDTree(boundary_nodes)
    # Find the distance from every inside_node to the nearest boundary_node
    distances, _ = boundary_tree.query(inside_nodes)
    # Keep only nodes that are NOT too close to the boundary
    inside_nodes = inside_nodes[distances > threshold]

    # Plot the nodes and boundary
    plt.figure(figsize=(12, 8))
    plt.plot(inside_nodes[:,0], inside_nodes[:,1], 'ko', ms=1)
    plt.plot(Bnodes[:, 0], Bnodes[:, 1], 'o', ms=1, color='tab:orange',
            label='Boundary')
    plt.axis('equal')
    plt.savefig(os.path.join(output_dir, 'nodes.png'), dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()

    # Combine boundary and inside nodes, with boundary nodes first
    all_nodes = np.vstack((Bnodes, inside_nodes))

    # Perform Delaunay triangulation
    tri = Delaunay(all_nodes)
    simp = tri.simplices # triangle simplices (element matrix) defined by indices of all_nodes
    
    # remove triangles whose centroids are outside the boundary
    centroids = np.mean(all_nodes[simp], axis=1)
    inside_centroid_mask = outer_path.contains_points(centroids)
    simp = simp[inside_centroid_mask]

    # Ensure counter-clockwise ordering of triangle nodes
    simp = make_ccw(all_nodes, simp)

    # Convert to 1-based indexing
    simp_1based = simp + 1

    # Save nodes coordinates
    np.savetxt(os.path.join(output_dir, "boundary_nodes.txt"), Bnodes, fmt='%.4f', delimiter=',')
    np.savetxt(os.path.join(output_dir, "inside_nodes.txt"), inside_nodes, fmt='%.4f', delimiter=',')
    np.savetxt(os.path.join(output_dir, "all_nodes.txt"), all_nodes, fmt='%.4f', delimiter=',')

    # Save elements with 1-based indexing
    np.savetxt(os.path.join(output_dir, "E.txt"), simp_1based, fmt='%d', delimiter=',')

    # Plot the triangulation
    plt.figure(figsize=(12, 8))
    plt.triplot(all_nodes[:, 0], all_nodes[:, 1], simp, color='gray', lw=0.5)
    #plt.triplot(all_nodes[:, 0], all_nodes[:, 1], tri.simplices, color='gray', lw=0.5)
    plt.savefig(os.path.join(output_dir, 'tri_mesh.png'), dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    plt.close()

    # return all nodes, boundary nodes, inside nodes, and triangle elements
    return all_nodes, Bnodes, inside_nodes, simp_1based


def save_gri_file(filename, nodes, elements, Btmn_idx, PG):
    """
    Saves mesh data into the .gri format.
    nodes: np.array of shape (nNode, dim)
    elements: list of lists containing node indices for each cell
    Btmn_idx: list of boundary terminal indices, last index is the total number of boundary nodes.
    PG: dictionary of periodic groups
    """
    nNode = nodes.shape[0]
    nElemTot = elements.shape[0]
    dim = nodes.shape[-1]
    
    with open(filename, 'w') as f:
        # Header: nNode nElemTot Dim
        f.write(f"{nNode} {nElemTot} {dim}\n")
        
        # Node Coordinates
        for coord in nodes:
            f.write(" ".join(map(str, coord)) + "\n")
            
        # Boundary Groups
        nBGroup = len(Btmn_idx) - 1
        f.write(f"{nBGroup}\n")
        for i in range(nBGroup):
            nBFace = Btmn_idx[i+1] - Btmn_idx[i]
            nf = 2 # assuming line segments for boundary faces
            title = f"BGroup{i+1}"
            f.write(f"{nBFace} {nf} {title}\n")
            
            # write 1-based indices of the two nodes for each boundary face
            for j in range(Btmn_idx[i], Btmn_idx[i+1]):
                node1 = j + 1
                node2 = j + 2 if j + 1 != Btmn_idx[-1] else 1 # wrap around to first node
                f.write(f"{node1} {node2}\n")

        # Element Groups (Simplified for one group)
        nElem = nElemTot
        order = 1              # Geometry order of elements
        basis = "TriLagrange"  # geometry interpolation basis
        f.write(f"{nElem} {order} {basis}\n")
        for j in range(nElem):
            f.write(" ".join(map(str, elements[j])) + "\n")

        # Peroidic Boundary Groups
        nPG = 2
        f.write(f"{nPG} PeriodicGroup\n")
        for i in sorted(PG.keys()):
            nPGNode = len(PG[i])
            periodicity = "Translational"
            f.write(f"{nPGNode} {periodicity}\n")
            for pair in PG[i]:
                f.write(f"{int(pair[0]+1)} {int(pair[1]+1)}\n") # 1-based indexing for nodes

if __name__ == "__main__":
    output_dir = '../../output/initial_mesh4'
    os.makedirs(output_dir, exist_ok=True)

    a = 17./8. # edge length of an equilateral triangle

    # Build domain and get boundary nodes
    boundary_nodes, boundary_terminal_idx, periodic_groups, extent = build_domain(a, os.path.join(output_dir, "domain_boundary.png"))

    # Generate unstructured mesh
    all_nodes, Bnodes, inside_nodes, elem_1based = generate_unstructured_mesh(boundary_nodes, extent, a, output_dir)

    # Save mesh to .gri file
    save_gri_file(os.path.join(output_dir, 'initial_mesh.gri'), all_nodes, elem_1based,
                  boundary_terminal_idx, periodic_groups)
