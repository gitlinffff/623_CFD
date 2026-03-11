import numpy as np
import os
import matplotlib.pyplot as plt
from projection import calcProjection, read_blade_segments
from projection import sizing_function_2 as sizing_function
from generate_matrices import readgri, generate_matrices
from scipy.interpolate import griddata

# get the blade segment coordinates for projection and sizing function evaluation
blade_seg_coords = read_blade_segments('../../data/bladeupper.txt', '../../data/bladelower.txt')

# get x coordinates of blade leading and trailing points
xL, xT = -9.461749, 9.342241

def plot_h_field_tri(input_gri_file, output_dir):
    mesh = readgri(input_gri_file)
   
    E = mesh['E']     # 0-based indexing
    V = mesh['V']     # 0-based indexing

    h_values = []
    for node in V:
        d, xb, proj_point = calcProjection(node, blade_seg_coords)
        h = sizing_function(d, xb, xL, xT)
        h_values.append(h)
    h_values = np.array(h_values)

    plt.figure(figsize=(12, 8))
    
    # 1. Plot the filled contours using the mesh connectivity
    # V[:, 0] is x, V[:, 1] is y, E is the (Ne, 3) connectivity
    tc = plt.tricontourf(V[:, 0], V[:, 1], E, h_values, levels=50, cmap='viridis')
    
    # 2. Add a colorbar
    plt.colorbar(tc, label='Sizing Function h_field')
    
    # Optional: Overlay the mesh lines in faint white to see the elements
    # plt.triplot(V[:, 0], V[:, 1], E, color='white', lw=0.2, alpha=0.3)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('Contour of h_field (Triangulation Method)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(os.path.join(output_dir, 'sizingfunction_test.png'), dpi=400, bbox_inches='tight')
    plt.close()

def plot_mesh_and_h_field(V, E, edge_midpoints, h_field, save_path):
    centroids = np.mean(V[E], axis=1)
    plt.figure(figsize=(12, 8))
#    plt.triplot(V[:, 0], V[:, 1], E, color='gray', lw=0.5)
    plt.scatter(edge_midpoints[:, 0], edge_midpoints[:, 1], 3, c=h_field, cmap='viridis')
    plt.colorbar(label='h_field')

#    for i, (x, y) in enumerate(centroids):
#        plt.text(x, y, str(i), fontsize=1.5, color='red', ha='center', va='center')

#    for i, (x, y) in enumerate(V):
#        plt.text(x, y, str(i), fontsize=2., color='red', ha='center', va='center')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('Sizing Function h_field')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(save_path, dpi=400, pad_inches=0.1, bbox_inches='tight')
    plt.close()

def plot_sizing_function(input_gri_file, output_dir):
    mesh = readgri(input_gri_file)
    mats = generate_matrices(input_gri_file)
   
    # extract the necessary information from the mesh and matrices
    E = mesh['E']     # 0-based indexing
    V = mesh['V']     # 0-based indexing
    I2E = mats['I2E'] # 1-based indexing
    B2E = mats['B2E'] # 1-based indexing

    # populate an array the same shape as E with -1 to flag the edges with new nodes index
    flag = np.ones(E.shape, dtype=int) * (-1)

    edge_midpoints = []
    h_field = []

    # loop over internal edges to flag those need refinement
    for row in I2E:
        elemL = row[0] - 1
        faceL = row[1]
        elemR = row[2] - 1
        faceR = row[3]
        
        # get the two node indices for this internal edge from the left element
        L_n1 = E[elemL, faceL-2]
        L_n2 = E[elemL, faceL-3]

        # get the two node indices for this internal edge from the right element
        R_n1 = E[elemR, faceR-2]
        R_n2 = E[elemR, faceR-3]

        # The two midpoints should be the same, unless this internal edge is actually a periodic boundary edge
        edge_midpoint1 = 0.5 * (V[L_n1] + V[L_n2])
        edge_midpoint2 = 0.5 * (V[R_n1] + V[R_n2])
        edge_midpoints.append(edge_midpoint1)

        d, xb, proj_point = calcProjection(edge_midpoint1, blade_seg_coords)
        h = sizing_function(d, xb, xL, xT)
        h_field.append(h)
    # loop over boundary edges to flag those need refinement
    for row in B2E:
        elem = row[0] - 1
        face = row[1]
        bgroup = row[2]

        n1 = E[elem, face-2]
        n2 = E[elem, face-3]

        edge_midpoint = 0.5 * (V[n1] + V[n2])
        edge_midpoints.append(edge_midpoint)
   
        d, xb, proj_point = calcProjection(edge_midpoint, blade_seg_coords)
        h = sizing_function(d, xb, xL, xT)
        h_field.append(h)
    
    edge_midpoints = np.array(edge_midpoints)
    h_field = np.array(h_field)

    plot_mesh_and_h_field(V, E, edge_midpoints, h_field, os.path.join(output_dir,
                                                                      'sizingfunction_test.png'))

if __name__ == "__main__":
    input_gri_file = "/Users/linfel/git_repos/623_CFD_codes/project1/output/refined_mesh_ver3/global_refine_2.gri"
    output_dir = "/Users/linfel/git_repos/623_CFD_codes/project1/output/refined_mesh_ver3"
    #plot_sizing_function(input_gri_file, output_dir)
    plot_h_field_tri(input_gri_file, output_dir)
