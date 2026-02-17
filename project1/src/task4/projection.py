import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy.interpolate import CubicSpline

# create blade geometry using cubic splines
# Coordinates start and end at the trailing tip
coords = np.loadtxt('../../data/blade.txt')
# Calculate the distance between consecutive points
ds = np.sqrt(np.sum(np.diff(coords, axis=0)**2, axis=1))
# Create the parameter s (cumulative distance) starting at 0
s = np.insert(np.cumsum(ds), 0, 0)
# Create splines for x and y
spline_x = CubicSpline(s, coords[:, 0], bc_type='periodic')
spline_y = CubicSpline(s, coords[:, 1], bc_type='periodic')


def read_gri(filename):
    with open(filename, 'r') as f:
        tokens = f.read().split()
    
    it = iter(tokens)
    try:
        n_node = int(next(it))
        n_elem_tot = int(next(it))
        dim = int(next(it))
        
        nodes = np.zeros((n_node, 2))
        for i in range(n_node):
            nodes[i, 0] = float(next(it))
            nodes[i, 1] = float(next(it))
            if dim == 3: next(it)

        n_bgroup = int(next(it))
        boundaries = {}
        for _ in range(n_bgroup):
            n_bface = int(next(it))
            nodes_per_face = int(next(it))
            title = next(it).strip('"')
            
            faces = []
            for _ in range(n_bface):
                face = [int(next(it)) - 1 for _ in range(nodes_per_face)]
                faces.append(face)
            boundaries[title] = faces

        elements = []
        elems_read = 0
        while elems_read < n_elem_tot:
            n_elem = int(next(it))
            order = int(next(it))
            basis = next(it).strip('"')
            
            nodes_per_elem = 3 if 'tri' in basis.lower() else 4
            
            for _ in range(n_elem):
                elem = [int(next(it)) - 1 for _ in range(nodes_per_elem)]
                elements.append(elem)
            elems_read += n_elem
            
        return nodes, boundaries, elements
    except StopIteration:
        return None, None, None

def get_blade_segments(nodes, boundaries):
    blade_segments = []
    target_groups = ["BGroup2", "BGroup6"] 
    
    for title, faces in boundaries.items():
        if title in target_groups:
            for face in faces:
                p1 = nodes[face[0]]
                p2 = nodes[face[1]]
                blade_segments.append((p1, p2))
                
    return np.array(blade_segments)


def calcDist(P, A, B):
    AB = B - A
    AP = P - A
    len_AB_sq = np.dot(AB, AB)
    if len_AB_sq == 0:
        return np.linalg.norm(AP), A
    t = np.dot(AP, AB) / len_AB_sq
    t_clamped = np.clip(t, 0.0, 1.0)
    closest_point = A + t_clamped * AB
    distance = np.linalg.norm(P - closest_point)
    return distance, closest_point, np.sqrt(len_AB_sq)

def calcProjection_multinodes(nodes, blade_segments):
    nNodes = len(nodes)
    d = np.zeros(nNodes)
    xb = np.zeros(nNodes)
    proj_points = np.zeros_like(nodes)
    for i in range(nNodes):
        P = nodes[i]
        min_dist = float('inf')
        best_proj = P
        for seg in blade_segments:
            dist, proj = calcDist(P, seg[0], seg[1])
            if dist < min_dist:
                min_dist = dist
                best_proj = proj
        d[i] = min_dist
        xb[i] = best_proj[0]
        proj_points[i] = best_proj
    return d, xb, proj_points

def read_blade_segments(bladeupper_filepath, bladelower_filepath):
    bladeupper_coords = np.loadtxt(bladeupper_filepath)
    bladelower_coords = np.loadtxt(bladelower_filepath)

    # shift the original upper blade down by 18 units to be the lower boundary
    bladeupper_coords[:, 1] -= 18

    N_upper = len(bladeupper_coords); N_lower = len(bladelower_coords)
   
    # create segments by connecting consecutive points
    a = np.linspace(0, N_upper-2, num=N_upper-1, dtype=int)
    b = np.linspace(1, N_upper-1, num=N_upper-1, dtype=int)
    c = np.vstack((a, b)).T
    upper_seg = bladeupper_coords[c]

    a = np.linspace(0, N_lower-2, num=N_lower-1, dtype=int)
    b = np.linspace(1, N_lower-1, num=N_lower-1, dtype=int)
    c = np.vstack((a, b)).T
    lower_seg = bladelower_coords[c]

    blade_seg = np.vstack((upper_seg, lower_seg))
    
    return blade_seg


def calcProjection(P, blade_segments):
    # find the closest point on the blade segments to P
    best_proj = np.array([np.nan, np.nan])
    min_dist = float('inf')
    lengths = []
    i_best = None
    for i, seg in enumerate(blade_segments):
        dist, proj, seg_len = calcDist(P, seg[0], seg[1])
        lengths.append(seg_len)
        if dist < min_dist:
            min_dist = dist
            best_proj = proj
            i_best = i
    # snap the projection point from the blade segments to the spline curve
    # sum the arc length all the way to the best projection point
    arc_len = sum(lengths[:i_best]) + np.linalg.norm(best_proj - blade_segments[i_best][0])
    xs, ys = spline_x(arc_len), spline_y(arc_len)
    # shift the y-coordinate down by 18 units to align with the lower boundary
    if abs(best_proj[1] - ys) > 17.: ys -= 18
    spl_proj = np.array([xs, ys])
    min_dist = np.linalg.norm(P - spl_proj)
    return min_dist, xs, spl_proj

# def sizing_function(d, xb):
#     hmax = 0.15
#     hmin = 0.1
#     x_le = np.min(xb)
#     x_te = np.max(xb)
#     growth_rate = 0.1
#     chord = x_te - x_le
#     dist_to_edge = np.minimum(np.abs(xb - x_le), np.abs(xb - x_te))
#     h_base = hmin + (hmax - hmin) * np.clip(
#         dist_to_edge / (0.2 * chord), 0.0, 1.0
#     )
#     h = h_base + growth_rate * d
#     return h

def sizing_function_1(d, xb, xL, xT):
    A = 0.47
    B = 0.47
    alpha = 0.5
    delta = 0.1
    hmin = 0.01
    
    term1 = A * (1 + np.exp(-((xL-xT)/delta)**2) - np.exp(-((xb-xL)/delta)**2) - np.exp(-((xb-xT)/delta)**2)) + np.sqrt(hmin)
    term2 = B * d**alpha + np.sqrt(hmin)

    return term1 * term2

def sizing_function_2(d, xb, xL, xT):
    hmin = 0.6
    sigma = 12
    delta = 11

    h = hmin * np.exp(-(xb-xL)*(xb-xT)/sigma**2) * np.exp(d/delta)
    return h

def sizing_function(d, xb):
    hmax = 0.12
    hmin = 0.1
    growth_rate = 0.4

    x_le = xb.min()
    x_te = xb.max()
    chord = x_te - x_le

    dist_to_edge = np.minimum(np.abs(xb - x_le), np.abs(xb - x_te))
    s = np.clip(dist_to_edge / (0.25 * chord), 0.0, 1.0)

    # cosine smoothstep
    smooth = 0.5 * (1.0 - np.cos(np.pi * s))
    h_base = hmin + (hmax - hmin) * smooth

    h = h_base + growth_rate * d
    return h

def plot_results(nodes, elements, d, h, proj_points, blade_segments):
    tri_elements = [e[:3] for e in elements]
    triangulation = mtri.Triangulation(nodes[:, 0], nodes[:, 1], tri_elements)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6,9))
    
    ax1.set_title("Projection Visualization (Sample Points)", fontsize=14)
    ax1.set_aspect('equal')
    for seg in blade_segments:
        ax1.plot([seg[0][0], seg[1][0]], [seg[0][1], seg[1][1]], 'k-', linewidth=2)
    
    indices = np.random.choice(len(nodes), 100, replace=False) 
    for idx in indices:
        P = nodes[idx]
        Proj = proj_points[idx]
        ax1.plot([P[0], Proj[0]], [P[1], Proj[1]], 'r--', alpha=0.3, linewidth=0.5)
        ax1.plot(P[0], P[1], 'b.', markersize=2)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.grid(True, linestyle=':', alpha=0.5)

    ax2.set_title("Mesh Sizing Function $h(d, x_b)$ Contours", fontsize=14)
    ax2.set_aspect('equal')
    for seg in blade_segments:
        ax2.plot([seg[0][0], seg[1][0]], [seg[0][1], seg[1][1]], 'k-', linewidth=1)
        
    levels = np.linspace(np.min(h), np.max(h), 20)
    contour = ax2.tricontourf(triangulation, h, levels=levels, cmap='viridis')
    cbar = plt.colorbar(contour, ax=ax2)
    cbar.set_label('Target Edge Length $h$', fontsize=12)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.tight_layout()
    plt.show()

def plot_and_export_results(nodes, elements, d, h, proj_points, blade_segments, output_dir="plots"):
    # Create directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    tri_elements = [e[:3] for e in elements]
    triangulation = mtri.Triangulation(nodes[:, 0], nodes[:, 1], tri_elements)
    
    # --- Figure 1: Projection Visualization ---
    plt.figure(figsize=(10, 7))
    plt.title("Projection Visualization", fontsize=14)
    plt.gca().set_aspect('equal')
    
    plt.plot([], [], 'k-', linewidth=1, label='Blade Geometry')
    for seg in blade_segments:
        plt.plot([seg[0][0], seg[1][0]], [seg[0][1], seg[1][1]], 'k-', linewidth=2)
    
    indices = np.random.choice(len(nodes), min(100, len(nodes)), replace=False) 
    first_sample = True
    for idx in indices:
        P = nodes[idx]
        Proj = proj_points[idx]
        if first_sample:
            plt.plot([P[0], Proj[0]], [P[1], Proj[1]], 'r--', alpha=0.5, linewidth=1, label='Projection Line')
            plt.plot(P[0], P[1], 'b.', markersize=4, label='Random Sample Point')
            first_sample = False
        else:
            plt.plot([P[0], Proj[0]], [P[1], Proj[1]], 'r--', alpha=0.5, linewidth=1)
            plt.plot(P[0], P[1], 'b.', markersize=4)
        
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend(loc='upper right', fontsize=10, frameon=True, shadow=True)
    
    # Export Figure 1
    plt.savefig(f"{output_dir}/projection_viz_hd.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/projection_viz_vector.pdf", bbox_inches='tight')
    plt.show()

    # --- Figure 2: Optimized Mesh Sizing Contours ---
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_title("Mesh Sizing Function $h(d, x_b)$ Distribution", fontsize=16, pad=20)
    ax.set_aspect('equal')
    
    levels = np.linspace(np.min(h), np.max(h), 50)
    contour = ax.tricontourf(triangulation, h, levels=levels, cmap='viridis')
    
    for seg in blade_segments:
        ax.plot([seg[0][0], seg[1][0]], [seg[0][1], seg[1][1]], 'w-', linewidth=2, alpha=0.8)
    
    cbar = fig.colorbar(contour, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label('Target Edge Length $h$ (mm)', fontsize=13, labelpad=10)
    
    ax.set_xlabel('x (mm)', fontsize=12)
    ax.set_ylabel('y (mm)', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.tick_params(direction='in', top=True, right=True)

    # Export Figure 2
    plt.tight_layout()
    plt.savefig(f"{output_dir}/sizing_contours_hd.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/sizing_contours_vector.pdf", bbox_inches='tight')
    plt.show()

    print(f"Success! High-resolution plots saved to the '{output_dir}' folder.")
def write_sizing_to_file(filename, nodes, h_field):
    with open(filename, 'w') as f:
        f.write("# node_id   x   y   h\n")
        for i, (xy, h) in enumerate(zip(nodes, h_field)):
            f.write(f"{i} {xy[0]:.8e} {xy[1]:.8e} {h:.8e}\n")


def get_h_field(nodes, boundaries):
    blade_segments = get_blade_segments(nodes, boundaries)
    d_field, xb_field, _ = calcProjection(nodes, blade_segments)
    h_field = sizing_function(d_field, xb_field)
    return h_field

if __name__ == "__main__":
    filepath = '/Users/linfel/git_repos/623_CFD_codes/project1/output/initial_mesh_3/initial_mesh.gri' 
    nodes, boundaries, elements = read_gri(filepath)
    
    if nodes is not None:
        blade_segments = get_blade_segments(nodes, boundaries)
        print(f"Extracted {len(blade_segments)} blade segments from BGroup2 and BGroup6.")
        
        d_field, xb_field, proj_points = calcProjection_multinodes(nodes, blade_segments)
        h_field = sizing_function(d_field, xb_field)
        # target_path = "/Users/zhangzhouyu/Desktop/26winter/AE623/Project-1/623_CFD-Task-5-and-Task-6/h.dat"

        # write_sizing_to_file(target_path, nodes, h_field)
        plot_and_export_results(nodes, elements, d_field, h_field, proj_points, blade_segments)
