import shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import linalg
from projection import calcProjection, read_blade_segments
from projection import sizing_function_2 as sizing_function
from generate_matrices import readgri, generate_matrices

# get the blade segment coordinates for projection and sizing function evaluation
blade_seg_coords = read_blade_segments('../../data/bladeupper.txt', '../../data/bladelower.txt')

# get x coordinates of blade leading and trailing points
xL, xT = -9.461749, 9.342241

def plot_mesh_and_h_field(V, E, edge_midpoints, h_field, save_path):
    centroids = np.mean(V[E], axis=1)
    plt.figure(figsize=(12, 8))
    plt.triplot(V[:, 0], V[:, 1], E, color='gray', lw=0.5)
#    plt.scatter(edge_midpoints[:, 0], edge_midpoints[:, 1], 3, c=h_field, cmap='viridis')
#    plt.colorbar(label='h_field')

#    for i, (x, y) in enumerate(centroids):
#        plt.text(x, y, str(i), fontsize=1.5, color='red', ha='center', va='center')

#    for i, (x, y) in enumerate(V):
#        plt.text(x, y, str(i), fontsize=2., color='red', ha='center', va='center')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('Sizing Function h_field on Edge Midpoints')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(save_path, dpi=400, pad_inches=0.1, bbox_inches='tight')

def plot_mesh(V, E, title, save_path=None):
    centroids = np.mean(V[E], axis=1)
    plt.figure(figsize=(12, 8))
    plt.triplot(V[:, 0], V[:, 1], E, color='gray', lw=0.5)
    for i, (x, y) in enumerate(V):
        plt.text(x, y, str(i), fontsize=1.5, color='red', ha='center', va='center')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.show()
    #plt.savefig(save_path, dpi=400, pad_inches=0.1, bbox_inches='tight')

def local_refinement(input_gri_file, output_gri_file='local_refined.gri'):
    print("# Performing local refinement on the mesh...", flush=True)
    # Read the .gri file. There are two different functions to read the .gri file. The
    # first one is from projection.py and the second one is from generate_matrices.py.
    # We need both of them to get the edge midpoints and lengths.
    mesh = readgri(input_gri_file)
    mats = generate_matrices(input_gri_file)
   
    # extract the necessary information from the mesh and matrices
    E = mesh['E']     # 0-based indexing
    V = mesh['V']     # 0-based indexing
    I2E = mats['I2E'] # 1-based indexing
    B2E = mats['B2E'] # 1-based indexing

    # get the blade edges from B2E
#    blade_edges = []
#    for row in B2E:
#        elem = row[0] - 1
#        face = row[1]
#        bgroup = row[2]
#        n1 = E[elem, face-2]
#        n2 = E[elem, face-3]
#        if (bgroup == 2) or (bgroup == 6):
#            blade_edges.append([n1, n2])
#    blade_edges = np.array(blade_edges, dtype=int)
#    blade_edge_coords = V[blade_edges]


    # populate an array the same shape as E with -1 to flag the edges with new nodes index
    flag = np.ones(E.shape, dtype=int) * (-1)

    edge_midpoints = []
    edge_lengths = []
    h_field = []
    next_node_index = len(V)  # start indexing new nodes from here
    new_nodes = []
    boundary_edges = {int(bgid): [] for bgid in range(1, 9)} # 8 BGroups
    periodic_pairs = {1:[], 2:[]} # use a dictionary to track 2 periodic groups of node pairs

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

        edge_length = np.linalg.norm(V[L_n2] - V[L_n1])
        
        edge_midpoints.append(edge_midpoint1)
        edge_lengths.append(edge_length)

        d, xb, proj_point = calcProjection(edge_midpoint1, blade_seg_coords)
        h = sizing_function(d, xb, xL, xT)
        h_field.append(h)

        # if it is ordinary internal edge
        if (L_n1 == R_n2) and (L_n2 == R_n1):
            # flag the edge if h < edge_length
            if h < edge_length:
                # add the midpoint to new_nodes and flag it for both elements
                new_nodes.append(edge_midpoint1)
                flag[elemL, faceL-1] = next_node_index
                flag[elemR, faceR-1] = next_node_index
                next_node_index += 1
        # if the internal edge is a periodic boundary edge
        else:
            if V[L_n1, 0] <= xL: # left Periodic Group, BGroup 1 and 7
                PG = 1
                if V[L_n1, 1] > V[R_n2, 1]:
                    bgroup_L = 1; bgroup_R = 7
                else:
                    bgroup_L = 7; bgroup_R = 1
            elif V[L_n1, 0] >= xT: # right Periodic Group, BGroup 3 and 5
                PG = 2
                if V[L_n1, 1] > V[R_n2, 1]:
                    bgroup_L = 3; bgroup_R = 5
                else:
                    bgroup_L = 5; bgroup_R = 3
            else: raise ValueError("Unexpected periodic node x coordinate.")

            # track the periodic node pairs 
            periodic_pairs[PG].append((np.minimum(L_n1, R_n2), np.maximum(L_n1, R_n2)))
            periodic_pairs[PG].append((np.minimum(L_n2, R_n1), np.maximum(L_n2, R_n1)))

            # flag the edge if h < edge_length
            if h < edge_length:
                # add 2 midpoints to new_nodes and flag them for both elements
                new_nodes.append(edge_midpoint1)
                flag[elemL, faceL-1] = next_node_index
                next_node_index += 1
                new_nodes.append(edge_midpoint2)
                flag[elemR, faceR-1] = next_node_index
                next_node_index += 1
                # track the boundary edges
                boundary_edges[bgroup_L].append([L_n1, flag[elemL, faceL-1]])
                boundary_edges[bgroup_L].append([flag[elemL, faceL-1], L_n2])
                boundary_edges[bgroup_R].append([R_n1, flag[elemR, faceR-1]])
                boundary_edges[bgroup_R].append([flag[elemR, faceR-1], R_n2])
                # track the newly added periodic node pairs
                periodic_pairs[PG].append((flag[elemL, faceL-1], flag[elemR, faceR-1]))
            else:
                # track the boundary edges without adding new nodes
                boundary_edges[bgroup_L].append([L_n1, L_n2])
                boundary_edges[bgroup_R].append([R_n1, R_n2])

    # loop over boundary edges to flag those need refinement
    for row in B2E:
        elem = row[0] - 1
        face = row[1]
        bgroup = row[2]

        n1 = E[elem, face-2]
        n2 = E[elem, face-3]

        edge_midpoint = 0.5 * (V[n1] + V[n2])
        edge_length = np.linalg.norm(V[n2] - V[n1])
        edge_midpoints.append(edge_midpoint)
        edge_lengths.append(edge_length)
   
        d, xb, proj_point = calcProjection(edge_midpoint, blade_seg_coords)
        h = sizing_function(d, xb, xL, xT)
        h_field.append(h)
        if bgroup in [2, 6]: # snap the edge midpoint to the blade if it is a blade edge
            edge_midpoint = proj_point  # use spline would be better

        # flag the edge if h < edge_length
        if h < edge_length:
            new_nodes.append(edge_midpoint)
            flag[elem, face-1] = next_node_index
            next_node_index += 1
            # track the boundary edges
            boundary_edges[bgroup].append([n1, flag[elem, face-1]])
            boundary_edges[bgroup].append([flag[elem, face-1], n2])
        else:
            # track the boundary edges without adding new nodes
            boundary_edges[bgroup].append([n1, n2])

    # if no edges are flagged for refinement, stop here.
    if len(new_nodes) == 0:
        print("No edges flagged for refinement. No new nodes added.")
        return boundary_edges, periodic_pairs, 0

    edge_midpoints = np.array(edge_midpoints)
    edge_lengths = np.array(edge_lengths)
    h_field = np.array(h_field)

    # remove duplicate periodic node pairs
    periodic_pairs[1] = list(set(periodic_pairs[1]))
    periodic_pairs[2] = list(set(periodic_pairs[2]))

    plot_mesh_and_h_field(V, E, edge_midpoints, h_field, 'initial_mesh.png')

    # add the new nodes to the list of vertices
    V = np.vstack((V, np.array(new_nodes)))

    new_elements = []
    unflagged_elements = []
    # loop over elements to create new elements for the refined mesh
    for i in range(flag.shape[0]):
        if np.count_nonzero(flag[i] == -1) == 0: # 3 edges flagged
            for j in range(flag.shape[1]):
                new_elements.append([E[i, j], flag[i, j-1], flag[i, j-2]])
            new_elements.append([flag[i, 0], flag[i, 1], flag[i, 2]])
        elif np.count_nonzero(flag[i] == -1) == 1: # 2 edges flagged
            cos_min = 2.0; j_larger_angle = np.nan; j_unflagged = np.nan
            for j in range(flag.shape[1]):
                if flag[i, j] == -1:
                    j_unflagged = j
                    continue
                vec1 = V[E[i, j-1]] - V[E[i, j]]
                vec2 = V[E[i, j-2]] - V[E[i, j]]
                cos = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * \
                                            np.linalg.norm(vec2))
                if cos < cos_min:
                    cos_min = cos; j_larger_angle = j
            # create new elements, make sure node ordering is correct
            # (counter-clockwise)
            new_elements.append([E[i, j_unflagged], flag[i, j_unflagged-1], flag[i, j_unflagged-2]])
            new_elements.append([E[i, j_larger_angle], flag[i, j_unflagged-2], flag[i, j_unflagged-1]])
            new_elements.append([E[i, j_unflagged-2], E[i, j_unflagged-1], flag[i, j_larger_angle]])
        elif np.count_nonzero(flag[i] == -1) == 2: # 1 edge flagged
            j_flagged = np.where(flag[i] != -1)[0][0]
            new_elements.append([E[i, j_flagged], E[i, j_flagged-2], flag[i, j_flagged]])
            new_elements.append([E[i, j_flagged], flag[i, j_flagged], E[i, j_flagged-1]])

        else:
            unflagged_elements.append(i)

    # combine the unflagged elements and the new elements to create the new element
    # array for the refined mesh
    E = np.vstack((E[unflagged_elements], np.array(new_elements)))

    plot_mesh_and_h_field(V, E, edge_midpoints, h_field, 'refined_mesh.png')

    # output new mesh in .gri
    save_gri_file(output_gri_file, V, E+1, boundary_edges, periodic_pairs) # E has to be converted back to 1-based indexing

    return boundary_edges, periodic_pairs, len(new_nodes)
#    return E, V, boundary_edges, periodic_pairs

def smooth_mesh(input_gri_file, BG, PG, w=0.8, output_gri_file='smoothed.gri'):
    print("# Performing smoothing on the mesh...", flush=True)
    mesh = readgri(input_gri_file)
    mats = generate_matrices(input_gri_file)
   
    # extract the necessary information from the mesh and matrices
    E = mesh['E']     # 0-based indexing
    V = mesh['V']     # 0-based indexing
    I2E = mats['I2E'] # 1-based indexing
    B2E = mats['B2E'] # 1-based indexing
    periodic_pairs = PG[1] + PG[2] # combine the two periodic groups of node pairs into one list

    adj_mtx = sparse.lil_matrix((len(V),len(V)), dtype=int) # to track the adjacency of vertices for smoothing later
    # the connectivity to the adjacent node will be flagged as 1 for ordinary internal edge,
    #                                                          2 for periodic boundary edge (BGroup 1,3 5,7),
    #                                                          3 for boundary edge (BGroup 2,4,6,8)
    #                                                          4 for the internal nodes across the periodic boundary

    # loop over internal edges to track the adjacency of vertices
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

        # if it is ordinary internal edge
        if (L_n1 == R_n2) and (L_n2 == R_n1):
            # flag the connectivity as 1
            adj_mtx[L_n1, L_n2] = 1; adj_mtx[R_n1, R_n2] = 1
        # if the internal edge is a periodic boundary edge
        else:
            # flag the connectivity as 2
            adj_mtx[L_n1, L_n2] = 2; adj_mtx[L_n2, L_n1] = 2
            adj_mtx[R_n1, R_n2] = 2; adj_mtx[R_n2, R_n1] = 2

    # loop over boundary edges to flag the connectivity as 3
    for row in B2E:
        elem = row[0] - 1
        face = row[1]
        bgroup = row[2]

        # get the two node indices for this boundary edge
        n1 = E[elem, face-2]
        n2 = E[elem, face-3]

        # flag the connectivity as 3
        adj_mtx[n1, n2] = 3; adj_mtx[n2, n1] = 3

    # flag the connectivity across the periodic boundary as 4
    for n1, n2 in periodic_pairs:
        if 3 in adj_mtx.data[n1]: continue
        for j, ctype in enumerate(adj_mtx.data[n1]):
            if ctype == 1: 
                type1_neighbor = adj_mtx.rows[n1][j]
                adj_mtx[n2, type1_neighbor] = 4

        for j, ctype in enumerate(adj_mtx.data[n2]):
            if ctype == 1: 
                type1_neighbor = adj_mtx.rows[n2][j]
                adj_mtx[n1, type1_neighbor] = 4
    adj_csr = adj_mtx.tocsr()

    for pair in periodic_pairs: # used for debugging
        n1, n2 = pair
        # if V[n1,0] not equal to  V[n2,0], assert error
        assert V[n1,0] == V[n2,0], f"Periodic node pair {n1} and {n2} do not have the same x coordinate."

    # smooth the mesh several times
    for iter in range(10):
        new_V = np.copy(V)
        # Loop over all nodes
        for i in range(adj_csr.shape[0]):
            # Get the start and end indices for node i's neighbors in the CSR structure
            start = adj_csr.indptr[i]
            end = adj_csr.indptr[i + 1]
            
            # Extract the neighbor indices and the edge types (1, 2, 3 or 4)
            neighbor_indices = adj_csr.indices[start:end]
            edge_types = adj_csr.data[start:end]
            
            # Skip nodes with no neighbors (isolated nodes)
            if len(neighbor_indices) == 0:
                continue
            
            if 3 in edge_types:
                # Case 1: This is boundary node (BGroup 2,4,6,8), do not move it
                continue
            else:
                if 2 in edge_types:
                    # Case 2: This is periodic boundary node (BGroup 1,3,5,7), only move it in x direction
                    continue # do this in the next part
                else:
                    # Case 3: This is a completely interior node, move it in both x and y directions
                    new_V[i] = (1 - w) * V[i] + w * np.mean(V[neighbor_indices], axis=0)

        # loop over periodic node pairs
        for n1, n2 in periodic_pairs:
            if 3 in adj_mtx.data[n1]: continue
            # move the periodic node pair together in x direction only
            neighbor_indices = adj_mtx.rows[n2]
            new_x = (1 - w) * V[n2, 0] + w * np.mean(V[neighbor_indices, 0])
            new_V[n2, 0] = new_x; new_V[n1, 0] = new_x

        V = new_V
        #plot_mesh(V, E, f'smooth mesh iter:{iter+1}')

    # output new mesh in .gri
    save_gri_file(output_gri_file, V, E+1, BG, PG) # E has to be converted back to 1-based indexing
    
    return
    #return E, V, I2E, B2E
    

def save_gri_file(filename, nodes, elements, BG, PG):
    """
    Saves mesh data into the .gri format.
    nodes: np.array of shape (nNode, dim)
    elements: np.array containing node indices for each cell
    BG: a dictionary of boundary groups, each group is a list of boundary faces (edges)
    PG: dictionary of periodic groups, each group is a list of node pairs (tuple) that are periodic counterparts.
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
            
        # Boundary Groups of boundary faces (edges)
        nBGroup = len(BG)
        f.write(f"{nBGroup}\n")
        for BGid, BFaces in sorted(BG.items()):
            nBFace = len(BFaces)
            nf = 2 # number of nodes per boundary face (edge)
            title = f"BGroup{BGid}"
            f.write(f"{nBFace} {nf} {title}\n")
            
            # write 1-based indices of the two nodes for each boundary face
            for face in BFaces:
                node1, node2 = face[0]+1, face[1]+1 # convert to 1-based indexing
                f.write(f"{node1} {node2}\n")

        # Element Groups (Simplified for one group)
        nElem = nElemTot
        order = 1              # Geometry order of elements
        basis = "TriLagrange"  # geometry interpolation basis
        f.write(f"{nElem} {order} {basis}\n")
        for j in range(nElem):
            f.write(" ".join(map(str, elements[j])) + "\n")

        # Peroidic Boundary Groups
        nPG = len(PG)
        f.write(f"{nPG} PeriodicGroup\n")
        for PGid, pairs in sorted(PG.items()):
            nPGNode = len(pairs)
            periodicity = "Translational"
            f.write(f"{nPGNode} {periodicity}\n")
            for pair in pairs:
                node1, node2 = pair[0]+1, pair[1]+1 # convert to 1-based indexing
                f.write(f"{node1} {node2}\n")

def check_mesh_consistency(coarsemesh_filepath, refinedmesh_filepath):    
    coarsemesh_filepath = '/Users/linfel/git_repos/623_CFD_codes/project1/output/initial_mesh_3/initial_mesh.gri' 
    refinedmesh_filepath = 'refined_mesh.gri'
    E1, V1, BG, PG = local_refinement(coarsemesh_filepath)
    E2, V2, I2E, B2E = smooth_mesh(refinedmesh_filepath)
    print(E1.shape, V1.shape)
    print(E2.shape, V2.shape)

    for i, j in zip(E1.flatten(), E2.flatten()):
        if i != j:
            print("Element arrays do not match!")
    
    for i, j in zip(V1.flatten(), V2.flatten()):
        if abs(i-j) > 1e-8:
            print("Vertex arrays do not match!")

    for key, value in BG.items():
        print(f"Boundary Group {key}: {len(value)} faces")
    print("# BE: ", len(B2E))


    boundary_edges = {int(bgid): [] for bgid in set(B2E[:, 2])}

    for row in B2E:
        elem = row[0] - 1
        face = row[1]
        bgroup = row[2]

        n1 = E2[elem, face-2]
        n2 = E2[elem, face-3]
        boundary_edges[bgroup].append([n1, n2])
    for key, value in boundary_edges.items():
        print(f"Boundary Group {key}: {len(value)} faces")
    
    periodic_pairs = {1:[], 2:[]}   
    for row in I2E:
        elemL = row[0] - 1
        faceL = row[1]
        elemR = row[2] - 1
        faceR = row[3]
        
        # get the two node indices for this internal edge from the left element
        L_n1 = E2[elemL, faceL-2]
        L_n2 = E2[elemL, faceL-3]

        # get the two node indices for this internal edge from the right element
        R_n1 = E2[elemR, faceR-2]
        R_n2 = E2[elemR, faceR-3]
        if not ((L_n1 == R_n2) and (L_n2 == R_n1)):
            if V2[L_n1, 0] <= 0.: PG = 1
            elif V2[L_n1, 0] >= 0.: PG = 2
            else: raise ValueError("Unexpected periodic node x coordinate.")
            periodic_pairs[PG].append((np.minimum(L_n1, R_n2), np.maximum(L_n1, R_n2)))
            periodic_pairs[PG].append((np.minimum(L_n2, R_n1), np.maximum(L_n2, R_n1)))

    periodic_pairs[1] = list(set(periodic_pairs[1]))
    periodic_pairs[2] = list(set(periodic_pairs[2]))
    for key, value in periodic_pairs.items():
        print(f"Periodic Group {key}:")
        for pair in value:
            print(f"  Node {pair[0]} <--> Node {pair[1]}")

def main():
    src = "../../output/initial_mesh4/initial_mesh.gri"
    filepath = "mesh.gri"
    shutil.copy(src, filepath)

    iter_count = 0
    while True:
        iter_count += 1
        BG, PG, N_newnodes = local_refinement(filepath, filepath)
        if N_newnodes == 0:
            print("No new nodes added. Refinement complete.")
            break
        smooth_mesh(filepath, BG, PG, output_gri_file=filepath)
        print(f"iter {iter_count}: Number of new nodes added: {N_newnodes}")

def test():
    src = "../../output/initial_mesh4/initial_mesh.gri"
    filepath = "mesh.gri"
    shutil.copy(src, filepath)

    for i in range(10):
        BG, PG, N_newnodes = local_refinement(filepath, filepath)
        print(f"iter {i+1}: Number of new nodes added: {N_newnodes}")
        if i==15:
            break
        smooth_mesh(filepath, BG, PG, output_gri_file=filepath)

def debug():
    input_gri_file = "mesh.gri"
    output_gri_file = "mesh1.gri"
    BG, PG, N_newnodes = local_refinement(input_gri_file, output_gri_file)

if __name__ == "__main__":
    main()
