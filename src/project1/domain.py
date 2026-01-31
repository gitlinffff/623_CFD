import numpy as np
import matplotlib.pyplot as plt

def resample_by_distance(coords, target_distance):
    # Compute distances between consecutive points
    ds = np.sqrt(np.sum(np.diff(coords, axis=0)**2, axis=1))

    # Compute cumulative distance (arc-length) along the curve
    s = np.insert(np.cumsum(ds), 0, 0)
    
    buckets = np.floor(s / target_distance)

    # Find the index of the first occurrence of each unique bucket value
    _, first_indices = np.unique(buckets, return_index=True)
    
    # Ensure the last point of the original coordinates is included
    if first_indices[-1] != (len(coords) - 1):
        first_indices = np.append(first_indices, len(coords) - 1)

    resampled_coords = coords[first_indices]

    return resampled_coords

def build_domain(a, figname=""):

    file_bladelower = "../../data/project1/bladelower.txt"
    file_bladeupper = "../../data/project1/bladeupper.txt"
    bladelower = np.loadtxt(file_bladelower)
    bladeupper = np.loadtxt(file_bladeupper)
    bladeupper[:,1] -= 18.0  # Adjust y-coordinates of the upper blade surface

    # x and y coordinates of the blade corners
    Ax, Ay = bladelower[0,0], bladelower[0,1]
    By = bladeupper[-1,1]
    Cy = bladeupper[0,1]
    Dx, Dy = bladelower[-1,0], bladelower[-1,1]

    # x grid for the periodic boundary group 1 and group 2
    nBNodes = 14
    PG1_up_xgrid = np.linspace(Ax-17.0, Ax, nBNodes) # x-nodes on periodic boundaries must match
    PG1_low_xgrid = PG1_up_xgrid[::-1]                    # reverse
    PG2_up_xgrid = np.linspace(Dx, Dx+17.0, nBNodes) # x-nodes on periodic boundaries must match
    PG2_low_xgrid = PG2_up_xgrid[::-1]                    # reverse

    # y grid for the periodic boundary group 1 and group 2
    PG1_up_ygrid  = np.full(nBNodes, Ay)
    PG1_low_ygrid = np.full(nBNodes, By)
    PG2_up_ygrid  = np.full(nBNodes, Dy)
    PG2_low_ygrid = np.full(nBNodes, Cy)

    # inflow and outflow boundaries
    nBNodes_inout = 15
    inflow_xgrid = np.full(nBNodes_inout, Ax-17.0)
    inflow_ygrid = np.linspace(By, Ay, nBNodes_inout)
    outflow_xgrid = np.full(nBNodes_inout, Dx+17.0)
    outflow_ygrid = np.linspace(Dy, Cy, nBNodes_inout)

    # Resample blade nodes
    bladelower_coarse = resample_by_distance(bladelower, 0.85*a)
    bladeupper_coarse = resample_by_distance(bladeupper, 0.85*a)

    # Construct a dictionary for boundary groups
    Bgroups = {}
    Bgroups[1] = np.column_stack((PG1_up_xgrid, PG1_up_ygrid))[:-1, :]
    Bgroups[2] = bladelower_coarse[:-1, :]
    Bgroups[3] = np.column_stack((PG2_up_xgrid, PG2_up_ygrid))[:-1, :]
    Bgroups[4] = np.column_stack((outflow_xgrid, outflow_ygrid))[:-1, :]
    Bgroups[5] = np.column_stack((PG2_low_xgrid, PG2_low_ygrid))[:-1, :]
    Bgroups[6] = bladeupper_coarse[:-1, :]
    Bgroups[7] = np.column_stack((PG1_low_xgrid, PG1_low_ygrid))[:-1, :]
    Bgroups[8] = np.column_stack((inflow_xgrid, inflow_ygrid))[:-1, :]

    # combine all boundaries into a single array
    Bnodes = np.vstack([Bgroups[i] for i in range(1, len(Bgroups) + 1)])

    # find the extents of the domain
    x_min, x_max = np.min(Bnodes[:, 0]), np.max(Bnodes[:, 0])
    y_min, y_max = np.min(Bnodes[:, 1]), np.max(Bnodes[:, 1])
    extent = [x_min, x_max, y_min, y_max]

    # Calculate terminal indices for boundary groups
    # The last value represents the total number of boundary nodes
    lengths = np.array([len(Bgroups[i]) for i in range(1, len(Bgroups) + 1)])
    Btmn_idx = np.cumsum(np.insert(lengths, 0, 0))
    print("Boundary Group Terminal Indices:", Btmn_idx)
    print("The last value represents the total number of boundary nodes.")

    # Periodic boundary groups
    PG = {}
    PG1_a = np.linspace(Btmn_idx[0], Btmn_idx[1], Btmn_idx[1]-Btmn_idx[0]+1)
    PG1_b = np.linspace(Btmn_idx[7], Btmn_idx[6], Btmn_idx[7]-Btmn_idx[6]+1)
    PG[1] = np.column_stack((PG1_a, PG1_b))
    PG2_a = np.linspace(Btmn_idx[2], Btmn_idx[3], Btmn_idx[3]-Btmn_idx[2]+1)
    PG2_b = np.linspace(Btmn_idx[5], Btmn_idx[4], Btmn_idx[5]-Btmn_idx[4]+1)
    PG[2] = np.column_stack((PG2_a, PG2_b))

    # Plotting the boundary nodes if a figure name is provided
    if (not not figname): 
        # Create the plot
        plt.figure(figsize=(12, 8))

        # Periodic Boundaries
        plt.plot(Bgroups[1][:, 0], Bgroups[1][:, 1], '.', ms=3, color='tab:green',  label='B1: PG1 Top')
        plt.plot(Bgroups[3][:, 0], Bgroups[3][:, 1], '.', ms=3, color='tab:purple', label='B3: PG2 Top')
        plt.plot(Bgroups[5][:, 0], Bgroups[5][:, 1], '.', ms=3, color='tab:pink',   label='B5: PG2 Bottom')
        plt.plot(Bgroups[7][:, 0], Bgroups[7][:, 1], '.', ms=3, color='tab:olive',  label='B7: PG1 Bottom')

        # Inflow and Outflow
        plt.plot(Bgroups[8][:, 0], Bgroups[8][:, 1], '.', ms=3, color='black',      label='B8: Inflow')
        plt.plot(Bgroups[4][:, 0], Bgroups[4][:, 1], '.', ms=3, color='tab:orange', label='B4: Outflow')

        # Blade Surfaces
        plt.plot(Bgroups[2][:, 0], Bgroups[2][:, 1], '.', ms=3, color='tab:red',    label='B2: Blade Lower')
        plt.plot(Bgroups[6][:, 0], Bgroups[6][:, 1], '.', ms=3, color='tab:blue',   label='B6: Blade Upper')

        # Formatting the plot for your report
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')
        plt.title('Discretized Domain Boundary Segments')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.axis('equal') 
        plt.grid(True, linestyle=':', alpha=0.4)
        plt.tight_layout()

        plt.savefig(figname, dpi=200, bbox_inches='tight', pad_inches=0.1)
    
    return Bnodes, Btmn_idx, PG, extent

if __name__ == "__main__":
    a = 17./13. # edge length of an equilateral triangle
    build_domain(a, "../../output/project1/domain_boundary.png")
    #build_domain()
