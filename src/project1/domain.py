import numpy as np
import matplotlib.pyplot as plt

def build_domain(figname=""):

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

    # Pick every x-th point
    bladelower_coarse = bladelower[::4, :]
    bladeupper_coarse = bladeupper[::4, :]

    # x-nodes on boundary 1 and 2 match with each other
    boundary_1 = np.column_stack((PG1_up_xgrid, PG1_up_ygrid))[:-1, :]
    boundary_2 = bladelower_coarse[:-1, :]
    boundary_3 = np.column_stack((PG2_up_xgrid, PG2_up_ygrid))[:-1, :]
    boundary_4 = np.column_stack((outflow_xgrid, outflow_ygrid))[:-1, :]
    boundary_5 = np.column_stack((PG2_low_xgrid, PG2_low_ygrid))[:-1, :]
    boundary_6 = bladeupper_coarse[:-1, :]
    boundary_7 = np.column_stack((PG1_low_xgrid, PG1_low_ygrid))[:-1, :]
    boundary_8 = np.column_stack((inflow_xgrid, inflow_ygrid))[:-1, :]

    # combine all boundaries into a single array if needed
    boundaries = np.vstack((boundary_1, boundary_2, boundary_3, boundary_4,
                            boundary_5, boundary_6, boundary_7, boundary_8))

    # find the extents of the domain
    x_min, x_max = np.min(boundaries[:, 0]), np.max(boundaries[:, 0])
    y_min, y_max = np.min(boundaries[:, 1]), np.max(boundaries[:, 1])
    extent = [x_min, x_max, y_min, y_max]

    # Plotting the domain boundaries if a figure name is provided
    if (not not figname): 
        # Create the plot
        plt.figure(figsize=(12, 8))

        # Periodic Boundaries
        plt.plot(boundary_1[:, 0], boundary_1[:, 1], '.', ms=3, color='tab:green',  label='B1: PG1 Top')
        plt.plot(boundary_3[:, 0], boundary_3[:, 1], '.', ms=3, color='tab:purple', label='B3: PG2 Top')
        plt.plot(boundary_5[:, 0], boundary_5[:, 1], '.', ms=3, color='tab:pink',   label='B5: PG2 Bottom')
        plt.plot(boundary_7[:, 0], boundary_7[:, 1], '.', ms=3, color='tab:olive',  label='B7: PG1 Bottom')

        # Inflow and Outflow
        plt.plot(boundary_8[:, 0], boundary_8[:, 1], '.', ms=3, color='black',      label='B8: Inflow')
        plt.plot(boundary_4[:, 0], boundary_4[:, 1], '.', ms=3, color='tab:orange', label='B4: Outflow')

        # Blade Surfaces
        plt.plot(boundary_2[:, 0], boundary_2[:, 1], '.', ms=3, color='tab:red',    label='B2: Blade Lower')
        plt.plot(boundary_6[:, 0], boundary_6[:, 1], '.', ms=3, color='tab:blue',   label='B6: Blade Upper')

        # Formatting the plot for your report
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')
        plt.title('Discretized Domain Boundary Segments')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.axis('equal') 
        plt.grid(True, linestyle=':', alpha=0.4)
        plt.tight_layout()

        plt.savefig(figname, dpi=200, bbox_inches='tight', pad_inches=0.1)
    
    return boundaries, extent

if __name__ == "__main__":
    build_domain("../../output/project1/domain_boundary.png")
    #build_domain()
