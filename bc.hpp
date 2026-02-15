#ifndef BC_HPP
#define BC_HPP

#include "readgri.hpp"

/** Boundary condition types. BGroup2,6->WALL; BGroup8->INFLOW; BGroup4->OUTFLOW. */
enum BCType {
    BC_WALL,
    BC_INFLOW,
    BC_OUTFLOW,
    BC_UNKNOWN
};

/** Get BC type for boundary face i from BGroup name. */
BCType get_bc_type(const GriMesh& mesh, int bface_idx);

#endif
