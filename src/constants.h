//
// Created by hugues on 12/12/24.
//

#ifndef THINSHELLS_CONSTANTS_H
#define THINSHELLS_CONSTANTS_H

#define BENDING_STIFFNESS  10.0
#define STRETCHING_STIFFNESS 1
#define AREA_PRESERVATION_STIFFNESS 10
#define DAMPING = 0.1

// Define the number of faces igl will decimate the mesh into.
#define N_FACES_MESH 400

#define TOLERANCE_NEWTON 1e-3
#define MAX_ITERATIONS_NEWTON 100

#define SIMULATION_DT 1
#define SIMULATION_NEWMARK_BETA 0.25
#define SIMULATION_NEWMARK_GAMMA 0.5

// Enable or disable stretching forces in the newmark algorithm. Useful to debug some stuff.
#define ENABLE_STRETCHING_FORCES true
#define ENABLE_AREA_PRESERVATION_FORCES true


#define APPLY_INITIAL_CHANGES(V) \
/*
    (V)->row(0)      += Eigen::Vector3d(0.0, 0.0, 30.0); \
    (V)->row(1)      += Eigen::Vector3d(0.0, 0.0, 30.0); \
    (V)->row(2)      += Eigen::Vector3d(0.0, 0.0, 30.0);
    (V)->row(14).z() += 30; \
    (V)->row(42)     += Eigen::Vector3d(0.0, 0.0, 30.0); \
    (V)->row(43)     += Eigen::Vector3d(0.0, 0.0, 30.0); \
*/


#endif //THINSHELLS_CONSTANTS_H
