//
// Created by hugues on 12/12/24.
//

#ifndef THINSHELLS_CONSTANTS_H
#define THINSHELLS_CONSTANTS_H

#define BENDING_STIFFNESS  20.0
#define STRETCHING_STIFFNESS 10
#define AREA_PRESERVATION_STIFFNESS 10
#define DAMPING = 0.1

// Define the number of faces igl will decimate the mesh into.
// #define DECIMATE_N_FACES_TARGET_MESH 400
// #define UPSAMPLE_FACTOR 4

#define TOLERANCE_NEWTON 1e-3
#define MAX_ITERATIONS_NEWTON 400

#define SIMULATION_DT 1
#define SIMULATION_NEWMARK_BETA 0.25
#define SIMULATION_NEWMARK_GAMMA 0.5

// Enable or disable stretching forces in the newmark algorithm. Useful to debug some stuff.
#define ENABLE_STRETCHING_FORCES true
#define ENABLE_AREA_PRESERVATION_FORCES true
#define ENABLE_BENDING_FORCES true

#define APPLY_INITIAL_CHANGES(V) \

#define SAVE_FRAMES 1

inline void add_F_ext_internal(Eigen::MatrixXd &_forces, int step) {
    // Just gravity for now
    // Substract gravity :
    Eigen::RowVector3d gravity(0, -9.81, 0);
    _forces.rowwise() += gravity / 6;
    // Wind only on vertex 42
    Eigen::RowVector3d wind(0, 50, 0);
    _forces.row(0) += wind;
}


#endif //THINSHELLS_CONSTANTS_H
