#ifndef DISCRETE_SHELL_H
#define DISCRETE_SHELL_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "Mesh.h"
#include <iostream>
#include <igl/massmatrix.h>


class DiscreteShell {
public:
    // Constructor
    DiscreteShell();
    // Destructor
    ~DiscreteShell();

    void initializeFromFile(const std::string& filename);

    // Advance one time step
    bool advanceOneStep(int step);

    // Given an index of a vertex, print some stuff. Very much useful for debugging
    void debug_vertex(int vertex_index) {
        std::cout << "Vertex " << vertex_index << std::endl;
        // Velocity
        std::cout << "Velocity : " << Velocity->row(vertex_index) << std::endl;
        // Position
        std::cout << "Position : " << V->row(vertex_index) << std::endl;
        // Forces
        std::cout << "Forces : " << forces.row(vertex_index) << std::endl;
    }

// Returns the positions to get drawn
// TODO : this should be a pointer
const Eigen::MatrixXd* getPositions();
const Eigen::MatrixXi* getFaces();

private:

    // Physical properties
    double dt; // Time step
    double simulation_duration; // Total simulation duration
    // double bending_stiffness; // Stiffness for bending energy

    // State variables
    Mesh deformedMesh;
    Mesh undeformedMesh;

    Eigen::MatrixXi *F; // Faces of the shell
    Eigen::MatrixXd *V; // Current mesh
    Eigen::MatrixXd V_rest; // Undeformed mesh
    Eigen::MatrixXd *Velocity; // Velocity of the shell (point-wise)
    Eigen::MatrixXi *E; // Edges of the shell
    Eigen::VectorXd E_length_rest ; // Rest length of edges
    Eigen::SparseMatrix<double> M; // Massmatrix
    Eigen::MatrixX3d forces; // Forces applied point-wise.
    Eigen::MatrixX3d bending_forces; // Bending forces applied point-wise.
    int k_membrane = 1000; // Membrane stiffness
    double stiffness_damping = 0.5; // Stiffness damping coefficient
    double mass_damping = 0.5; // Mass damping coefficient

    // Time integration (Newmark scheme)
    void updateDynamicStates();
    bool linearSolve(Eigen::SparseMatrix<double>& K, const Eigen::VectorXd& residual, Eigen::VectorXd& du);

    // Newmark-specific parameters
    double beta = 0.25; // Newmark parameter (default: 0.25 for implicit integration)
    double gamma = 0.5; // Newmark parameter (default: 0.5 for implicit integration)

    // Helper function to build system matrix
    void buildSystemMatrix(Eigen::SparseMatrix<double>& K);

    void computeStrechingForces(Eigen::MatrixX3d& forces); // Compute stretching forces
    void computeBendingForces(Eigen::MatrixX3d& bending_forces); // Compute bending forces
    void computeDampingForces(Eigen::MatrixX3d& damping_forces); // Compute damping forces
    var totalBendingEnergy();
    var BendingEnergy(int i);


};

#endif // DISCRETE_SHELL_H
