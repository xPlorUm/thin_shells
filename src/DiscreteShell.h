#ifndef DISCRETE_SHELL_H
#define DISCRETE_SHELL_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>

class DiscreteShell {
public:
    // Constructor
    DiscreteShell();
    // Destructor
    ~DiscreteShell();

    // Initialize from an OBJ file
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
    double bending_stiffness; // Stiffness for bending energy

    // State variables
    Eigen::VectorXd deformed; // Deformed configuration (vertex positions)
    Eigen::VectorXd undeformed; // Undeformed configuration (reference positions)
    Eigen::VectorXd external_force; // External forces applied to the shell
    Eigen::VectorXd u; // Displacement vector
    Eigen::VectorXd vn; // Velocity vector
    Eigen::VectorXd xn; // Previous position vector (Newmark integration)

    Eigen::MatrixXi *F; // Faces of the shell
    Eigen::MatrixXd *V;
    Eigen::MatrixXd *Velocity; // Velocity of the shell (point-wise)
    Eigen::MatrixXi *E; // Edges of the shell
    Eigen::VectorXd E_length_rest ; // Rest length of edges
    // Mass matrix.
    Eigen::MatrixXd M_inv;
    // Forces applied point-wise.
    Eigen::MatrixX3d forces;

    // Energy and force computation
    double computeTotalEnergy();
    void addShellBendingEnergy(double& energy);
    void addShellBendingForce(Eigen::VectorXd& residual);
    void addShellBendingHessian(Eigen::SparseMatrix<double>& K);

    // Time integration (Newmark scheme)
    void updateDynamicStates();
    bool linearSolve(Eigen::SparseMatrix<double>& K, const Eigen::VectorXd& residual, Eigen::VectorXd& du);

    // Newmark-specific parameters
    double beta; // Newmark parameter (default: 0.25 for implicit integration)
    double gamma; // Newmark parameter (default: 0.5 for implicit integration)

    // Helper function to build system matrix
    void buildSystemMatrix(Eigen::SparseMatrix<double>& K);

    void computeStrechingForces(Eigen::MatrixX3d &forces); // Compute stretching forces


};

#endif // DISCRETE_SHELL_H
