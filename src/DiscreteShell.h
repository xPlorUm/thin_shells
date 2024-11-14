#ifndef DISCRETE_SHELL_H
#define DISCRETE_SHELL_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class DiscreteShell {
public:
    // Constructor
    DiscreteShell::DiscreteShell();
    void initializeMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXd& F);

    // Advance one time step
    bool advanceOneStep(int step);

private:
    // Numerical stability
    double Epsilon = 1e-4;

    // Physical properties
    double dt; // Time step
    double simulation_duration; // Total simulation duration
    double bending_stiffness; // Stiffness for bending energy

    // State variables
    Eigen::MatrixXd V_deformed; // Deformed configuration (vertex positions)
    Eigen::MatrixXd V_undeformed; // Undeformed configuration (reference positions)
    Eigen::VectorXd external_force; // External forces applied to the shell
    Eigen::VectorXd u; // Displacement vector
    Eigen::VectorXd vn; // Velocity vector
    Eigen::VectorXd xn; // Previous position vector (Newmark integration)


    // Energy and force computation
    double computeTotalEnergy();
    void addShellBendingForce(Eigen::VectorXd& residual);
    void addShellBendingHessian(Eigen::SparseMatrix<double>& K);

    double totalBendingEnergy();
    double edgeBendingEnergy(int edgeIndex);

    // Time integration (Newmark scheme)
    void updateDynamicStates();
    bool linearSolve(Eigen::SparseMatrix<double>& K, const Eigen::VectorXd& residual, Eigen::VectorXd& du);

    // Newmark-specific parameters
    double beta; // Newmark parameter (default: 0.25 for implicit integration)
    double gamma; // Newmark parameter (default: 0.5 for implicit integration)

    // Helper function to build system matrix
    void buildSystemMatrix(Eigen::SparseMatrix<double>& K);

};

#endif // DISCRETE_SHELL_H
