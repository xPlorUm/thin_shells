/*
Boilerplate code cleaned up from the sacred texts
Maybe we should change the function names since they're currently the same as the sacred texts...
*/

#include "DiscreteShell.h"
#include <iostream>
#include <igl/readOBJ.h>

// Constructor
DiscreteShell::DiscreteShell()
    : dt(0.01), simulation_duration(10.0), bending_stiffness(1.0),
      beta(0.25), gamma(0.5)
{
    // ChatGPT gave me these values idk
}

// Initialize from an OBJ file
// Might do away with this and just pass in the mesh from main
void DiscreteShell::initializeFromFile(const std::string& filename) {
    // Load mesh (vertices and faces)
    igl::readOBJ(filename, undeformed, faces); 
    deformed = undeformed;
    vn = Eigen::VectorXd::Zero(deformed.size()); // Initial velocity
    xn = undeformed; // Initial previous position
    u = Eigen::VectorXd::Zero(deformed.size()); // Initial displacement
    external_force = Eigen::VectorXd::Zero(deformed.size()); // No external force initially
}

// Advance one time step
bool DiscreteShell::advanceOneStep(int step) {
    std::cout << "Time step: " << step * dt << "s" << std::endl;

    Eigen::VectorXd residual = external_force;
    Eigen::SparseMatrix<double> K; // Stiffness matrix

    // Compute bending energy and forces
    double energy = 0.0;
    addShellBendingEnergy(energy);
    addShellBendingForce(residual);
    addShellBendingHessian(K);

    // Solve for displacement increment du
    Eigen::VectorXd du;
    bool success = linearSolve(K, residual, du);
    if (!success) {
        std::cout << "Linear solve failed." << std::endl;
        return false;
    }

    // Newmark Integration
    deformed = xn + dt * vn + beta * dt * dt * du; // Position update
    vn = vn + gamma * dt * du; // Velocity update

    updateDynamicStates();
    return (step * dt > simulation_duration); // End simulation after duration
}

// Compute total energy (only bending energy in this case)
double DiscreteShell::computeTotalEnergy() {
    double energy = 0.0;
    addShellBendingEnergy(energy);
    return energy;
}

// TODO Add bending energy
void DiscreteShell::addShellBendingEnergy(double& energy) {
    // Calculate bending energy based on deformed configuration
}

// TODO Add bending force to residual
void DiscreteShell::addShellBendingForce(Eigen::VectorXd& residual) {
    // Add bending forces to the residual vector
}

// TODO Add bending Hessian (stiffness matrix entries)
void DiscreteShell::addShellBendingHessian(Eigen::SparseMatrix<double>& K) {
    // Populate the stiffness matrix K with bending Hessian entries
}

// Linear solve for displacement increment (du)
// Here, you can experiment with different solvers by changing this function.
// Thanks ChatGPT for the placeholder code
bool DiscreteShell::linearSolve(Eigen::SparseMatrix<double>& K, const Eigen::VectorXd& residual, Eigen::VectorXd& du) {
    // TODO: Choose and set up a solver here, based on your needs
    // Example: Using Conjugate Gradient solver with placeholder code (replace with actual solver if desired)

    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    // solver.setTolerance(1e-8); // Adjust tolerance for convergence
    // solver.compute(K);
    // if (solver.info() != Eigen::Success) {
    //     std::cout << "Solver initialization failed." << std::endl;
    //     return false;
    // }
    // du = solver.solve(residual);
    // if (solver.info() != Eigen::Success) {
    //     std::cout << "Solver failed to converge." << std::endl;
    //     return false;
    // }
    
    // Placeholder return value; change this once the solver is implemented.
    return true;
}

// Build the system matrix (stiffness matrix with only bending contributions)
void DiscreteShell::buildSystemMatrix(Eigen::SparseMatrix<double>& K) {
    addShellBendingHessian(K);
}

// Update dynamic states (velocity and previous position) after each time step
void DiscreteShell::updateDynamicStates() {
    vn = (deformed - xn) / dt; // Update velocity
    xn = deformed; // Update previous position for the next time step
}
