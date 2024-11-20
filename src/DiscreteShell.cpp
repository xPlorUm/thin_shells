/*
Boilerplate code cleaned up from the sacred texts
Maybe we should change the function names since they're currently the same as the sacred texts...
*/

#include "DiscreteShell.h"
#include "igl/read_triangle_mesh.h"
#include <iostream>
#include <igl/readOBJ.h>
#include <filesystem>

// Constructor
DiscreteShell::DiscreteShell()
    : dt(0.01), simulation_duration(10.0), bending_stiffness(1.0),
      beta(0.25), gamma(0.5)
{
    F = new Eigen::MatrixXi(0, 3);
    V = new Eigen::MatrixXd(0, 3);
}

// Initialize from an OBJ file
// Might do away with this and just pass in the mesh from main
void DiscreteShell::initializeFromFile(const std::string& filename) {
    if (filename.empty()) {
        std::cerr << "Error: No file provided" << std::endl;
        return;
    }
    if (!std::filesystem::exists(filename)) {
        std::cerr << "Error: File does not exist" << std::endl;
        std::cerr << "You should provide a valid file as first argument !" << std::endl;
        return;
    }
    // check if filename is a real file
    // if not, return
    // Load mesh (vertices and faces)
    igl::read_triangle_mesh(filename, *V, *F);
    std::cout << "Loaded mesh " << filename << std::endl;
    std::cout << "V : (" << V->rows() << ", " << V->cols() << ")" << std::endl;
    std::cout << "F : (" << F->rows() << ", " << F->cols() << ")" << std::endl;
}

// Advance one time step
bool DiscreteShell::advanceOneStep(int step) {
    // std::cout << "Step " << step << std::endl;
    // Scale down V
    *V *= 0.9;
    return false;
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

const Eigen::MatrixXd* DiscreteShell::getPositions() {
    return V;
}

const Eigen::MatrixXi* DiscreteShell::getFaces() {
    return F;
}