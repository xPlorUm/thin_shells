/*
Boilerplate code cleaned up from the sacred texts
Maybe we should change the function names since they're currently the same as the sacred texts...
*/

#include "DiscreteShell.h"
#include <iostream>
#include <igl/readOBJ.h>

Mesh undeformedMesh;
Mesh deformedMesh;

// Constructor
DiscreteShell::DiscreteShell()
    : dt(0.01), simulation_duration(10.0), bending_stiffness(1.0),
      beta(0.25), gamma(0.5)
{
    // ChatGPT gave me these values idk
}

// Initialize from an OBJ file
// Might do away with this and just pass in the mesh from main
void initializeFromFile(const std::string& filename) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXi E;
    if (!igl::readOBJ(filename, V, F)) {
        std::cerr << "Failed to read OBJ file." << std::endl;
        return;
    }

    igl::edges(F, E);

    // Initialize the undeformed mesh
    undeformedMesh = Mesh(V, F, E);
    deformedMesh = Mesh(V, F, E);

    vn = Eigen::MatrixXd::Zero(V.rows(), 3); // Initialize velocity
    u = Eigen::MatrixXd::Zero(V.rows(), 3); // Initialize displacement
    external_force = Eigen::MatrixXd::Zero(V.rows(), 3); // No external force initially
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
    totalBendingEnergy();
    return energy;
}



// TODO Add bending energy
void DiscreteShell::totalBendingEnergy() {
    // Calculate bending energy based on deformed configuration
    double energy = 0.0;

    // Loop over every edge
    for (int i = 0; i < undeformedMesh.edgeList.size(); ++i) {
        // Calculate bending energy for each edge
        energy += edgeBendingEnergy(i);
    }
}

int findOppositeVertex(const Eigen::MatrixXi& triangle, int v1, int v2) {
    for (int i = 0; i < 3; ++i) {
        int vertex = triangle(0, i);
        if (vertex != v1 && vertex != v2) {
            return vertex;
        }
    }
    return -1; // Should never happen if the mesh is correctly formed
}

double computeHeight(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
    Eigen::Vector3d edge = v1 - v0;  // Edge vector
    Eigen::Vector3d vertex = v2 - v0;  // Vertex vector

    // Projection of vertex onto edge
    Eigen::Vector3d proj = (vertex.dot(edge) / edge.squaredNorm()) * edge;

    // Perpendicular component from edge to vertex
    Eigen::Vector3d perp = vertex - proj;

    // Height is the norm of the perpendicular component
    return perp.norm();
}


double DiscreteShell::edgeBendingEnergy(int edgeIndex) {
    // Retrieve the edge from both undeformed and deformed mesh configurations
    Edge& edge_undeformed = undeformedMesh.edgeList[edgeIndex];
    Edge& edge_deformed = deformedMesh.edgeList[edgeIndex];

    // Calculate the length of the edge in the undeformed configuration
    Eigen::Vector3d v1 = undeformedMesh.vertices.row(edge_undeformed.v1);
    Eigen::Vector3d v2 = undeformedMesh.vertices.row(edge_undeformed.v2);
    double edgeLength = (v2 - v1).norm();

    // Assuming both triangles are available and correctly referenced
    int oppositeVertex1 = findOppositeVertex(undeformedMesh.faces, edge_undeformed.v1, edge_undeformed.v2);
    int oppositeVertex2 = findOppositeVertex(deformedMesh.faces, edge_deformed.v1, edge_deformed.v2);

    double h1 = computeHeight(v1, v2, undeformedMesh.vertices.row(oppositeVertex1));
    double h2 = computeHeight(v1, v2, deformedMesh.vertices.row(oppositeVertex2));
    double he = (h1 + h2) / 6.0;

    // Calculate the difference in dihedral angles between deformed and undeformed configurations
    double angleDifference = edge_deformed.dihedralAngle - edge_undeformed.dihedralAngle;

    // Compute the bending energy contribution of this edge
    double bendingEnergy = pow(angleDifference, 2) * (edgeLength / he);

    return bendingEnergy;
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

