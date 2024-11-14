/*
Boilerplate code cleaned up from the sacred texts
Maybe we should change the function names since they're currently the same as the sacred texts...
*/

#include "DiscreteShell.h"
#include <iostream>
#include <igl/readOBJ.h>
#include <igl/edges.h>
#include <autodiff/forward/dual/dual.hpp>

#include "Mesh.h"



Mesh undeformedMesh = Mesh(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXi::Zero(0, 3), Eigen::MatrixXi::Zero(0, 2));
Mesh deformedMesh = Mesh(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXi::Zero(0, 3), Eigen::MatrixXi::Zero(0, 2));

// Constructor
DiscreteShell::DiscreteShell()
    : dt(0.01), simulation_duration(10.0), bending_stiffness(1.0),
      beta(0.25), gamma(0.5)
{
}


// Helper Methods
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






// Initialize from an OBJ file
// Might do away with this and just pass in the mesh from main
void DiscreteShell::initializeFromFile(const std::string& filename) {
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
// Initialize from given Vertex
void DiscreteShell::initializeMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXd &F) {
    // Load mesh (vertices and faces)
    V_deformed = V;
    V_undeformed = V;

    vn = Eigen::VectorXd::Zero(V_deformed.size()); //initialize velocity vector
    xn = V_undeformed; // Initial previous position
    u = Eigen::VectorXd::Zero(V_deformed.size()); // Initial displacement
    external_force = Eigen::VectorXd::Zero(V_deformed.size()); // No external force initially
}


// Advance one time step
bool DiscreteShell::advanceOneStep(int step) {
    std::cout << "Time step: " << step * dt << "s" << std::endl;

    Eigen::VectorXd residual = external_force;
    Eigen::SparseMatrix<double> K; // Stiffness matrix

    // Compute bending energy and forces
    double energy = 0.0;
    energy = computeTotalEnergy();
    addShellBendingForce(residual);
    addShellBendingHessian(K);

    // Solve for displacement increment du
    Eigen::VectorXd du;
    bool success = linearSolve(K, residual, du);
    if (!success) {
        std::cout << "Linear solve failed." << std::endl;
        return false;
    }

    // given displacement increment du, dt, V_undeformed=xn
    // calculate acceleration
    // use auto differentiation?
    //Eigen::VectorXd acceleration = K.inverse() * du;

    // Newmark Integration
    //V_deformed = xn + dt * vn + dt * dt * ((0.5 - beta) * acceleration + beta * (K * du));

    // Velocity update
    //vn = vn + dt * ((1 - gamma) * acceleration + gamma * (K * du));


    //TODO return deformed/Vnew to main program

    updateDynamicStates();
    return (step * dt > simulation_duration); // End simulation after duration
}

// Compute total energy (only bending energy in this case)
double DiscreteShell::computeTotalEnergy() {
    double energy = 0.0f;
    energy += totalBendingEnergy();
    return energy;
}



double DiscreteShell::totalBendingEnergy() {
    // Calculate bending energy based on deformed configuration
    double energy = 0.0;
    // Loop over every edge
    for (int i = 0; i < undeformedMesh.edgeList.size(); ++i) {
        // Calculate bending energy for each edge
        energy += edgeBendingEnergy(i);
    }
    return energy;
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

     Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
     solver.setTolerance(Epsilon); // Adjust tolerance for convergence
     solver.compute(K);
     if (solver.info() != Eigen::Success) {
         std::cout << "Solver initialization failed." << std::endl;
         return false;
     }
     du = solver.solve(residual);
     if (solver.info() != Eigen::Success) {
         std::cout << "Solver failed to converge." << std::endl;
         return false;
     }
    
    // Placeholder return value; change this once the solver is implemented.
    return true;
}

// Build the system matrix (stiffness matrix with only bending contributions)
void DiscreteShell::buildSystemMatrix(Eigen::SparseMatrix<double>& K) {
    addShellBendingHessian(K);
}

// Update dynamic states (velocity and previous position) after each time step
void DiscreteShell::updateDynamicStates() {
    vn = (V_deformed - xn) / dt; // Update velocity
    xn = V_deformed; // Update previous position for the next time step
}

