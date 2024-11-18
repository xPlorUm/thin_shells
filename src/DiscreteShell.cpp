/*
Boilerplate code cleaned up from the sacred texts
Maybe we should change the function names since they're currently the same as the sacred texts...
*/

#include "DiscreteShell.h"
#include <iostream>
#include <igl/readOBJ.h>
#include <igl/edges.h>




// Constructor
DiscreteShell::DiscreteShell()
    : dt(0.01), simulation_duration(10.0), bending_stiffness(1.0),
      beta(0.25), gamma(0.5)
{
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


    // Initialize the undeformed mesh
    undeformedMesh = Mesh(V, F);
    deformedMesh = Mesh(V, F);

    vn = Eigen::MatrixXd::Zero(V.rows(), 3); // Initialize velocity
    xn = V;
    u = Eigen::MatrixXd::Zero(V.rows(), 3); // Initialize displacement
    external_force = Eigen::MatrixXd::Zero(V.rows(), 3); // No external force initially
}

// Initialize from given Vertex
void DiscreteShell::initializeMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    undeformedMesh = Mesh(V, F);
    deformedMesh = Mesh(V, F);

    vn = Eigen::MatrixXd::Zero(V.rows(), 3); // Initialize velocity
    xn = V;
    u = Eigen::MatrixXd::Zero(V.rows(), 3); // Initialize displacement
    external_force = Eigen::MatrixXd::Zero(V.rows(), 3); // No external force initially
}


// Advance one time step
bool DiscreteShell::advanceOneStep(int step) {
    std::cout << "Time step: " << step * dt << "s" << std::endl;

    Eigen::VectorXd residual = external_force;

    int nV = deformedMesh.V.rows();
    // Compute bending energy and forces
    //dual energy = totalBendingEnergy();

    //TODO Compute derivative of energy
    auto energy_func = [this]() {
        return totalBendingEnergy();
        };
    DualVector gx = gradient(energy_func, wrt(deformedMesh.V), at(deformedMesh.V));




    //addShellBendingForce(residual);
    //addShellBendingHessian(K);

    // Solve for displacement increment du
    //Eigen::VectorXd du;
    //bool success = linearSolve(K, residual, du);
    //if (!success) {
    //    std::cout << "Linear solve failed." << std::endl;
    //    return false;
    //}

    //solve ODE of motion





    // given displacement increment du, dt, V_undeformed=xn
    // calculate acceleration

    //Eigen::VectorXd acceleration = K.inverse() * du;

    // Newmark Integration
    //V_deformed = xn + dt * vn + dt * dt * ((0.5 - beta) * acceleration + beta * (K * du));

    // Velocity update
    //vn = vn + dt * ((1 - gamma) * acceleration + gamma * (K * du));


    //TODO return deformed/Vnew to main program

    updateDynamicStates();
    return (step * dt > simulation_duration); // End simulation after duration
}



dual DiscreteShell::totalBendingEnergy() {
    int nE = undeformedMesh.uE.rows();
    assert(undeformedMesh.uE.rows() == deformedMesh.uE.rows());

    DualVector angles_u(nE), angles_d(nE);
    DualVector stiffness_u(nE), stiffness_d(nE);
    DualVector heights(nE), norms(nE);
    DualVector flex(nE); //flexural energy per undirected edge


    // Calculate the dihedral angles and the stiffness of the mesh
    undeformedMesh.calculateDihedralAngle(angles_u, stiffness_u);
    deformedMesh.calculateDihedralAngle(angles_d, stiffness_d);


    deformedMesh.computeAverageHeights(heights);

    deformedMesh.computeEdgeNorms(norms);


    flex = ((angles_u - angles_d).array().square() * norms.array()) / heights.array();

    // Handle case when heights entries are 0
    flex = flex.unaryExpr([](dual val) { return std::isfinite(val.val) ? val : dual(0.0); });


    return flex.sum();
}





//dual DiscreteShell::edgeBendingEnergy(int edgeIndex) {
//    // Retrieve the edge from both undeformed and deformed mesh configurations
//    Edge& edge_undeformed = undeformedMesh.edgeList[edgeIndex];
//    Edge& edge_deformed = deformedMesh.edgeList[edgeIndex];
//
//    // Calculate the length of the edge in the undeformed configuration
//    Eigen::Matrix<dual, 3, 1> v1 = undeformedMesh.vertices.row(edge_undeformed.v1).cast<dual>();
//    Eigen::Matrix<dual, 3, 1> v2 = undeformedMesh.vertices.row(edge_undeformed.v2).cast<dual>();
//    dual edgeLength = (v2 - v1).norm();
//
//    // Retrieve opposite vertices in undeformed and deformed configurations
//    int oppositeVertex1 = findOppositeVertex(undeformedMesh.faces, edge_undeformed.v1, edge_undeformed.v2);
//    int oppositeVertex2 = findOppositeVertex(deformedMesh.faces, edge_deformed.v1, edge_deformed.v2);
//
//    // Calculate heights in undeformed and deformed configurations
//    dual h1 = computeHeight(v1, v2, undeformedMesh.vertices.row(oppositeVertex1).cast<dual>());
//    dual h2 = computeHeight(v1, v2, deformedMesh.vertices.row(oppositeVertex2).cast<dual>());
//    dual he = (h1 + h2) / 6.0;
//
//    // Calculate the difference in dihedral angles between deformed and undeformed configurations
//    dual angleDifference = edge_deformed.dihedralAngle - edge_undeformed.dihedralAngle.cast<dual>();
//
//    // Compute the bending energy contribution of this edge
//    dual bendingEnergy = pow(angleDifference, 2) * (edgeLength / he);
//
//    return bendingEnergy;
//}



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
    //vn = (V_deformed - xn) / dt; // Update velocity
    //xn = V_deformed; // Update previous position for the next time step
}


