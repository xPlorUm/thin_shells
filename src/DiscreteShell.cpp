/*
Boilerplate code cleaned up from the sacred texts
Maybe we should change the function names since they're currently the same as the sacred texts...
*/

#include "DiscreteShell.h"
#include "igl/read_triangle_mesh.h"
#include "common.h"
#include <iostream>
#include <igl/edges.h> // Include this header to use igl::edges
#include <filesystem>


// Constructor
DiscreteShell::DiscreteShell()
    : dt(0.03), simulation_duration(10.0), bending_stiffness(1.0),
      beta(0.25), gamma(0.5)
{
    F = new Eigen::MatrixXi(0, 3);
    V = new Eigen::MatrixXd(0, 3);
    Velocity = new Eigen::MatrixXd(0, 3);
    E = new Eigen::MatrixXi(0, 2);
    forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
}

// Destructor
DiscreteShell::~DiscreteShell() {
    delete F;
    delete V;
    delete Velocity;
    delete E;
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

    // Load mesh (vertices and faces)
    igl::read_triangle_mesh(filename, *V, *F);
    V_rest = *V;
    undeformedMesh = Mesh(V_rest, *F);
    deformedMesh = undeformedMesh;

    // Set massmatrix
    igl::massmatrix(V_rest, *F, igl::MASSMATRIX_TYPE_VORONOI, M);


    // Get all edges
    std::cout << "Loaded mesh " << filename << std::endl;
    std::cout << "V : (" << V->rows() << ", " << V->cols() << ")" << std::endl;
    std::cout << "F : (" << F->rows() << ", " << F->cols() << ")" << std::endl;
    std::cout << "Faces : " << F->rows() << std::endl;

    // Compute the edges
    igl::edges(*F, *E);
    std::cout << "E : (" << E->rows() << ", " << E->cols() << ")" << std::endl;

    // Compute the rest length of the edges
    E_length_rest.resize(E->rows());
    for (int i = 0; i < E->rows(); i++) {
        // Get the indices of the two vertices of the edge
        int v1 = (*E)(i, 0);
        int v2 = (*E)(i, 1);
        // Get the positions of the two vertices
        Eigen::Vector3d p1 = V->row(v1);
        Eigen::Vector3d p2 = V->row(v2);
        // Compute the rest length of the edge
        E_length_rest(i) = (p2 - p1).norm();
    }

    // Set velocity to zero everywhere
    Velocity->resize(V->rows(), 3);
    Velocity->setZero();

}

// Advance one time step
bool DiscreteShell::advanceOneStep(int step) {
    // Reset forces
    forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
    computeStrechingForces(forces);

    bending_forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
    computeBendingForces(bending_forces);
    forces += bending_forces;


    // Update positions with Newtons law
    // Iterate over all velocities and update
    // Set first vertice force as -inf
    for (int i = 0; i < Velocity->rows(); i++) {
        // Gets the force for the vertex I
        Eigen::Vector3d force = forces.row(i);
        // Gets the velocity for the vertex I
        double mass_i = 1;
        if (i == 42) {
            mass_i = 0;
        }
        // Adds gravity as force
        force += Eigen::Vector3d(0, -9.81, 0);
        Eigen::Vector3d acceleration = mass_i * force;
        // add gravity
        Velocity->row(i) += dt * acceleration;
        // Update the position
        V->row(i) += dt * Velocity->row(i);
    }

    // Update deformed mesh
    deformedMesh.V = *V;

    return false;
}


void DiscreteShell::computeStrechingForces(Eigen::MatrixX3d &forces) {
    // TODO : that should use baraff triangle methods.
    // Compute stretching forces
    // Iterate over all edges
    // Compute the stretching force for each edge
    // Add the stretching force to the forces vector
    forces.setZero();
    for (int i = 0; i < E->rows(); i++) {
        // Get the indices of the two vertices of the edge
        int v1 = E->operator()(i, 0);
        int v2 = E->operator()(i, 1);
        // Get the positions of the two vertices
        Eigen::Vector3d p1 = V->row(v1);
        Eigen::Vector3d p2 = V->row(v2);
        // Compute the stretching force
        double current_length = (p2 - p1).norm();
        double rest_length = E_length_rest(i);
        double force_magnitude = 0.5 * 1000 *  (current_length - rest_length);
        auto direction = (p2 - p1).normalized();
        Eigen::Vector3d fi = force_magnitude * direction;
        forces.row(v1) += fi;
        forces.row(v2) -= fi;
    }
}

void DiscreteShell::computeBendingForces(Eigen::MatrixX3d& bending_forces) {
    bending_forces.setZero();
    // Add bending forces
    auto energy_f = [this]() -> var {
        return totalBendingEnergy();
        };

    // Derivative and forward pass of Bending energy function
    // Shape of derivative is (#V, 3)
    var energy = energy_f(); // forward
    DualVector x = Eigen::Map<DualVector>(deformedMesh.V.data(), deformedMesh.V.size(), 1);
    bending_forces = -Eigen::Map<Eigen::MatrixXd>(
        gradient(energy, x).data(), deformedMesh.V.rows(), deformedMesh.V.cols()
    );
    std::cout << x << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "Jacobian of energy with respect to x = " << bending_forces << std::endl;

}

var DiscreteShell::totalBendingEnergy() {
    int nE = undeformedMesh.uE.rows();
    assert(undeformedMesh.uE.rows() == deformedMesh.uE.rows());

    DualVector angles_u(nE), angles_d(nE);
    DualVector stiffness_u(nE), stiffness_d(nE);
    DualVector heights(nE), norms(nE);


    // Calculate the dihedral angles and the stiffness of the mesh
    undeformedMesh.calculateDihedralAngle(angles_u, stiffness_u);
    //std::cout << angles_u.sum() << std::endl;

    deformedMesh.calculateDihedralAngle(angles_d, stiffness_d);
    deformedMesh.computeAverageHeights(heights);
    //deformedMesh.computeEdgeNorms(norms);


    DualVector flex(nE); //flexural energy per undirected edge
    auto mask = (heights.array() != 0).cast<var>();
    flex = mask * (((angles_u - angles_d).array().square() * norms.array()) / heights.array());

    //std::cout << flex.sum() <<  std::endl;
    //TODO multiply flex with Stiffness Matrix

    return flex.sum();
}


const Eigen::MatrixXd* DiscreteShell::getPositions() {
    return V;
}

const Eigen::MatrixXi* DiscreteShell::getFaces() {
    return F;
}


// Linear solve for displacement increment (du)
// Here, you can experiment with different solvers by changing this function.
// Thanks ChatGPT for the placeholder code
//bool DiscreteShell::linearSolve(Eigen::SparseMatrix<double>& K, const Eigen::VectorXd& residual, Eigen::VectorXd& du) {
//
//     Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
//     solver.setTolerance(Epsilon); // Adjust tolerance for convergence
//     solver.compute(K);
//     if (solver.info() != Eigen::Success) {
//         std::cout << "Solver initialization failed." << std::endl;
//         return false;
//     }
//     du = solver.solve(residual);
//     if (solver.info() != Eigen::Success) {
//         std::cout << "Solver failed to converge." << std::endl;
//         return false;
//     }
//    
//    // Placeholder return value; change this once the solver is implemented.
//    return true;
//}
//
//// Build the system matrix (stiffness matrix with only bending contributions)
//void DiscreteShell::buildSystemMatrix(Eigen::SparseMatrix<double>& K) {
//    addShellBendingHessian(K);
//}
//
//// Update dynamic states (velocity and previous position) after each time step
//void DiscreteShell::updateDynamicStates() {
//    //vn = (V_deformed - xn) / dt; // Update velocity
//    //xn = V_deformed; // Update previous position for the next time step
//}
