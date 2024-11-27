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
#include <igl/decimate.h>
#include <fstream>
#include <iostream>

// Define a macro to enable/disable debug information
// #define DEBUG_VERTEX 0

int num_steps = 0;

// Constructor
DiscreteShell::DiscreteShell()
    : dt(0.03), simulation_duration(10.0),
    // bending_stiffness(1.0),
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

    Eigen::VectorXi J; // Output mapping from F to simplified faces

    // Specify the target number of faces
    int target_faces = 100;

    igl::decimate(*V, *F, target_faces, *V, *F, J);

    V_rest = Eigen::MatrixXd(*V);
    deformedMesh = Mesh(V_rest, *F);

    // Set only the y-coordinates to 0 for the undeformed mesh
    Eigen::MatrixXd V_undef = V_rest; // Copy vertices
    for (int i = 0; i < V_undef.rows(); ++i) {
        V_undef(i, 1) = 0.0; // Set the y-coordinate (index 1) to 0
    }

    undeformedMesh = Mesh(V_undef, *F);


    // Set massmatrix
    igl::massmatrix(V_rest, *F, igl::MASSMATRIX_TYPE_VORONOI, M);
    // to avoid instability
    M *= 500;

    // Get all edges
    std::cout << "Loaded mesh " << filename << std::endl;
    std::cout << "V : (" << V->rows() << ", " << V->cols() << ")" << std::endl;
    std::cout << "F : (" << F->rows() << ", " << F->cols() << ")" << std::endl;
    std::cout << "Faces : " << F->rows() << std::endl;

    // Compute the edges
    igl::edges(*F, *E);
    std::cout << "E : (" << E->rows() << ", " << E->cols() << ")" << std::endl;

    // Compute the rest length of the edges
    // This assumes that the undeformed mesh is the initial mesh
    // For paper, might be an okay assumption
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

// bool DiscreteShell::advanceOneStep(int step) {

//     // Reset forces
//     forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
//     computeStrechingForces(forces);

//     bending_forces = Eigen::MatrixX3d::Zero(V->rows(), 3);

//     //only make part time bending energy because it is so slow now
//     //if (step > 40) {
//     computeBendingForces(bending_forces);
//     // std::cout << "Bending force calculated" << std::endl;
//     forces += bending_forces;
//     //}

//     Eigen::MatrixX3d damping_forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
//     computeDampingForces(damping_forces);
//     forces += damping_forces;

//     // Open a file to log data (appending mode)
//     std::ofstream dataFile;
//     dataFile.open("vertex_data.csv", std::ios_base::app);

//     // Update positions with Newton's law
//     for (int i = 0; i < Velocity->rows(); i++) {
//         Eigen::Vector3d force = forces.row(i);
//         Eigen::Vector3d velocity = Velocity->row(i);

//         // If we are at vertex 0, log the velocity and force
//         // Log time, velocity, force, and position for vertex 0
//         if (i == 0) {
//             dataFile << num_steps++ << ","
//                     << V->row(0).x() << "," << V->row(0).y() << "," << V->row(0).z() << ","
//                     << velocity.x() << "," << velocity.y() << "," << velocity.z() << ","
//                     << force.x() << "," << force.y() << "," << force.z() << ","
//                     << bending_forces.row(0).x() << "," << bending_forces.row(0).y() << "," << bending_forces.row(0).z() << ","
//                     << damping_forces.row(0).x() << "," << damping_forces.row(0).y() << "," << damping_forces.row(0).z() << ","
//                     << "\n";
//         }

//         double mass_i = M.coeff(i, i);
//         Eigen::Vector3d acceleration = 1 / mass_i * force;
//         // actually scaling by the mass works better because it's less unstable but whatever this is "more correct"
//         Velocity->row(i) += dt * acceleration;
//         V->row(i) += dt * Velocity->row(i);
//     }

//     // Update deformed mesh
//     deformedMesh.V = *V;

//     // Close the file after logging
//     dataFile.close();

//     return false;
// }


bool DiscreteShell::advanceOneStep(int step) {
    // Reset forces
    forces.setZero(V->rows(), 3);

    // Compute forces
    computeStrechingForces(forces);
    computeBendingForces(bending_forces);
    forces += bending_forces;

    // Compute damping forces
    Eigen::MatrixX3d damping_forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
    computeDampingForces(damping_forces);
    forces += damping_forces;

    // Open a file to log data (appending mode)

    // Compute acceleration
    Eigen::MatrixX3d acceleration = Eigen::MatrixX3d::Zero(V->rows(), 3);
    for (int i = 0; i < Velocity->rows(); ++i) {
        double mass_i = M.coeff(i, i); // Mass for particle i
        acceleration.row(i) = forces.row(i) / mass_i;

        // If we are at vertex 0, log
        Eigen::Vector3d force = forces.row(i);
        Eigen::Vector3d velocity = Velocity->row(i);

#ifdef DEBUG_VERTEX
        // If we are at vertex 0, log the velocity and force
        // Log time, velocity, force, and position for vertex 0
        if (i == 0) {
            fileDebugger.getStream() << num_steps++ << ","
                    << V->row(0).x() << "," << V->row(0).y() << "," << V->row(0).z() << ","
                    << velocity.x() << "," << velocity.y() << "," << velocity.z() << ","
                    << force.x() << "," << force.y() << "," << force.z() << ","
                    << bending_forces.row(0).x() << "," << bending_forces.row(0).y() << "," << bending_forces.row(0).z() << ","
                    << damping_forces.row(0).x() << "," << damping_forces.row(0).y() << "," << damping_forces.row(0).z() << ","
                    << "\n";
        }
#endif

    }

    // Newmark integration loop
    for (int i = 0; i < V->rows(); ++i) {
        // Update position
        V->row(i) += dt * Velocity->row(i) +
                     dt * dt * ((1.0 - beta) * acceleration.row(i) + beta * acceleration.row(i));

        // Update velocity
        Velocity->row(i) += dt * ((1.0 - gamma) * acceleration.row(i) + gamma * acceleration.row(i));
    }

    // Update the deformed mesh
    deformedMesh.V = *V;

    return false; // Continue simulation
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
        double force_magnitude = k_membrane *  (current_length - rest_length);
        auto direction = (p2 - p1).normalized();
        Eigen::Vector3d fi = force_magnitude * direction;
        forces.row(v1) += fi;
        forces.row(v2) -= fi;
    }
}


void DiscreteShell::computeBendingForces(Eigen::MatrixX3d& bending_forces) {

    // std::cout << "Computing bending forces" << std::endl;
    
    bending_forces.setZero(V_rest.rows(), 3);

    for (int i = 0; i < V_rest.rows(); i++) {
        // Add bending forces
        auto energy_f = [this](int i) -> var {
            return BendingEnergy(i);
            };
        
        // Derivative and forward pass of Bending energy function
        // Shape of derivative is (#V, 3)
        var energy = energy_f(i); // forward
        DualVector x = deformedMesh.V.row(i);
        auto grad = gradient(energy, x);
        bending_forces.row(i) = -Eigen::Map<Eigen::Vector3d>(grad.data());
    }

}

var DiscreteShell::BendingEnergy(int i) {

    // std::cout << "Calculating bending energy for vertex " << i << std::endl;
    DualVector angles_d; 
    Eigen::VectorXd angles_u;
    // DualVector stiffness_u, stiffness_d;
    DualVector heights, norms;

    undeformedMesh.getDihedralAngles(i, angles_u); // don't need calculation for that

    deformedMesh.calculateDihedralAngles(i, angles_d);
    deformedMesh.computeAverageHeights(i, heights);
    deformedMesh.computeEdgeNorms(i, norms);
    //flexural energy per undirected edge
    DualVector flex = (((angles_u - angles_d).array().square() * norms.array()) / heights.array());
    // if (i == 0) std::cout << "Flexural energy for vertex " << i << " : " << flex.sum() << std::endl;
    return flex.sum();
}

void DiscreteShell::computeDampingForces(Eigen::MatrixX3d &damping_forces) {
    // Rayleigh type per-edge damping

    // Stiffness damping
    damping_forces.setZero();
    for (int i = 0; i < E->rows(); i++) {
        // Get the indices of the two vertices of the edge
        int v1 = E->operator()(i, 0);
        int v2 = E->operator()(i, 1);
        
        // Normalized direction of the edge
        Eigen::Vector3d edge = V->row(v2) - V->row(v1);
        Eigen::Vector3d edge_dir = edge.normalized();

        // Get the relaitve velocity of the two vertices
        Eigen::Vector3d vel1 = Velocity->row(v1);
        Eigen::Vector3d vel2 = Velocity->row(v2);

        double edge_stifness = deformedMesh.stiffness(i);

        // Compute the damping force
        Eigen::VectorXd damping_force = stiffness_damping * edge_stifness * (vel2 - vel1).dot(edge_dir) * edge_dir;

        // Add the damping force to the forces vector
        damping_forces.row(v1) += damping_force;
        damping_forces.row(v2) -= damping_force;
    }

    // Mass damping
    damping_forces += -1 * mass_damping * M * (*Velocity);

    // Damping force for vertex 0
    //td::cout << "Damping force for vertex 0 : " << damping_forces.row(0) << std::endl;
}

//void DiscreteShell::computeBendingForces(Eigen::MatrixX3d& bending_forces) {
//    bending_forces.setZero();
//
//    // Add bending forces
//    auto energy_f = [this]() -> var {
//        return totalBendingEnergy();
//        };
//
//    // Derivative and forward pass of Bending energy function
//    // Shape of derivative is (#V, 3)
//    var energy = energy_f(); // forward
//    DualVector x = Eigen::Map<DualVector>(deformedMesh.V.data(), deformedMesh.V.size(), 1);
//    bending_forces = -Eigen::Map<Eigen::MatrixXd>(
//        gradient(energy, x).data(), deformedMesh.V.rows(), deformedMesh.V.cols()
//    );
//    std::cout << x << std::endl;
//    std::cout << "x = " << x << std::endl;
//    std::cout << "Jacobian of energy with respect to x = " << bending_forces << std::endl;
//
//}

//var DiscreteShell::totalBendingEnergy() {
//    int nE = undeformedMesh.uE.rows();
//    assert(undeformedMesh.uE.rows() == deformedMesh.uE.rows());
//
//    DualVector angles_u(nE), angles_d(nE);
//    DualVector stiffness_u(nE), stiffness_d(nE);
//    DualVector heights(nE), norms(nE);
//
//
//    // Calculate the dihedral angles and the stiffness of the mesh
//    undeformedMesh.calculateDihedralAngle(angles_u, stiffness_u);
//    //std::cout << angles_u.sum() << std::endl;
//
//    deformedMesh.calculateDihedralAngle(angles_d, stiffness_d);
//    //deformedMesh.computeAverageHeights(heights);
//    //deformedMesh.computeEdgeNorms(norms);
//
//
//    DualVector flex(nE); //flexural energy per undirected edge
//    auto mask = (heights.array() != 0).cast<var>();
//    flex = mask * (((angles_u - angles_d).array().square() * norms.array()) / heights.array());
//
//    //std::cout << flex.sum() <<  std::endl;
//    //TODO multiply flex with Stiffness Matrix
//
//    return flex.sum();
//}


const Eigen::MatrixXd* DiscreteShell::getPositions() {
    return V;
}

const Eigen::MatrixXi* DiscreteShell::getFaces() {
    return F;
}


