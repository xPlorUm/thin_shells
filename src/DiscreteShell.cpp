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



#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>



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

    // Specify the target number of faces
    int target_faces = 100;

    //Eigen::VectorXi J; // Output mapping from F to simplified faces

    //igl::decimate(*V, *F, target_faces, *V, *F, J);

    V_rest = Eigen::MatrixXd(*V);
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

    //only make part time bending energy because it is so slow now
    computeBendingForces(bending_forces);
    forces += bending_forces;

    //add gravity
    forces.col(1) += Eigen::VectorXd::Constant(forces.rows(), - 9.81);

    Eigen::VectorXd mass = M.diagonal();
    mass = mass.cwiseInverse();
    //Eigen::VectorXd mass = Eigen::VectorXd::Constant(forces.rows(), 1.0);
    mass(42) = 0.0f; //handle

    Eigen::MatrixXd acceleration(forces.rows(), forces.cols());
    acceleration = forces.array().colwise() * mass.array();
    *Velocity = Velocity->array() + dt * acceleration.array();
    *V = V->array() + dt * Velocity->array();


    // Update deformed mesh
    deformedMesh.V = *V;

    return false;
}

//bool DiscreteShell::advanceOneStep(int step) {
//    // Reset forces
//    forces.setZero(V->rows(), 3);
//
//    // Compute forces
//    computeStrechingForces(forces);
//    computeBendingForces(bending_forces);
//    forces += bending_forces;
//
//    // Compute damping forces
//    Eigen::MatrixX3d damping_forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
//    computeDampingForces(damping_forces);
//    forces += damping_forces;
//
//    // Open a file to log data (appending mode)
//
//    // Compute acceleration
//    Eigen::MatrixX3d acceleration = Eigen::MatrixX3d::Zero(V->rows(), 3);
//    for (int i = 0; i < Velocity->rows(); ++i) {
//        double mass_i = M.coeff(i, i); // Mass for particle i
//        acceleration.row(i) = forces.row(i) / mass_i;
//
//        // If we are at vertex 0, log
//        Eigen::Vector3d force = forces.row(i);
//        Eigen::Vector3d velocity = Velocity->row(i);
//
//#ifdef DEBUG_VERTEX
//        // If we are at vertex 0, log the velocity and force
//        // Log time, velocity, force, and position for vertex 0
//        if (i == 0) {
//            fileDebugger.getStream() << num_steps++ << ","
//                    << V->row(0).x() << "," << V->row(0).y() << "," << V->row(0).z() << ","
//                    << velocity.x() << "," << velocity.y() << "," << velocity.z() << ","
//                    << force.x() << "," << force.y() << "," << force.z() << ","
//                    << bending_forces.row(0).x() << "," << bending_forces.row(0).y() << "," << bending_forces.row(0).z() << ","
//                    << damping_forces.row(0).x() << "," << damping_forces.row(0).y() << "," << damping_forces.row(0).z() << ","
//                    << "\n";
//        }
//#endif
//
//    }
//
//    // Newmark integration loop
//    for (int i = 0; i < V->rows(); ++i) {
//        // Update position
//        V->row(i) += dt * Velocity->row(i) +
//                     dt * dt * ((1.0 - beta) * acceleration.row(i) + beta * acceleration.row(i));
//
//        // Update velocity
//        Velocity->row(i) += dt * ((1.0 - gamma) * acceleration.row(i) + gamma * acceleration.row(i));
//    }
//
//    // Update the deformed mesh
//    deformedMesh.V = *V;
//
//    return false; // Continue simulation
//}

void DiscreteShell::computeStrechingForces(Eigen::MatrixX3d& forces) {
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
    bending_forces.setZero(V_rest.rows(), 3);
    //Eigen::VectorXd stiffness_vec = Eigen::VectorXd::Zero(V_rest.rows());

    //3 variables per vertex. Objective per edge, using 3 vertices each :
    auto func = TinyAD::scalar_function<3>(TinyAD::range(deformedMesh.V.rows()));

    func.add_elements<2>(TinyAD::range(deformedMesh.uE.rows()), [&](auto& element)->TINYAD_SCALAR_TYPE(element) {
        using T = TINYAD_SCALAR_TYPE(element);

        // Variables associated with the edge's vertices
        Eigen::Index edge_idx = element.handle;
        Eigen::RowVector3<T> v0 = element.variables(deformedMesh.uE(edge_idx, 0));
        Eigen::RowVector3<T> v1 = element.variables(deformedMesh.uE(edge_idx, 1));


        // Compute dihedral angle, height, and norm
        double angle_u = undeformedMesh.dihedralAngles(edge_idx);

        T angle_d = deformedMesh.computeDihedralAngle(v0, v1, edge_idx);
        if (angle_d > 2)
            std::cout << angle_d << std::endl;
        angle_d = deformedMesh.computeDihedralAngle(v0, v1, edge_idx);

        deformedMesh.dihedralAngles(edge_idx) = TinyAD::to_passive(angle_d);
        // TODO: change the resting angle of the undeformed mesh

        T height = deformedMesh.computeHeight(v0, v1, edge_idx);
        T norm = (v1 - v0).norm();
        double stiffness = deformedMesh.stiffness(edge_idx);

        // Edge Case
        if (height <= Epsilon || norm <= Epsilon)
            return (T)0.0;

        // Compute bending energy for this edge
        T bending_energy = ((angle_u - angle_d) * (angle_u - angle_d) * norm * stiffness) / height;

        return bending_energy;
    });

    
    Eigen::VectorXd x = func.x_from_data([&](int vertex_idx) -> Eigen::VectorXd{
        return deformedMesh.V.row(vertex_idx);
        });

    auto [f, g] = func.eval_with_gradient(x);

    for (int vertex_idx = 0; vertex_idx < g.size() / 3; vertex_idx++) {
        bending_forces.row(vertex_idx) += g.segment<3>(3 * vertex_idx);
    }
    

}

const Eigen::MatrixXd* DiscreteShell::getPositions() {
    return V;
}

const Eigen::MatrixXi* DiscreteShell::getFaces() {
    return F;
}

const Eigen::MatrixXd DiscreteShell::getNormals() {
    return deformedMesh.FN;
}



