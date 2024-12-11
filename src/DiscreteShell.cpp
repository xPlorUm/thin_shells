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
    : m_dt(0.03), simulation_duration(10.0),
    // bending_stiffness(1.0),
      m_beta(0.25), m_gamma(0.5)
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

    // Itialize the solver
    m_solver = new Solver(m_beta, m_gamma, m_dt, &M);
    m_solver->setDiscreteShell(this);

}


bool DiscreteShell::advanceOneStep(int step) {
    // Reset forces
    forces.setZero(V->rows(), 3);
    // Compute forces
    addStrechingForcesTo(forces, V);
    // addBendingForcesTo(forces, V, E);

    // Compute damping forces
    Eigen::MatrixX3d damping_forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
    computeDampingForces(damping_forces);
    forces += damping_forces;

    // Open a file to log data (appending mode)

    // Compute acceleration
    Eigen::MatrixXd acceleration = Eigen::MatrixXd::Zero(V->rows(), 3);
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
    m_solver->solve(V, Velocity, &acceleration);

    // Newmark integration loop
    for (int i = 0; i < V->rows(); ++i) {
        // Update position
        V->row(i) += m_dt * Velocity->row(i) +
                     m_dt * m_dt * ((1.0 - m_beta) * acceleration.row(i) + m_beta * acceleration.row(i));

        // Update velocity
        Velocity->row(i) += m_dt * ((1.0 - m_gamma) * acceleration.row(i) + m_gamma * acceleration.row(i));
    }

    // Update the deformed mesh
    deformedMesh.V = *V;

    return false; // Continue simulation
}


/**
 * Given the edges and vertices positions, adds the stretchings forces to the passed vector.
 */
void DiscreteShell::addStrechingForcesTo(Eigen::MatrixX3d &_forces, const Eigen::MatrixXd *V_) {
    // TODO : make this static or something
    Eigen::MatrixX3d forces_2 = Eigen::MatrixX3d::Zero(V_->rows(), 3);
    // Iterate over each edges
    for (int i = 0; i < E->rows(); i++) {
        // Get the indices of the two vertices of the edge
        int v1 = E->operator()(i, 0);
        int v2 = E->operator()(i, 1);
        // Get the positions of the two vertices
        Eigen::Vector3d p1 = V_->row(v1);
        Eigen::Vector3d p2 = V_->row(v2);
        // Compute the stretching force
        double current_length = (p2 - p1).norm();
        double rest_length = E_length_rest(i);
        double force_magnitude = k_membrane * (current_length - rest_length);
        auto direction = (p2 - p1).normalized();
        Eigen::Vector3d fi = force_magnitude * direction;
        forces_2.row(v1) += fi;
        forces_2.row(v2) -= fi;
    }
    _forces += forces_2;
}


// FULL CHATGPT GENERATED <3
void DiscreteShell::addStretchingHessianTo(
        std::vector<Eigen::Triplet<double>> &triplets,
        const Eigen::MatrixXd *V_
) {
    // Each vertex has 3 DOFs: the global index for dof (x,y,z) of vertex v is:
    // global_dof = v * 3 + component, where component in {0,1,2} for x,y,z.

    for (int i = 0; i < E->rows(); i++) {
        int v1 = (*E)(i, 0);
        int v2 = (*E)(i, 1);

        Eigen::Vector3d p1 = V_->row(v1);
        Eigen::Vector3d p2 = V_->row(v2);

        Eigen::Vector3d d = p2 - p1;
        double r = d.norm();
        double L0 = E_length_rest(i);
        double t = (r - L0);
        double inv_r = 1.0 / r;

        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d ddT = (d * d.transpose()) * (inv_r * inv_r);

        double factor = k_membrane;
        double t_over_r = t * inv_r;

        // Hessian block for a spring:
        // K_block = factor * [ (1 - t/r)*ddT + (t/r)*I ]
        Eigen::Matrix3d K_block = factor * ((1.0 - t_over_r)*ddT + t_over_r*I);

        // Blocks:
        // K_ii = K_block
        // K_jj = K_block
        // K_ij = -K_block
        // K_ji = -K_block
        Eigen::Matrix3d K_ii = K_block;
        Eigen::Matrix3d K_jj = K_block;
        Eigen::Matrix3d K_ij = -K_block;
        Eigen::Matrix3d K_ji = -K_block;

        int row_i = v1 * 3;
        int row_j = v2 * 3;
        int col_i = v1 * 3;
        int col_j = v2 * 3;

        // Helper lambda to add a 3x3 block to the triplets
        auto addBlock = [&](int base_row, int base_col, const Eigen::Matrix3d &M) {
            for (int rr = 0; rr < 3; rr++) {
                for (int cc = 0; cc < 3; cc++) {
                    double val = M(rr, cc);
                    if (val != 0.0) {
                        triplets.emplace_back(base_row + rr, base_col + cc, val);
                    }
                }
            }
        };

        // Add the blocks to the triplet list
        addBlock(row_i, col_i, K_ii);
        addBlock(row_i, col_j, K_ij);
        addBlock(row_j, col_i, K_ji);
        addBlock(row_j, col_j, K_jj);
    }
}



void DiscreteShell::addBendingForcesTo(Eigen::MatrixX3d &forces, const Eigen::MatrixXd *_V, const Eigen::MatrixXi *_E) {

    // std::cout << "Computing bending forces" << std::endl;
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
        forces.row(i) = -Eigen::Map<Eigen::Vector3d>(grad.data());
    }
}

/**
 * Calculates the bending energy for vertex i
 */
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

/**
 * Given a vertice configuration, computes the internal forces.
 */
void DiscreteShell::compute_F_int(Eigen::MatrixX3d &_forces, const Eigen::MatrixXd *_V) const {
    // addStrechingForcesTo(_forces, _V, E);
    // addBendingForces(forces);
}

const Eigen::MatrixXd* DiscreteShell::getPositions() {
    return V;
}

const Eigen::MatrixXi* DiscreteShell::getFaces() {
    return F;
}


