/*
Boilerplate code cleaned up from the sacred texts
Maybe we should change the function names since they're currently the same as the sacred texts...
*/

#include "DiscreteShell.h"
#include "igl/read_triangle_mesh.h"
#include "common.h"
#include "constants.h"
#include <iostream>
#include <igl/edges.h> // Include this header to use igl::edges
#include <filesystem>
#include <igl/decimate.h>
#include <fstream>
#include <iostream>

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

// Define a macro to enable/disable debug information
// #define DEBUG_VERTEX 0

int num_steps = 0;

// Constructor
DiscreteShell::DiscreteShell()
    : m_dt(SIMULATION_DT), simulation_duration(10.0),
    // bending_stiffness(1.0),
      m_beta(SIMULATION_NEWMARK_BETA), m_gamma(SIMULATION_NEWMARK_GAMMA)
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
    igl::decimate(*V, *F, N_FACES_MESH, *V, *F, J);
    V_rest = Eigen::MatrixXd(*V);
    deformedMesh = Mesh(V_rest, *F);

    // Set only the y-coordinates to 0 for the undeformed mesh
    Eigen::MatrixXd V_undef = V_rest; // Copy vertices
    for (int i = 0; i < V_undef.rows(); ++i) {
        V_undef(i, 1) = 0.0; // Set the y-coordinate (index 1) to 0
    }
    N_VERTICES = V->rows();
    // TODO REMOVE THAT SHIT
    undeformedMesh = Mesh(V_undef, *F);
    // Set massmatrix
    igl::massmatrix(V_rest, *F, igl::MASSMATRIX_TYPE_VORONOI, M);
    // M.setIdentity();

    // Set up the inverse of M
    // to avoid instability
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

    // ------ Artifical little change to the mesh -----------
    V->row(42)  += Eigen::Vector3d(0.0, 0, 30);
    V->row(43)  += Eigen::Vector3d(0.0, 0, 30);
    V->row(0)   += Eigen::Vector3d(0.0, 0, 30);
    V->row(1)   += Eigen::Vector3d(0.0, 0, 30);
    V->row(2)   += Eigen::Vector3d(0.0, 0, 30);
    // M.coeffRef(8, 8) = 1e50;
    // ------------------------------------------------------

    // Set velocity to zero everywhere
    Velocity->resize(V->rows(), 3);
    Velocity->setZero();

    // Initialize the solver
    m_solver = new Solver(m_beta, m_gamma, m_dt, &M);
    m_solver->setDiscreteShell(this);
}

bool DiscreteShell::advanceOneStep(int step) {
    // Reset forces
    forces.setZero(V->rows(), 3);
    // Compute forces
    // addStrechingForcesTo(forces, V);
    addBendingForcesTo(forces, V);
    // Compute damping forces
    Eigen::MatrixX3d damping_forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
    computeDampingForces(damping_forces);
    // forces += damping_forces;

    // Add external forces (only gravity for now)
    // TODO /!\ : FIgure out the + and - shit
    add_F_ext(forces);

    // Compute acceleration
    Eigen::MatrixXd acceleration = Eigen::MatrixXd::Zero(V->rows(), 3);
    for (int i = 0; i < Velocity->rows(); ++i) {
        double mass_i = M.coeff(i, i); // Mass for particle i
        // acceleration.row(i) = forces.row(i);// / mass_i;
        // Eigen::Vector3d force = forces.row(i);
        // Eigen::Vector3d velocity = Velocity->row(i);

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
    Eigen::MatrixXd new_position = Eigen::MatrixXd::Zero(V->rows(), 3);
    // PRINT_VECTOR(acceleration)
    m_solver->solve(new_position, Velocity, &acceleration, V);

    Eigen::MatrixXd a_new_i =
            (new_position - *V - m_dt * *Velocity - (m_dt * m_dt) * (0.5 - m_beta) * acceleration) / (m_beta * m_dt * m_dt);
    Eigen::MatrixXd v_new = *Velocity + m_dt * ((1.0 - m_gamma) * acceleration + m_gamma * a_new_i);

    *V = new_position;
    *Velocity = v_new;
    acceleration = a_new_i;

    // Update the deformed mesh
    deformedMesh.V = *V;

    return false; // Continue simulation
}


/**
 * Given the edges and vertices positions, adds the stretchings forces to the passed vector.
 */
void DiscreteShell::addStrechingForcesTo(Eigen::MatrixXd &_forces, const Eigen::MatrixXd *V_) {
    // TODO : make this static or something
    Eigen::MatrixXd forces_2 = Eigen::MatrixXd::Zero(V_->rows(), 3);
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
        double force_magnitude = STRETCHING_STIFFNESS * (current_length - rest_length);
        Eigen::VectorXd direction = (p2 - p1).normalized();
        Eigen::VectorXd fi = force_magnitude * direction;
        _forces.row(v1) += fi;
        _forces.row(v2) -= fi;
    }
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
}



void DiscreteShell::addBendingForcesAndHessianTo_internal(Eigen::MatrixXd &bending_forces, Eigen::SparseMatrix<double> &H,
                                                          const Eigen::MatrixXd *V, bool computeHessian) {
    //3 variables per vertex. Objective per edge, using 3 vertices each :
    auto func = TinyAD::scalar_function<3>(TinyAD::range(deformedMesh.V.rows()));
    func.add_elements<2>(TinyAD::range(deformedMesh.uE.rows()), [&](auto &element) ->
    TINYAD_SCALAR_TYPE(element)
    {
        using T =
        TINYAD_SCALAR_TYPE(element);

        // Variables associated with the edge's vertices
        Eigen::Index edge_idx = element.handle;
        Eigen::RowVector3<T> v0 = element.variables(deformedMesh.uE(edge_idx, 0));
        Eigen::RowVector3<T> v1 = element.variables(deformedMesh.uE(edge_idx, 1));
        // Compute dihedral angle, height, and norm
        double angle_u = deformedMesh.getDihedralAngles(edge_idx);
        T angle_d = deformedMesh.computeDihedralAngle(v0, v1, edge_idx);
        angle_d = deformedMesh.computeDihedralAngle(v0, v1, edge_idx);
        T height = deformedMesh.computeHeight(v0, v1, edge_idx);
        T norm = (v1 - v0).norm();
        double stiffness = 1.0f;
        // Edge Case
        if (height <= Epsilon)
            return (T) 0.0;
        // Compute bending energy for this edge
        T bending_energy = ((angle_u - angle_d) * (angle_u - angle_d) * norm * stiffness) / height;
        return bending_energy;
    });


    Eigen::VectorXd x = func.x_from_data([&](int vertex_idx) -> Eigen::VectorXd {
        return deformedMesh.V.row(vertex_idx);
    });


    double f;                          // Assuming f is a double
    Eigen::VectorXd g;                 // Assuming g is an Eigen vector
    if (computeHessian) {
        // Structured binding to unpack the tuple
        auto [f_, g_, H_t] = func.eval_with_hessian_proj(x);
        f = f_;
        g = g_;
        H = H_t;
    }
    else {
        // Structured binding for the gradient-only case
        auto [f_, g_] = func.eval_with_gradient(x);
        // Assign to the pre-declared variables
        f = f_;
        g = g_;
    }
    // We get the output force f, the gradient g and the Hessian H
    // With the Hessian H it is way slower now
    //auto [f, g, H] = func.eval_with_hessian_proj(x);
    for (int vertex_idx = 0; vertex_idx < g.size() / 3; vertex_idx++) {
        bending_forces.row(vertex_idx) += g.segment<3>(3 * vertex_idx);
    }
}

void DiscreteShell::addBendingForcesTo(Eigen::MatrixXd &bending_forces, const Eigen::MatrixXd *V) {
    // Dummy parameter
    Eigen::SparseMatrix<double> t = Eigen::SparseMatrix<double>(V->rows() * 3, V->rows() * 3);
    addBendingForcesAndHessianTo_internal(bending_forces, t, V, false);
}

void DiscreteShell::addBendingForcesAndHessianTo(Eigen::MatrixXd &bending_forces,
                                                          Eigen::SparseMatrix<double> &H, const Eigen::MatrixXd *V) {
    addBendingForcesAndHessianTo_internal(bending_forces, H, V, true);
}



void DiscreteShell::add_F_ext(Eigen::MatrixXd &_forces) {
    // Just gravity for now
    // Substract gravity :
    Eigen::RowVector3d gravity(0, -9.81, 0);
    for (int i = 0; i < V->rows(); i++) {
        // _forces.row(i).y() -= 9.81;
    }
    // Wind only on vertex 42
    Eigen::RowVector3d wind(0, 0, 100  );
    // _forces.row(42) += wind;
}

const Eigen::MatrixXd* DiscreteShell::getPositions() {
    return V;
}

const Eigen::MatrixXi* DiscreteShell::getFaces() {
    return F;
}



