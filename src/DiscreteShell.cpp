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
#include <igl/decimate.h>
#include <filesystem>
#include <igl/upsample.h>
#include <igl/write_triangle_mesh.h>

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>
#include <Eigen/Core>

// Define a macro to enable/disable debug information
// #define DEBUG_VERTEX 0

int num_steps = 0;

// Constructor
DiscreteShell::DiscreteShell()
        : m_dt(SIMULATION_DT), simulation_duration(10.0),
        // bending_stiffness(1.0),
          m_beta(SIMULATION_NEWMARK_BETA), m_gamma(SIMULATION_NEWMARK_GAMMA) {
    F = new Eigen::MatrixXi(0, 3);
    V = new Eigen::MatrixXd(0, 3);
    Velocity = new Eigen::MatrixXd(0, 3);
    Acceleration = new Eigen::MatrixXd(0, 3);
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
void DiscreteShell::initializeFromFile(const std::string &filename) {
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
#ifdef DECIMATE_N_FACES_TARGET_MESH
    igl::decimate(*V, *F, DECIMATE_N_FACES_TARGET_MESH, *V, *F, J);
#endif
#ifdef UPSAMPLE_FACTOR
    // Apply Loop subdivision instead of decimation
    const int NUM_SUBDIVISIONS = 4; // Number of subdivision steps
    Eigen::MatrixXd V_sub;
    Eigen::MatrixXi F_sub;

    // Apply Loop subdivision
    igl::upsample(*V, *F, V_sub, F_sub, UPSAMPLE_FACTOR);
    // Update V and F with the subdivided mesh
    *V = V_sub;
    *F = F_sub;

    // Write to file
    std::string output_file = "../data/rectangle_folded_upsampled_" + std::to_string(UPSAMPLE_FACTOR) + ".off";
    igl::write_triangle_mesh(output_file, V_sub, F_sub);
    std::cout << "Subdivided mesh written to " << output_file << std::endl;
#endif
    V_rest = Eigen::MatrixXd(*V);

    N_VERTICES = V->rows();
    // TODO REMOVE THAT SHIT
    // Set massmatrix
    igl::massmatrix(V_rest, *F, igl::MASSMATRIX_TYPE_VORONOI, M);
    M.setIdentity();

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

    deformedMesh = Mesh(V_rest, *F, *E);
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

    APPLY_INITIAL_CHANGES(V);

    // Set velocity to zero everywhere
    Velocity->resize(V->rows(), 3);
    Velocity->setZero();
    Acceleration->resize(V->rows(), 3);
    Acceleration->setZero();

    // Initialize the solver
    m_solver = new Solver(m_beta, m_gamma, m_dt, &M);
    m_solver->setDiscreteShell(this);
}

bool DiscreteShell::advanceOneStep(int step) {
    // Compute acceleration
    Eigen::MatrixXd new_position = Eigen::MatrixXd::Zero(V->rows(), 3);
    // PRINT_VECTOR(acceleration)
    m_solver->solve(new_position, Velocity, Acceleration, V);

    Eigen::MatrixXd a_new_i =
            (new_position - *V - m_dt * *Velocity - (m_dt * m_dt) * (0.5 - m_beta) * *Acceleration) /
            (m_beta * m_dt * m_dt);
    Eigen::MatrixXd v_new = *Velocity + m_dt * ((1.0 - m_gamma) * *Acceleration + m_gamma * a_new_i);

    *V = new_position;
    *Velocity = v_new;
    *Acceleration = a_new_i;

    // Update the deformed mesh
    deformedMesh.V = *V;

    return false; // Continue simulation
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

void DiscreteShell::addStretchingForcesAndHessianTo_AD_internal(Eigen::MatrixXd &forces, Eigen::SparseMatrix<double> &H,
                                                                const Eigen::MatrixXd *V, bool computeHessian) {
    //3 variables per vertex. Objective per edge, using 3 vertices each :
    // Upda the deformed mesh vertices to updated vertices. This could be a massive
    auto originalV = deformedMesh.V;
    deformedMesh.V = *V;

    //3 variables per vertex. Objective per edge, using 3 vertices each :
    auto func = TinyAD::scalar_function<3>(TinyAD::range(deformedMesh.V.rows()));
    func.add_elements<2>(TinyAD::range(deformedMesh.uE.rows()), [&](auto &element) ->
            TINYAD_SCALAR_TYPE(element) {
        using T = TINYAD_SCALAR_TYPE(element);
        // Variables associated with the edge's vertices
        Eigen::Index edge_idx = element.handle;
        Eigen::RowVector3<T> v0 = element.variables(deformedMesh.E(edge_idx, 0));
        Eigen::RowVector3<T> v1 = element.variables(deformedMesh.E(edge_idx, 1));
        // Compute current length e
        Eigen::Matrix<T, 3, 1> d = (v1 - v0).transpose();
        T current_length = d.norm();
        T rest_length = T(deformedMesh.E_resting_lengths(edge_idx));
        T ratio = current_length / rest_length;
        // Compute the stretching energy, accordin to paper's formula.
        T energy = STRETCHING_STIFFNESS * (1 - ratio) * (1 - ratio) * rest_length;
        return energy;
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
    } else {
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
        forces.row(vertex_idx) += g.segment<3>(3 * vertex_idx);
    }

    // Sets back the original vertices. This is a very dirty design but I'm too scared to modify anything at this point
    deformedMesh.V = originalV;
}


void DiscreteShell::addAreaPreserationForcesAndHessianTo(Eigen::MatrixXd &forces, Eigen::SparseMatrix<double> &H,
                                                         const Eigen::MatrixXd *V) {
    //3 variables per vertex. Objective per edge, using 3 vertices each :
    // Upda the deformed mesh vertices to updated vertices. This could be a massive
    auto originalV = deformedMesh.V;
    deformedMesh.V = *V;

    //3 variables per vertex. Objective per edge, using 3 vertices each :
    auto func = TinyAD::scalar_function<3>(TinyAD::range(deformedMesh.V.rows()));
    func.add_elements<3>(TinyAD::range(deformedMesh.F.rows()), [&](auto &element) ->
            TINYAD_SCALAR_TYPE(element) {
        using T = TINYAD_SCALAR_TYPE(element);
        // Variables associated with the edge's vertices
        Eigen::Index face_idcx = element.handle;
        Eigen::RowVector3<T> v0 = element.variables(deformedMesh.F(face_idcx, 0));
        Eigen::RowVector3<T> v1 = element.variables(deformedMesh.F(face_idcx, 1));
        Eigen::RowVector3<T> v2 = element.variables(deformedMesh.F(face_idcx, 2));
        // Compute current length e
        T rest_area = deformedMesh.F_resting_areas(face_idcx);
        T current_area = computeArea(v0, v1, v2);
        T ratio = current_area / rest_area;
        // Compute the stretching energy, accordin to paper's formula.
        T energy = AREA_PRESERVATION_STIFFNESS * (1 - ratio) * (1 - ratio) * rest_area;
        return energy;
    });

    Eigen::VectorXd x = func.x_from_data([&](int vertex_idx) -> Eigen::VectorXd {
        return deformedMesh.V.row(vertex_idx);
    });

    // Structured binding to unpack the tuple
    auto [f_, g, H_t] = func.eval_with_hessian_proj(x);
    H = H_t;

    for (int vertex_idx = 0; vertex_idx < g.size() / 3; vertex_idx++) {
        forces.row(vertex_idx) += g.segment<3>(3 * vertex_idx);
    }
    // Sets back the original vertices. This is a very dirty design but I'm too scared to modify anything at this point
    deformedMesh.V = originalV;
}


void
DiscreteShell::addBendingForcesAndHessianTo_internal(Eigen::MatrixXd &bending_forces, Eigen::SparseMatrix<double> &H,
                                                     const Eigen::MatrixXd *V, bool computeHessian) {
    // Upda the deformed mesh vertices to updated vertices. This could be a massive
    auto originalV = deformedMesh.V;
    deformedMesh.V = *V;
    //3 variables per vertex. Objective per edge, using 3 vertices each :
    auto func = TinyAD::scalar_function<3>(TinyAD::range(deformedMesh.V.rows()));
    func.add_elements<4>(TinyAD::range(deformedMesh.uE.rows()), [&](auto &element) ->
            TINYAD_SCALAR_TYPE(element) {
        using T =
                TINYAD_SCALAR_TYPE(element);
        // Variables associated with the edge's vertices
        Eigen::Index edge_idx = element.handle;

        Eigen::RowVector3<T> opp0 = element.variables(
                deformedMesh.F(deformedMesh.EF(edge_idx, 0), deformedMesh.EI(edge_idx, 0)));
        Eigen::RowVector3<T> opp1 = element.variables(
                deformedMesh.F(deformedMesh.EF(edge_idx, 1), deformedMesh.EI(edge_idx, 1)));
        Eigen::RowVector3<T> v0 = element.variables(deformedMesh.uE(edge_idx, 0));
        Eigen::RowVector3<T> v1 = element.variables(deformedMesh.uE(edge_idx, 1));
        // Compute dihedral angle, height, and norm
        double angle_u = deformedMesh.getDihedralAngles(edge_idx);
        T angle_d = deformedMesh.computeDihedralAngle(v0, v1, edge_idx, opp0, opp1);
        T height = deformedMesh.computeHeight(v0, v1, edge_idx, opp0, opp1);
        T norm = (v1 - v0).norm();
        // Edge Case
        if (height <= Epsilon)
            return (T) 0;
        // Compute bending energy for this edge
        T bending_energy = BENDING_STIFFNESS * ((angle_u - angle_d) * (angle_u - angle_d) * norm * deformedMesh.stiffness(edge_idx)) / height;
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
    } else {
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

    // Sets back the original vertices. This is a very dirty design
    deformedMesh.V = originalV;
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


void DiscreteShell::add_F_ext(Eigen::MatrixXd &_forces, int step) {
    add_F_ext_internal(_forces, step);
}

const Eigen::MatrixXd *DiscreteShell::getPositions() {
    return V;
}

const Eigen::MatrixXi *DiscreteShell::getFaces() {
    return F;
}



