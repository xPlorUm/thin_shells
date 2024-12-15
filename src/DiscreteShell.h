#ifndef DISCRETE_SHELL_H
#define DISCRETE_SHELL_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "Mesh.h"
#include <iostream>
#include <igl/massmatrix.h>
#include <fstream>
#include "Solver.h"

class FileDebugger {
private :
    std::ofstream dataFile;
public:
    FileDebugger(const std::string &filename) {
        dataFile.open(filename);
    }

    // Get the stream so one can write on it
    std::ofstream &getStream() {
        return dataFile;
    }

    void close() {
        dataFile.close();
    }


    ~FileDebugger() {
        dataFile.close();
    }
};

class DiscreteShell {
public:
    // Constructor
    DiscreteShell();

    // Destructor
    ~DiscreteShell();

    int N_VERTICES;

    void initializeFromFile(const std::string &filename);

    // Advance one time step
    bool advanceOneStep(int step);

    // Given an index of a vertex, print some stuff. Very much useful for debugging
    void debug_vertex(int vertex_index) {
        std::cout << "Vertex " << vertex_index << std::endl;
        // Velocity
        std::cout << "Velocity : " << Velocity->row(vertex_index) << std::endl;
        // Position
        std::cout << "Position : " << V->row(vertex_index) << std::endl;

        // Compute the bending forces
        Eigen::MatrixXd bending_forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
        addBendingForcesTo(bending_forces, V);
        std::cout << "Bending Forces : " << bending_forces.row(vertex_index) << std::endl;
        Eigen::SparseMatrix<double> H = Eigen::SparseMatrix<double>(V->rows() * 3, V->rows() * 3);
        Eigen::MatrixXd stretching_forces = Eigen::MatrixX3d::Zero(V->rows(), 3);
        addStretchingForcesAndHessianTo_AD(stretching_forces, H, V);
        std::cout << "Stretching forces with AD : " << stretching_forces.row(vertex_index) << std::endl;

    }

    // Returns the positions to get drawn
    // TODO : this should be a pointer
    const Eigen::MatrixXd *getPositions();

    const Eigen::MatrixXi *getFaces();

    void compute_F_int(Eigen::MatrixX3d &_forces, const Eigen::MatrixXd *_V) const;

    void add_F_ext(Eigen::MatrixXd &_forces);

    void addBendingForcesTo(Eigen::MatrixXd &bending_forces, const Eigen::MatrixXd *V);

    void
    addBendingForcesAndHessianTo(Eigen::MatrixXd &bending_forces, Eigen::SparseMatrix<double> &H,
                                 const Eigen::MatrixXd *V);

    void addStretchingForcesTo_AD(Eigen::MatrixXd &forces, const Eigen::MatrixXd *V) {
        Eigen::SparseMatrix<double> H = Eigen::SparseMatrix<double>(V->rows() * 3, V->rows() * 3);
        addStretchingForcesAndHessianTo_AD_internal(forces, H, V, false);
    }

    void addStretchingForcesAndHessianTo_AD(Eigen::MatrixXd &forces, Eigen::SparseMatrix<double> &H,
                                            const Eigen::MatrixXd *V) {
        addStretchingForcesAndHessianTo_AD_internal(forces, H, V, true);
    }

private:

    // Physical properties
    double m_dt; // Time step
    double simulation_duration; // Total simulation duration
    // double bending_stiffness; // Stiffness for bending energy

    // State variables
    Mesh deformedMesh;
    Mesh undeformedMesh;

    Eigen::MatrixXi *F; // Faces of the shell
    Eigen::MatrixXd *V; // Current mesh,
    Eigen::MatrixXd V_rest; // Undeformed mesh
    Eigen::MatrixXd *Velocity; // Velocity of the shell (point-wise)
    Eigen::MatrixXd *Acceleration; // Acceleration of the shell (point-wise)
    Eigen::MatrixXi *E; // Edges of the shell
    Eigen::VectorXd E_length_rest; // Rest length of edges
    Eigen::SparseMatrix<double> M; // Massmatrix
    Eigen::SparseMatrix<double> M_i;
    Eigen::MatrixXd forces; // Forces applied point-wise.
    Eigen::MatrixX3d bending_forces; // Bending forces applied point-wise.
    int k_membrane = 10; // Membrane stiffness
    double stiffness_damping = 0.5; // Stiffness damping coefficient
    double mass_damping = 0.5; // Mass damping coefficient


    // Newmark-specific parameters
    double m_beta = 0.25; // Newmark parameter (default: 0.25 for implicit integration)
    double m_gamma = 0.5; // Newmark parameter (default: 0.5 for implicit integration)

    void computeDampingForces(Eigen::MatrixX3d &damping_forces); // Compute damping forces

    FileDebugger fileDebugger = FileDebugger("vertex_data.csv");

    Solver *m_solver;

    void addBendingForcesAndHessianTo_internal(Eigen::MatrixXd &bending_forces, Eigen::SparseMatrix<double> &H,
                                               const Eigen::MatrixXd *V, bool computeHessian = false);


    void addStretchingForcesAndHessianTo_AD_internal(Eigen::MatrixXd &forces, Eigen::SparseMatrix<double> &H,
                                                     const Eigen::MatrixXd *V, bool computeHessian);

};


#endif // DISCRETE_SHELL_H
