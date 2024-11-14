#ifndef MESH_H
#define MESH_H

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <utility>
#include <autodiff/forward/dual.hpp>


// Mesh class representing a 3D mesh structure
class Mesh {
public:
    //dynamic
    Eigen::MatrixXd V;  // Matrix storing vertex positions (#V, 3)
    Eigen::MatrixXd FN; // Normals of each face (#F, 3)

    // Constructor to initialize the mesh with vertices, faces
    Mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_);

    // Returns the dihedral angles of the current mesh (#E, 1) and also the Stiffness Matrix (#E, #E)
    void calculateDihedralAngle(Eigen::VectorXd& angles, Eigen::MatrixXd& K)

private:
    //static
    Eigen::MatrixXi F;  // Matrix storing indices of vertices forming each face (#F, 3)
    Eigen::MatrixXi E;  // Edges (#E, 2)
    Eigen::MatrixXi EF; // Face-to-edge incidence matrix (#E, 2)
    Eigen::SparseMatrix<int> IN; // Incidence matrix (#V, #E)
};

#endif // MESH_H
