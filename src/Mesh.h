#ifndef MESH_H
#define MESH_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>
#include <utility>
//#include <autodiff/forward/dual.hpp>


// Mesh class representing a 3D mesh structure
class Mesh {
public:
    //dynamic
    Eigen::MatrixXd V;  // Matrix storing vertex positions (#V, 3)
    Eigen::MatrixXd FN; // Normals of each face (#F, 3)
    Eigen::SparseMatrix<double> M; // Massmatrix

    //static
    Eigen::MatrixXi F;  // Matrix storing indices of vertices forming each face (#F, 3)
    Eigen::MatrixXi uE;  // undirected Edges for each face (#uE, 2)
    Eigen::VectorXi EMAP; // maps each row from E to uE (#F*3, 1)
    Eigen::MatrixXi EF; // edge-to-face incidence matrix (#uE, 2)
    Eigen::MatrixXi EI; // edge-to-vertex incidence matrix (#uE, 2)
    //Eigen::SparseMatrix<double> IN;// Incidence matrix (#V, #E)


    Mesh::Mesh();

    // Constructor to initialize the mesh with vertices, faces
    Mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_);

    // Computes and saves the dihedral angles of the current mesh (#uE, 1) and also the Stiffness Matrix (#uE, 1)
    void calculateDihedralAngle(Eigen::VectorXd& angles, Eigen::VectorXd& stiffness);

    // Computes and saves the third of the average heights of the unique mesh edges into a vector(#uE, 1)
    void Mesh::computeAverageHeights(Eigen::VectorXd& heights);

    // Saves the norms of the undirected Edges of the mesh and saves them in a vector (#uE, 1)
    void Mesh::computeEdgeNorms(Eigen::VectorXd& norms);

private:
    // Computes the height of one face given the indec of the corner {0, 1, 2}
    double computeFaceHeight(const Eigen::RowVector3i& face, const int corner);

    // Computes the incidence matrix of the undirected Edges
    void incidenceMatrix(const Eigen::MatrixXi& F, const Eigen::MatrixXi& uE);
};

#endif // MESH_H
