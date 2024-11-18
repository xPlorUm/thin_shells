#ifndef MESH_H
#define MESH_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>
#include <utility>


// For AD
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
using namespace autodiff;
typedef Eigen::Matrix<dual, Eigen::Dynamic, Eigen::Dynamic> DualMatrix;
typedef Eigen::Matrix<dual, Eigen::Dynamic, 1> DualVector;
typedef Eigen::Matrix<dual, 1, 3> Dual3DVector;


// Mesh class representing a 3D mesh structure
class Mesh {
public:
    //dynamic
    DualMatrix V;  // Matrix storing vertex positions (#V, 3)
    DualMatrix FN; // Normals of each face (#F, 3)
    Eigen::SparseMatrix<dual> M; // Massmatrix

    //static
    Eigen::MatrixXi F;  // Matrix storing indices of vertices forming each face (#F, 3)
    Eigen::MatrixXi uE;  // undirected Edges for each face (#uE, 2)
    Eigen::VectorXi EMAP; // maps each row from E to uE (#F*3, 1)
    Eigen::MatrixXi EF; // edge-to-face incidence matrix (#uE, 2)
    Eigen::MatrixXi EI; // edge-to-vertex incidence matrix (#uE, 2)
    //Eigen::SparseMatrix<double> IC;// Incidence matrix (#V, #E)


    Mesh::Mesh();

    // Constructor to initialize the mesh with vertices, faces
    Mesh(const DualMatrix& V_, const Eigen::MatrixXi& F_);

    // Computes and saves the dihedral angles of the current mesh (#uE, 1) and also the Stiffness Matrix (#uE, 1)
    void calculateDihedralAngle(DualVector& angles, DualVector& stiffness);

    // Computes and saves the third of the average heights of the unique mesh edges into a vector(#uE, 1)
    void Mesh::computeAverageHeights(DualVector& heights);

    // Saves the norms of the undirected Edges of the mesh and saves them in a vector (#uE, 1)
    void Mesh::computeEdgeNorms(DualVector& norms);

private:
    // Computes the height of one face given the indec of the corner {0, 1, 2}
    dual computeFaceHeight(const Eigen::RowVector3i& face, const int corner);

    // Computes the incidence matrix of the undirected Edges
    //void incidenceMatrix(const Eigen::MatrixXi& F, const Eigen::MatrixXi& uE);
};

#endif // MESH_H
