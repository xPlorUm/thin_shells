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
#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
using namespace autodiff;
typedef Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> DualMatrix;
typedef Eigen::Matrix<var, Eigen::Dynamic, 1> DualVector;
typedef Eigen::Matrix<var, 1, 3> Dual3DVector;

// Mesh class representing a 3D mesh structure
class Mesh {
public:
    //dynamic
    DualMatrix V;  // Matrix storing vertex positions (#V, 3)
    DualMatrix FN; // Normals of each face (#F, 3)
    Eigen::MatrixXd V_d;
    Eigen::MatrixXd FN_d;
    Eigen::VectorXd stiffness; // Stiffness of each edge (#uE, 1)

    var Epsilon = 1e-4;

    //static
    Eigen::MatrixXi F;  // Matrix storing indices of vertices forming each face (#F, 3)
    Eigen::MatrixXi uE;  // Undirected Edges for each face (#uE, 2)
    Eigen::VectorXi EMAP; // Maps each row from E to uE (#F*3, 1)
    Eigen::MatrixXi EF; // Edge-to-face incidence matrix (#uE, 2)
    Eigen::MatrixXi EI; // Edge-to-vertex incidence matrix (#uE, 2)

    //face adjacency
    std::vector<std::vector<int>> VF; // List of faces incident to each vertex
    std::vector<std::vector<int>> VFi; // Indices of corners of faces incident to each vertex
    std::vector<std::vector<int>> VE; // VE[i] will contain the edges incident to vertex i


    Eigen::VectorXd dihedralAngles;

    Mesh();

    // Constructor to initialize the mesh with vertices, faces
    Mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_);

    void calculateDihedralAngles(int i, DualVector& angles);
    void computeAverageHeights(int i, DualVector& heights);
    void computeEdgeNorms(int i, DualVector& norms);

    // Computes and saves the dihedral angles of the current mesh (#uE, 1) and also the Stiffness Matrix (#uE, 1)
    //void calculateAllDihedralAngles(DualVector& angles, DualVector& stiffness);
    //void getDihedralAngles(int i, DualVector& angles, DualVector& stiffness);
    void calculateAllDihedralAngles(Eigen::VectorXd& angles);
    // Returns all precomputed dihedral angles of the ith vertex to all incident edges
    void getDihedralAngles(int i, Eigen::VectorXd& angles);


    void calculateDihedralAngles(DualVector& angles);
    // Computes and saves the third of the average heights of the unique mesh edges into a vector(#uE, 1)
    void computeAverageHeights(DualVector& heights);
    // Saves the norms of the undirected Edges of the mesh and saves them in a vector (#uE, 1)
    void computeEdgeNorms(DualVector& norms);

private:
    // Computes the height of one face given the indec of the corner {0, 1, 2}
    var computeFaceHeight(const Eigen::RowVector3i& face, const int corner);
    void computeFaceNormal(int faceI, Dual3DVector& n);
    static constexpr double plastic_deformation_threshold = M_PI / 4.0;

};

#endif // MESH_H
