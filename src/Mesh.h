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

constexpr double PI = 3.14159265358979323846;
constexpr double Epsilon = 1e-4;

// Mesh class representing a 3D mesh structure
class Mesh {
public:
    //dynamic
    DualMatrix V;  // Matrix storing vertex positions (#V, 3)
    DualMatrix FN; // Normals of each face (#F, 3)
    Eigen::MatrixXd V_d;
    Eigen::MatrixXd FN_d;
    Eigen::VectorXd stiffness; // Stiffness of each edge (#uE, 1)
    // Resting angles of the mesh.
    Eigen::VectorXd restDihedralAngles;


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

    Mesh();

    // Constructor to initialize the mesh with vertices, faces
    Mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_);

    void calculateDihedralAngles(int vertex_i, DualVector& angles);
    void computeAverageHeights(int i, DualVector& heights);
    void computeEdgeNorms(int i, DualVector& norms);

    // Computes and saves the dihedral angles of the current mesh (#uE, 1) and also the Stiffness Matrix (#uE, 1)
    //void calculateAllDihedralAngles(DualVector& angles, DualVector& stiffness);
    //void getDihedralAngles(int i, DualVector& angles, DualVector& stiffness);
    void calculateAllDihedralAngles(Eigen::VectorXd& angles);
    // Returns all precomputed dihedral angles of the ith vertex to all incident edges
    void getRestingDihedralAngles(int i, Eigen::VectorXd& angles);
private:
    // Computes the height of one face given the indec of the corner {0, 1, 2}
    var computeFaceHeight(const Eigen::RowVector3i& face, const int corner);
    static constexpr double plastic_deformation_threshold = PI / 4.0;

// Templated computeFaceNormal function
    template <typename VectorType>
    void computeFaceNormal(int faceI, VectorType& n) const {
        // Ensure that faceI is within bounds
        if (faceI < 0 || faceI >= F.rows()) {
            throw std::out_of_range("Face index out of range.");
        }
        // Extract the vertex indices for the face
        int vi0 = F(faceI, 0);
        int vi1 = F(faceI, 1);
        int vi2 = F(faceI, 2);

        // Ensure vertex indices are within bounds
        if (vi0 < 0 || vi0 >= V.rows() ||
            vi1 < 0 || vi1 >= V.rows() ||
            vi2 < 0 || vi2 >= V.rows()) {
            throw std::out_of_range("Vertex index out of range.");
        }

        // Extract vertex positions
        VectorType v0 = V.row(vi0).template cast<typename VectorType::Scalar>();
        VectorType v1 = V.row(vi1).template cast<typename VectorType::Scalar>();
        VectorType v2 = V.row(vi2).template cast<typename VectorType::Scalar>();

        // Compute the normal
        n = (v1 - v0).cross(v2 - v0).normalized();
    }

};

#endif // MESH_H
