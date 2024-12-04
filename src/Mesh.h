#ifndef MESH_H
#define MESH_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>
#include <utility>




// For AD
constexpr double PI = 3.14159265358979323846;

constexpr double Epsilon = 1e-6;

// Mesh class representing a 3D mesh structure
class Mesh {
public:

    //dynamic
    Eigen::MatrixXd V;  // Matrix storing vertex positions (#V, 3)
    Eigen::MatrixXd FN; // Normals of each face (#F, 3)


    Eigen::VectorXd stiffness; // Stiffness of each edge (#uE, 1)

    //static
    Eigen::MatrixXi F;  // Matrix storing indices of vertices forming each face (#F, 3)
    Eigen::MatrixXi uE;  // Undirected Edges for each face (#uE, 2)
    Eigen::VectorXi EMAP; // Maps each row from E to uE (#F*3, 1)
    Eigen::MatrixXi EF; // Edge-to-face incidence matrix (#uE, 2)
    Eigen::MatrixXi EI; // Edge-to-vertex incidence matrix (#uE, 2)


    Mesh();

    // Constructor to initialize the mesh with vertices, faces
    Mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_);

    Eigen::VectorXd precomputedDihedralAngles; // Precomputed in Mesh() dihedralAngles of size (#uE, 1)


    template <typename Scalar>
    Scalar computeHeight(Eigen::Matrix<Scalar, 1, 3>& v0, Eigen::Matrix<Scalar, 1, 3>& v1, int e_idx) {
        if (EF(e_idx, 1) == -1 || EF(e_idx, 0) == -1) { //boundary edge
            return 0.0f;
        }

        // Get the indices of the adjacent faces
        int face0 = EF(e_idx, 0);
        int face1 = EF(e_idx, 1);

        // Get the coner point adjacent to the edge
        Eigen::RowVector3d corner0 = V.row(F(face0, EI(e_idx, 0)));
        Eigen::RowVector3d corner1 = V.row(F(face1, EI(e_idx, 1)));

        Scalar height0 = computeFaceHeight(v0, v1, corner0);
        Scalar height1 = computeFaceHeight(v0, v1, corner1);

        // Only compute a third of the average
        return (height0 + height1) / 6.f;
    }

    template <typename Scalar>
    Scalar computeFaceHeight(Eigen::Matrix<Scalar, 1, 3>& v0, Eigen::Matrix<Scalar, 1, 3>& v1,
        Eigen::RowVector3d& corner) {
        // Base of the triangle not including the corner
        Eigen::Matrix<Scalar, 1, 3> edge = v1 - v0;

        // Compute the height using the triangle area formula:
        // Area = 0.5 * |edge x (v2 - v0)|, Height = 2 * Area / |edge|
        Scalar area = 0.5 * edge.cross(corner - v0).norm();
        Scalar norm = edge.norm();
        if (norm < Epsilon) 
            return 0;
        
        return (2.0 * area) / norm;
    }


    template <typename Scalar>
    Scalar computeDihedralAngle(Eigen::Matrix<Scalar, 1, 3>& v0,
        Eigen::Matrix<Scalar, 1, 3>& v1, int e_idx) {
        if (EF(e_idx, 1) == -1 || EF(e_idx, 0) == -1) { //boundary edge
            return 0.0f;
        }

        // Get the indices of the adjacent faces
        int face0 = EF(e_idx, 0);
        int face1 = EF(e_idx, 1);

        // Get the coner point adjacent to the edge
        Eigen::RowVector3d corner0 = V.row(F(face0, EI(e_idx, 0)));
        Eigen::RowVector3d corner1 = V.row(F(face1, EI(e_idx, 1)));

        Eigen::Matrix<Scalar, 1, 3> n0, n1;
        computeNormal(v0, v1, corner0, n0);
        computeNormal(v0, v1, corner1, n1);

        Scalar cos = n0.dot(n1);

        if (cos >= 1.0f) cos = 1.0f - Epsilon;
        if (cos <= -1.0f) cos = -1.0f + Epsilon;

        // TODO not so sure if correct what about degree of 180+?
        Scalar angle = acos(cos);
        // Determine sign of the angle using the cross product of normals
        Eigen::Matrix<Scalar, 1, 3> edge_dir = (v1 - v0).normalized();
        Eigen::Matrix<Scalar, 1, 3> cross_n0_n1 = n0.cross(n1);
        if (cross_n0_n1.dot(edge_dir) < 0) {
            angle = -angle; // Negative angle for >180°
        }

        return angle;
    }

    template <typename Scalar>
    void computeNormal(Eigen::Matrix<Scalar, 1, 3>& v0, Eigen::Matrix<Scalar, 1, 3>& v1,
        Eigen::Matrix<double, 1, 3>& corner, Eigen::Matrix<Scalar, 1, 3>& res) {

        Eigen::Matrix<Scalar, 1, 3> edge1 = v1 - v0;
        Eigen::Matrix<Scalar, 1, 3> edge2 = corner - v0;

        res = edge1.cross(edge2);
        Scalar norm = res.norm();
        if (norm > Epsilon) {
            res /= norm;
        }
        else {
            // Handle degenerate triangles (zero-area)
            res.setZero();
        }
    }


    // Computes and saves the dihedral angles of the current mesh (#uE, 1) and also the Stiffness Matrix (#uE, 1)
    void calculateAllDihedralAngles(Eigen::VectorXd& angles);

private:
    // Computes the height of one face given the index of the corner {0, 1, 2}
    static constexpr double plastic_deformation_threshold = PI / 4.0;

};

#endif // MESH_H
