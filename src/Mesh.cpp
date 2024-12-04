#include "Mesh.h"

#include <igl/facet_adjacency_matrix.h>
#include <igl/per_face_normals.h>
#include <igl/edge_flaps.h>
#include <igl/edges.h>


//template <typename T>
Mesh::Mesh() {
}

// Constructor for Mesh class
Mesh::Mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_)
    : V(V_), F(F_) {
    igl::edge_flaps(F, uE, EMAP, EF, EI); // Compute the edges and the edge-face incidence
    calculateAllDihedralAngles(precomputedDihedralAngles);


    // set stiffness to 1 for all edges
    stiffness = Eigen::VectorXd(uE.rows());
    stiffness.setOnes();
}


// For undeformed Mesh
void Mesh::calculateAllDihedralAngles(Eigen::VectorXd& angles) {
    angles.resize(uE.rows());

    igl::per_face_normals(V, F, FN); // Compute per face normals have to be normalized
    for (int i = 0; i < EF.rows(); i++) {
        if (EF(i, 1) == -1 || EF(i, 0) == -1) { //boundary edge
            angles(i) = 0.0f;
            continue;
        }

        Eigen::Vector3d n0 = FN.row(EF(i, 0));
        Eigen::Vector3d n1 = FN.row(EF(i, 1));
        double cos = n0.dot(n1);
        if (cos >= 1.0f) cos = 1.0f - Epsilon;
        if (cos <= -1.0f) cos = -1.0f + Epsilon;

        double angle = acos(cos);
        angles(i) = angle;

    }
}













