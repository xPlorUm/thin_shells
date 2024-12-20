#include "Mesh.h"
#include "common.h"

#include <igl/per_face_normals.h>
#include <igl/edge_flaps.h>

#include <utility>


//template <typename T>
Mesh::Mesh() {
}

// Constructor for Mesh class
Mesh::Mesh(Eigen::MatrixXd V_, const Eigen::MatrixXi &F_, const Eigen::MatrixXi &E_)
        : V(std::move(V_)), F(F_), E(E_) {
    igl::edge_flaps(F, uE, EMAP, EF, EI); // Compute the edges and the edge-face incidence
    VVE = Eigen::MatrixXi::Zero(V.rows(), V.rows());
    for (int i = 0; i < uE.rows(); i++) {
        VVE(uE(i, 0), uE(i, 1)) = i;
        VVE(uE(i, 1), uE(i, 0)) = i;
    }
    stiffness = Eigen::VectorXd(uE.rows());
    stiffness.setOnes();
    calculateAllDihedralAngles(dihedralAngles);
    // set stiffness to 1 for all edges
    // Edge length :
    igl::edge_lengths(V, E, E_resting_lengths);
    // Compute all resting faces area
    F_resting_areas.resize(F.rows());
    for (int i = 0; i < F.rows(); i++) {
        F_resting_areas(i) = computeArea(V.row(F(i, 0)), V.row(F(i, 1)), V.row(F(i, 2)));
    }

}


// For undeformed Mesh
void Mesh::calculateAllDihedralAngles(Eigen::VectorXd &angles) {
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

/**
 * Get the precomputed dihedral angles of the edge at idx
 */
double Mesh::getDihedralAngles(int idx) {
    return dihedralAngles[idx];
}
