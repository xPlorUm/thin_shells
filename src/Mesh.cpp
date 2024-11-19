#include "Mesh.h"

#include <igl/massmatrix.h>
#include <igl/facet_adjacency_matrix.h>
#include <igl/per_face_normals.h>
#include <igl/edge_flaps.h>


Mesh::Mesh() {
}

// Constructor for Mesh class
Mesh::Mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_)
    : V(V_), F(F_) {
    igl::edge_flaps(F, uE, EMAP, EF, EI); // Compute the edges and the edge-face incidence
    igl::massmatrix(V_, F, igl::MASSMATRIX_TYPE_DEFAULT, M); // Compute the mass matrix
}


void Mesh::calculateDihedralAngle(DualVector& angles, DualVector& stiffness) {
    // Initialize them to the right size
    angles.resize(uE.rows());
    stiffness.resize(uE.rows());

    igl::per_face_normals(V, F, FN); // Compute per face normals
    for (int i = 0; i < EF.rows(); i++) {
        if (EF(i, 1) == -1 || EF(i, 0) == -1) { //boundary edge
            angles(i) = 0.0f;
            continue;
        }
        Dual3DVector n0 = FN.row(EF(i, 0));
        Dual3DVector n1 = FN.row(EF(i, 1));
        var angle = acos(n0.dot(n1));
        angles(i) = angle;

        //also determine stiffness
        if (abs(angle) > 0.0) { // Apply a threshold to determine if the edge is a crease
            stiffness(i) = var(0.5f);
        }
    }
}


void Mesh::computeAverageHeights(DualVector& heights) {
    // Initialize to the right size
    heights.resize(uE.rows());

    for (int i = 0; i < EF.rows(); i++) {
        if (EF(i, 1) == -1 || EF(i, 0) == -1) { //boundary edge
            heights(i) = 0.0f;
            continue;
        }
        // Get the indices of the adjacent faces
        int face0 = EF(i, 0);
        int face1 = EF(i, 1);

        // Get the coner point adjacent to the edge
        int corner0 = EI(i, 0);
        int corner1 = EI(i, 1);

        var height0 = computeFaceHeight(F.row(face0), corner0);
        var height1 = computeFaceHeight(F.row(face1), corner1);

        // Only compute a third of the average
        heights(i) = (height0 + height1) / 6.f;
    }
}

var Mesh::computeFaceHeight(const Eigen::RowVector3i& face, const int corner) {
    Dual3DVector v0 = V.row(face((corner + 1) % 3));
    Dual3DVector v1 = V.row(face((corner + 2) % 3));
    Dual3DVector v2 = V.row(face(corner));

    // Base of the triangle not including the corner
    Dual3DVector edge = v1 - v0;

    // Compute the height using the triangle area formula:
    // Area = 0.5 * |edge x (v2 - v0)|, Height = 2 * Area / |edge|
    var area = 0.5 * edge.cross(v2 - v0).norm();
    var height = (2.0 * area) / edge.norm();

    return height;
}

void Mesh::computeEdgeNorms(DualVector& norms) {
    for (int i = 0; i < uE.rows(); i++) {
        int v0 = uE(i, 0);
        int v1 = uE(i, 1);

        Dual3DVector edge = V.row(v1) - V.row(v0);

        norms(i) = edge.norm();
    }
}






