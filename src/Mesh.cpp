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
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M); // Compute the mass matrix
    //incidenceMatrix(F, uE);
}



void Mesh::calculateDihedralAngle(Eigen::VectorXd& angles, Eigen::VectorXd& stiffness) {
    // Initialize them to the right size
    angles.resize(uE.rows());
    stiffness.resize(uE.rows());

    igl::per_face_normals(V, F, FN); // Compute per face normals
    for (int i = 0; i < EF.rows(); i++) {
        if (EF(i, 1) == -1 || EF(i, 0) == -1) { //boundary edge
            angles(i) = 0.0f;
            continue;
        }
        Eigen::Vector3d n0 = FN.row(EF(i, 0));
        Eigen::Vector3d n1 = FN.row(EF(i, 1));
        double angle = acos(n0.dot(n1));
        angles(i) = angle;

        //also determine stiffness
        if (fabs(angle) > 0.0) { // Apply a threshold to determine if the edge is a crease
            stiffness(i) = 0.5;
        }
    }
}


void Mesh::computeAverageHeights(Eigen::VectorXd& heights) {
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

        double height0 = computeFaceHeight(F.row(face0), corner0);
        double height1 = computeFaceHeight(F.row(face1), corner1);

        // Only compute a third of the average
        heights(i) = (height0 + height1) / 6.f;
    }
}

double Mesh::computeFaceHeight(const Eigen::RowVector3i& face, const int corner) {
    Eigen::Vector3d v0 = V.row(face((corner + 1) % 3));
    Eigen::Vector3d v1 = V.row(face((corner + 2) % 3));
    Eigen::Vector3d v2 = V.row(face(corner));

    // Base of the triangle
    Eigen::Vector3d edge = v1 - v0;

    // Compute the height using the triangle area formula:
    // Area = 0.5 * |edge x (v2 - v0)|, Height = 2 * Area / |edge|
    Eigen::Vector3d crossProduct = edge.cross(v2 - v0);
    double area = 0.5 * crossProduct.norm();
    double height = (2.0 * area) / edge.norm();

    return height;
}

void Mesh::computeEdgeNorms(Eigen::VectorXd& norms) {
    for (int i = 0; i < uE.rows(); i++) {
        int v0 = uE(i, 0);
        int v1 = uE(i, 1);

        Eigen::Vector3d edge = V.row(v1) - V.row(v0);

        norms(i) = edge.norm();
    }
}


//void Mesh::incidenceMatrix(const Eigen::MatrixXi& F, const Eigen::MatrixXi& uE) {
//    int nV = F.maxCoeff() + 1;
//    int nE = uE.rows();     
//
//    IN = Eigen::SparseMatrix<double>(nV, nE);
//
//    for (int i = 0; i < nE; ++i) {
//        int v1 = uE(i, 0);
//        int v2 = uE(i, 1);
//
//        IN.insert(v1, i) = 1.0f;
//        IN.insert(v2, i) = 1.0f;
//    }
//}





