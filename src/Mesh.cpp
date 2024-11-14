#include "Mesh.h"

#include <igl/massmatrix.h>
#include <igl/facet_adjacency_matrix.h>
#include <igl/per_face_normals.h>
#include <igl/edge_flaps.h>
#include <igl/adjacency_matrix.h>



// Constructor for Mesh class
Mesh::Mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_)
    : V(V_), F(F_) {

    Eigen::MatrixXi EI; Eigen::VectorXd EMAP;
    igl::edge_flaps(F_, E, EMAP, EF, EI); // Compute the edges and the edge-face incidence
    igl::massmatrix(vertices, faces, igl::MASSMATRIX_TYPE_DEFAULT, M); // Compute the mass matrix
    incidenceMatrix(F, E, IN);
    
}


void Mesh::calculateDihedralAngle(Eigen::MatrixXd& K, Eigen::VectorXd& angles) {
    igl::per_face_normals(V, F, FN); // Compute per face normals
    for (int i = 0; i < EF.rows(); i++) {
        if (EF(i, 1) == -1) //boundary edge
            angles(i) = 0.0f;
        Eigen::Vector3d v0 = FN.row(EF(i, 0));
        Eigen::Vector3d v1 = FN.row(EF(i, 1));
        angle = acos(v0.dot(v1));
        angles(i) = angle;

        //also determine stiffness
        if (fabs(edge.dihedralAngle) > 0.0) { // Apply a threshold to determine if the edge is a crease
            K(i, i) = 0.5;
        }
    }

}


void incidenceMatrix(const Eigen::MatrixXi& F, const Eigen::MatrixXi& E, Eigen::SparseMatrix<int>& IN) {
    int nV = F.maxCoeff() + 1;
    int nE = E.rows();     

    IN.resize(nV, nE);
    IN.setZero();

    for (int i = 0; i < nE; ++i) {
        int v1 = E(i, 0);
        int v2 = E(i, 1);

        IN.insert(v1, i) = 1;
        IN.insert(v2, i) = 1;
    }
}





