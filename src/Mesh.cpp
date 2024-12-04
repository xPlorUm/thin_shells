#include "Mesh.h"

#include <igl/facet_adjacency_matrix.h>
#include <igl/per_face_normals.h>
#include <igl/edge_flaps.h>
#include <igl/edges.h>


Mesh::Mesh() {
}

// Constructor for Mesh class
Mesh::Mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_)
    : V(V_), F(F_), V_d(V_) {
    igl::edge_flaps(F, uE, EMAP, EF, EI); // Compute the edges and the edge-face incidence
    calculateAllDihedralAngles(preComputedDihedralAngles);


    // Prepare the adjacency list
    VE = std::vector<std::vector<int>>(V.rows());
    for (int e = 0; e < uE.rows(); ++e) {
        int v0 = uE(e, 0);
        int v1 = uE(e, 1);

        VE[uE(e, 0)].push_back(e);
        VE[uE(e, 1)].push_back(e);
    }

    // set stiffness to 1 for all edges
    stiffness = Eigen::VectorXd(uE.rows());
    stiffness.setOnes();
}






// For undeformed Mesh

void Mesh::calculateAllDihedralAngles(Eigen::VectorXd& angles) {
    angles.resize(uE.rows());

    igl::per_face_normals(V_d, F, FN_d); // Compute per face normals have to be normalized
    for (int i = 0; i < EF.rows(); i++) {
        if (EF(i, 1) == -1 || EF(i, 0) == -1) { //boundary edge
            angles(i) = 0.0f;
            continue;
        }

        Eigen::Vector3d n0 = FN_d.row(EF(i, 0));
        Eigen::Vector3d n1 = FN_d.row(EF(i, 1));
        double cos = n0.dot(n1);
        if (cos >= 1.0f) cos = 1.0f - Epsilon;
        if (cos <= -1.0f) cos = -1.0f + Epsilon;

        double angle = acos(cos);
        angles(i) = angle;

    }
}

void Mesh::getDihedralAngles(int i, Eigen::VectorXd& angles) {
    angles.resize(VE[i].size());
    int ni = 0;
    for (int n : VE[i]) {
        angles(ni) = dihedralAngles(n);
        ni++;
    }
}



// For all vertices
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


void Mesh::computeEdgeNorms(Eigen::VectorXd& norms) {
    for (int i = 0; i < uE.rows(); i++) {
        int v0 = uE(i, 0);
        int v1 = uE(i, 1);

        Eigen::Vector3d edge = V.row(v1) - V.row(v0);

        norms(i) = edge.norm();
    }
}

double Mesh::computeFaceHeight(const Eigen::RowVector3i& face, const int corner) {
    Eigen::Vector3d v0 = V.row(face((corner + 1) % 3));
    Eigen::Vector3d v1 = V.row(face((corner + 2) % 3));
    Eigen::Vector3d v2 = V.row(face(corner));

    // Base of the triangle not including the corner
    Eigen::Vector3d edge = v1 - v0;

    // Compute the height using the triangle area formula:
    // Area = 0.5 * |edge x (v2 - v0)|, Height = 2 * Area / |edge|
    double area = 0.5 * edge.cross(v2 - v0).norm();
    double norm = edge.norm();
    if (norm < Epsilon) return Epsilon;
    else return (2.0 * area) / norm;
}





