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
    calculateAllDihedralAngles(dihedralAngles);


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

//Vertex wise computation
void Mesh::calculateDihedralAngles(int i, Eigen::VectorXd& angles) {
    std::vector<int> incEdges = VE[i];
    angles.setZero(incEdges.size());
    int ni = 0;
    for (int edgeI : incEdges) {
        if (EF(edgeI, 1) == -1 || EF(edgeI, 0) == -1) { //boundary edge
            ni++;
            continue;
        }
         
        Eigen::Vector3d n0, n1;
        computeFaceNormal(EF(edgeI, 0), n0);
        computeFaceNormal(EF(edgeI, 1), n1);
        double cos = n0.dot(n1);
        // Clamp cos value not derivable at cos == 1.0 or -1.0f
        if (cos >= val(1.0f)) cos = 1.0f - Epsilon;
        if (cos <= val(-1.0f)) cos = -1.0f + Epsilon;
        
        angles(ni) = acos(cos);

        //also determine stiffness
        if (abs(angles(ni)) > plastic_deformation_threshold) { // Apply a threshold to determine if the edge is a crease
            stiffness(ni) = 0.5;
            // TODO: change rest angle
        }
        ni++;
    }
}

void Mesh::computeFaceNormal(int faceI, Eigen::Vector3d& n) {
    Eigen::Vector3d v0 = V.row(F(faceI, 0));
    Eigen::Vector3d v1 = V.row(F(faceI, 1));
    Eigen::Vector3d v2 = V.row(F(faceI, 2));
    n = (v1 - v0).cross(v2 - v0).normalized();
}

void Mesh::computeAverageHeights(int i, Eigen::VectorXd& heights) {
    std::vector<int> incEdges = VE[i];
    heights.setZero(incEdges.size());
    int ni = 0;
    for (int edgeI : incEdges) {
        if (EF(edgeI, 1) == -1 || EF(edgeI, 0) == -1) { //boundary edge
            heights(ni) = Epsilon;
            ni++;
            continue;
        }
        // Get the indices of the adjacent faces
        int face0 = EF(edgeI, 0);
        int face1 = EF(edgeI, 1);

        // Get the coner point adjacent to the edge
        int corner0 = EI(edgeI, 0);
        int corner1 = EI(edgeI, 1);

        double height0 = computeFaceHeight(F.row(face0), corner0);
        double height1 = computeFaceHeight(F.row(face1), corner1);

        // Only compute a third of the average
        double res = (height0 + height1) / 6.f;
        heights(ni) = res;
        ni++;
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

void Mesh::computeEdgeNorms(int i, Eigen::VectorXd& norms) {
    std::vector<int> incEdges = VE[i];
    norms.setZero(incEdges.size());
    int ni = 0;
    for (int edgeI : incEdges) {
        int v0 = uE(edgeI, 0);
        int v1 = uE(edgeI, 1);

        Eigen::Vector3d edge = V.row(v1) - V.row(v0);
        norms(ni) = edge.norm() * stiffness(edgeI);
        ni++;
    }
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

// For AD
//void Mesh::calculateDihedralAngles(int i, DualVector& angles) {
//    std::vector<int> incEdges = VE[i];
//    angles.setZero(incEdges.size());
//    int ni = 0;
//    for (int edgeI : incEdges) {
//        if (EF(edgeI, 1) == -1 || EF(edgeI, 0) == -1) { //boundary edge
//            ni++;
//            continue;
//        }
//
//        Dual3DVector n0, n1;
//        computeFaceNormal(EF(edgeI, 0), n0);
//        computeFaceNormal(EF(edgeI, 1), n1);
//        var cos = n0.dot(n1);
//        // Clamp cos value not derivable at cos == 1.0 or -1.0f
//        if (cos >= val(1.0f)) cos = 1.0f - Epsilon;
//        if (cos <= val(-1.0f)) cos = -1.0f + Epsilon;
//
//        angles(ni) = acos(cos);
//
//        //also determine stiffness
//        if (abs(angles(ni)) > plastic_deformation_threshold) { // Apply a threshold to determine if the edge is a crease
//            stiffness(ni) = 0.5;
//            // TODO: change rest angle
//        }
//        ni++;
//    }
//}
//
//void Mesh::computeFaceNormal(int faceI, Dual3DVector& n) {
//    Dual3DVector v0 = V.row(F(faceI, 0));
//    Dual3DVector v1 = V.row(F(faceI, 1));
//    Dual3DVector v2 = V.row(F(faceI, 2));
//    n = (v1 - v0).cross(v2 - v0).normalized();
//}
//
//void Mesh::computeAverageHeights(int i, DualVector& heights) {
//    std::vector<int> incEdges = VE[i];
//    heights.setZero(incEdges.size());
//    int ni = 0;
//    for (int edgeI : incEdges) {
//        if (EF(edgeI, 1) == -1 || EF(edgeI, 0) == -1) { //boundary edge
//            heights(ni) = Epsilon;
//            ni++;
//            continue;
//        }
//        // Get the indices of the adjacent faces
//        int face0 = EF(edgeI, 0);
//        int face1 = EF(edgeI, 1);
//
//        // Get the coner point adjacent to the edge
//        int corner0 = EI(edgeI, 0);
//        int corner1 = EI(edgeI, 1);
//
//        var height0 = computeFaceHeight(F.row(face0), corner0);
//        var height1 = computeFaceHeight(F.row(face1), corner1);
//
//        // Only compute a third of the average
//        var res = (height0 + height1) / 6.f;
//        heights(ni) = res;
//        ni++;
//    }
//}
//
//var Mesh::computeFaceHeight(const Eigen::RowVector3i& face, const int corner) {
//    Dual3DVector v0 = V.row(face((corner + 1) % 3));
//    Dual3DVector v1 = V.row(face((corner + 2) % 3));
//    Dual3DVector v2 = V.row(face(corner));
//
//    // Base of the triangle not including the corner
//    Dual3DVector edge = v1 - v0;
//
//    // Compute the height using the triangle area formula:
//    // Area = 0.5 * |edge x (v2 - v0)|, Height = 2 * Area / |edge|
//    var area = 0.5 * edge.cross(v2 - v0).norm();
//    var norm = edge.norm();
//    if (norm < Epsilon) return Epsilon;
//    else return (2.0 * area) / norm;
//}
//
//void Mesh::computeEdgeNorms(int i, DualVector& norms) {
//    std::vector<int> incEdges = VE[i];
//    norms.setZero(incEdges.size());
//    int ni = 0;
//    for (int edgeI : incEdges) {
//        int v0 = uE(edgeI, 0);
//        int v1 = uE(edgeI, 1);
//
//        Dual3DVector edge = V.row(v1) - V.row(v0);
//        norms(ni) = edge.norm() * stiffness(edgeI);
//        ni++;
//    }
//}


// For all vertices

//void Mesh::computeAverageHeights(DualVector& heights) {
//    // Initialize to the right size
//    heights.resize(uE.rows());
//
//    for (int i = 0; i < EF.rows(); i++) {
//        if (EF(i, 1) == -1 || EF(i, 0) == -1) { //boundary edge
//            heights(i) = 0.0f;
//            continue;
//        }
//        // Get the indices of the adjacent faces
//        int face0 = EF(i, 0);
//        int face1 = EF(i, 1);
//
//        // Get the coner point adjacent to the edge
//        int corner0 = EI(i, 0);
//        int corner1 = EI(i, 1);
//
//        var height0 = computeFaceHeight(F.row(face0), corner0);
//        var height1 = computeFaceHeight(F.row(face1), corner1);
//
//        // Only compute a third of the average
//        heights(i) = (height0 + height1) / 6.f;
//    }
//}


//void Mesh::computeEdgeNorms(DualVector& norms) {
//    for (int i = 0; i < uE.rows(); i++) {
//        int v0 = uE(i, 0);
//        int v1 = uE(i, 1);
//
//        Dual3DVector edge = V.row(v1) - V.row(v0);
//
//        norms(i) = edge.norm();
//    }
//}






