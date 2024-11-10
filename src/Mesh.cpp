#include "Mesh.h"

// Constructor for Edge structure
Edge::Edge(int vertex1, int vertex2, double stiff, double angle)
    : v1(vertex1), v2(vertex2), stiffness(stiff), dihedralAngle(angle) {}

// Constructor for Mesh class
Mesh::Mesh(const Eigen::MatrixXd& v, const Eigen::MatrixXi& f, const Eigen::MatrixXi& e)
    : vertices(v), faces(f) {
    computeFaceNormals();  // Compute the normals for each face
    buildEdges(e);         // Build edge data structures
    buildAdjacencyMatrix(); // Build adjacency matrix for face connectivity
}

// Computes the normals for each face in the mesh
void Mesh::computeFaceNormals() {
    faceNormals.resize(faces.rows());
    Eigen::Vector3d referenceNormal(0, 0, 1);  // Use the global up direction as reference

    for (int i = 0; i < faces.rows(); ++i) {
        Eigen::Vector3d v0 = vertices.row(faces(i, 0));
        Eigen::Vector3d v1 = vertices.row(faces(i, 1));
        Eigen::Vector3d v2 = vertices.row(faces(i, 2));

        // Calculate two edges of the face
        Eigen::Vector3d edge1 = v1 - v0;
        Eigen::Vector3d edge2 = v2 - v0;
        // Cross product of edges gives the normal
        Eigen::Vector3d normal = edge1.cross(edge2);

        // Ensure normal points in the correct direction by comparing with reference
        if (normal.dot(referenceNormal) < 0) {
            std::swap(faces(i, 1), faces(i, 2));  // Correct vertex order if normal is inverted
            normal = -normal;  // Update the normal
        }

        faceNormals[i] = normal.normalized();  // Normalize and store the normal
    }
}

// Build the edges from a given matrix of edge indices
void Mesh::buildEdges(const Eigen::MatrixXi& edges) {
    edgeList.reserve(edges.rows());  // Pre-allocate memory for efficiency
    for (int i = 0; i < edges.rows(); ++i) {
        int v1 = edges(i, 0);
        int v2 = edges(i, 1);
        edgeList.emplace_back(v1, v2);
    }

    // Populate adjacentFaces for each edge by mapping edges to their indices
    std::unordered_map<std::pair<int, int>, int> edgeMap;
    for (int i = 0; i < edgeList.size(); ++i) {
        edgeMap[{edgeList[i].v1, edgeList[i].v2}] = i;
        edgeMap[{edgeList[i].v2, edgeList[i].v1}] = i;
    }

    // Assign adjacent faces to each edge
    for (int i = 0; i < faces.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = faces(i, j);
            int v2 = faces(i, (j + 1) % 3);
            int edgeIndex = edgeMap[{v1, v2}];
            edgeList[edgeIndex].adjacentFaces.push_back(i);
        }
    }

    // Calculate the dihedral angles for each edge using the populated adjacent faces
    for (Edge& edge : edgeList) {
        edge.dihedralAngle = calculateDihedralAngle(edge, faceNormals);
        if (fabs(edge.dihedralAngle) > 0.0) { // Apply a threshold to determine if the edge is a crease
            edge.stiffness = 0.5;
        }
    }
}

// Calculates the dihedral angle between the two faces adjacent to an edge
double Mesh::calculateDihedralAngle(Edge& edge, const std::vector<Eigen::Vector3d>& faceNormals) {
    if (edge.adjacentFaces.size() < 2) return 0.0;  // Return 0 for boundary edges

    const Eigen::Vector3d& n1 = faceNormals[edge.adjacentFaces[0]];
    const Eigen::Vector3d& n2 = faceNormals[edge.adjacentFaces[1]];

    double cosAngle = n1.dot(n2);
    return acos(std::clamp(cosAngle, -1.0, 1.0));  // Return the angle in radians
}

// Builds a sparse matrix representing face adjacency
void Mesh::buildAdjacencyMatrix() {
    igl::facet_adjacency_matrix(faces, adjacencyMatrix);
}
