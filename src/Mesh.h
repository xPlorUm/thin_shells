#ifndef MESH_H
#define MESH_H

#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <igl/facet_adjacency_matrix.h>
#include <cmath>

// Define a structure for an edge in the mesh
struct Edge {
    int v1, v2;  // Vertex indices of the edge's endpoints
    double stiffness;  // Stiffness coefficient used for simulations or rendering effects
    double dihedralAngle;  // Dihedral angle between the faces adjacent to the edge
    std::vector<int> adjacentFaces;  // Indices of faces adjacent to this edge

    // Constructor for initializing an edge with default or specified values
    Edge(int vertex1, int vertex2, double stiff = 1.0, double angle = 0.0);
};

// Mesh class representing a 3D mesh structure
class Mesh {
public:
    Eigen::MatrixXd vertices;  // Matrix storing vertex positions
    Eigen::MatrixXi faces;     // Matrix storing indices of vertices forming each face
    std::vector<Edge> edgeList;  // List of edges in the mesh
    Eigen::SparseMatrix<int> adjacencyMatrix;  // Adjacency matrix of faces
    std::vector<Eigen::Vector3d> faceNormals;  // Normals of each face

    // Constructor to initialize the mesh with vertices, faces, and edges
    Mesh(const Eigen::MatrixXd& v, const Eigen::MatrixXi& f, const Eigen::MatrixXi& e);

    // Computes the normals for each face in the mesh
    void computeFaceNormals();

    // Build the edges from a given matrix of edge indices
    void buildEdges(const Eigen::MatrixXi& edges);

    // Calculates the dihedral angle between the two faces adjacent to an edge
    double calculateDihedralAngle(Edge& edge, const std::vector<Eigen::Vector3d>& faceNormals);

    // Builds a sparse matrix representing face adjacency
    void buildAdjacencyMatrix();
};

#endif // MESH_H
