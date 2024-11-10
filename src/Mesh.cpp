#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <igl/facet_adjacency_matrix.h>

// Define a structure for an edge in the mesh
struct Edge {
    int v1, v2;  // Vertex indices of the edge's endpoints
    double stiffness;  // Stiffness coefficient used for simulations or rendering effects
    double dihedralAngle;  // Dihedral angle between the faces adjacent to the edge
    std::vector<int> adjacentFaces;  // Indices of faces adjacent to this edge

    // Constructor for initializing an edge with default or specified values
    Edge(int vertex1, int vertex2, double stiff = 1.0, double angle = 0.0)
        : v1(vertex1), v2(vertex2), stiffness(stiff), dihedralAngle(angle) {}
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
    Mesh(const Eigen::MatrixXd& v, const Eigen::MatrixXi& f, const Eigen::MatrixXi& e)
        : vertices(v), faces(f) {
        computeFaceNormals();  // Compute the normals for each face
        buildEdges(e);         // Build edge data structures
        buildAdjacencyMatrix(); // Build adjacency matrix for face connectivity
    }

    void computeFaceNormals() {
        faceNormals.resize(faces.rows());

        // Step 1: Compute the centroid of the entire mesh
        Eigen::Vector3d meshCentroid = vertices.colwise().mean();

        for (int i = 0; i < faces.rows(); ++i) {
            // Get vertices of the face
            Eigen::Vector3d v0 = vertices.row(faces(i, 0));
            Eigen::Vector3d v1 = vertices.row(faces(i, 1));
            Eigen::Vector3d v2 = vertices.row(faces(i, 2));

            // Step 2: Calculate two edges of the face
            Eigen::Vector3d edge1 = v1 - v0;
            Eigen::Vector3d edge2 = v2 - v0;
            
            // Cross product of edges gives the normal
            Eigen::Vector3d normal = edge1.cross(edge2).normalized();

            // Step 3: Calculate the centroid of the face
            Eigen::Vector3d faceCentroid = (v0 + v1 + v2) / 3.0;

            // Step 4: Check the alignment of the normal with respect to the mesh centroid
            Eigen::Vector3d directionToFace = faceCentroid - meshCentroid;
            if (normal.dot(directionToFace) < 0) {
                // If the normal is pointing inward, flip it by swapping vertices and inverting the normal
                std::swap(faces(i, 1), faces(i, 2));  // Reverse winding order
                normal = -normal;  // Flip normal direction
            }

            // Store the normalized normal
            faceNormals[i] = normal;
        }
    }


    // Build the edges from a given matrix of edge indices
    void buildEdges(const Eigen::MatrixXi& edges) {
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
    double calculateDihedralAngle(Edge& edge, const std::vector<Eigen::Vector3d>& faceNormals) {
        if (edge.adjacentFaces.size() < 2) return 0.0;  // Return 0 for boundary edges

        const Eigen::Vector3d& n1 = faceNormals[edge.adjacentFaces[0]];
        const Eigen::Vector3d& n2 = faceNormals[edge.adjacentFaces[1]];

        double cosAngle = n1.dot(n2);
        return acos(std::clamp(cosAngle, -1.0, 1.0));  // Return the angle in radians
    }

    // Builds a sparse matrix representing face adjacency
    void buildAdjacencyMatrix() {
        igl::facet_adjacency_matrix(faces, adjacencyMatrix);
    }
};
