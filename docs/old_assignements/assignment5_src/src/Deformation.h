#ifndef ex6_Solution_h
#define ex6_Solution_h

#include <Eigen/Core>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <map>
#include <igl/slice.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/grad.h>
#include <igl/doublearea.h>


class Deformation
{
public:
    // DO NOT change the signature of any public members.
    // DO NOT add/remove any new public members.
    Eigen::MatrixXi F;   // Faces of the original mesh
    void set_initial_mesh(const Eigen::MatrixXd& V_, const Eigen::MatrixXi& F_) {
        V_original = V_;
        F = F_;
    }
    void update_handle_vertex_selection(const Eigen::VectorXi&, const Eigen::VectorXi&);
    void get_smooth_mesh(Eigen::MatrixXd&);
    void get_deformed_smooth_mesh(const Eigen::MatrixXd&, Eigen::MatrixXd&);
    // TODO change back 
    //void get_deformed_mesh(const Eigen::MatrixXd&, Eigen::MatrixXd&);
    void get_deformed_mesh(const Eigen::MatrixXd&, Eigen::MatrixXd&);
    void get_deformed_mesh_deformation_transfer(const Eigen::MatrixXd&, Eigen::MatrixXd&);

private:
    Eigen::MatrixXd V_original;  // Vertices of the original mesh
    // handle_id : V x1 Vector with entries -1: Free or i = indexOfHandle for ith handle
    Eigen::VectorXi handle_id;
    // handle_vertices : HVx1 = number of vertices which are part of a handle times vertex index <= V.rows()
    Eigen::VectorXi handle_vertices;
    Eigen::VectorXi free_vertices;

    Eigen::SparseMatrix<double> Aff, Afc;

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> solver;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> Deformationsolver;
};

#endif
