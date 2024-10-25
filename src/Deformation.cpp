#include "Deformation.h"




void Deformation::update_handle_vertex_selection(const Eigen::VectorXi &new_handle_id,
                                                 const Eigen::VectorXi &new_handle_vertices)
{
    // handle_id: the vertex-to-handle index, #V x1 (-1 if vertex is free)
    // handle_vertices: list of all vertices belonging to handles, #HV x1
    handle_id = new_handle_id;
    handle_vertices = new_handle_vertices;
    int N = V_original.rows();
    int numFree = N - handle_vertices.rows();
    int numHandle = handle_vertices.rows();

    free_vertices = Eigen::VectorXi(numFree);
    int count = 0;
    for (int i = 0; i < N; i++) {
        if (handle_id[i] == -1) {
            free_vertices(count) = i;
            count++;
        }
    }

    Eigen::SparseMatrix<double> L, M, A;
    igl::MassMatrixType type = igl::MASSMATRIX_TYPE_VORONOI;
    igl::cotmatrix(V_original, F, L);
    igl::massmatrix(V_original, F, type, M);
    A = L * (M.cwiseInverse()) * L;

    igl::slice(A, free_vertices.array(), free_vertices.array(), Aff);
    igl::slice(A, free_vertices.array(), handle_vertices.array(), Afc);

    solver.compute(Aff);

    if (solver.info() != Eigen::Success) {
        // decomposition failed
        std::cerr << "Decomposition failed!" << std::endl;
    }

}

void Deformation::get_smooth_mesh(Eigen::MatrixXd &V_res) {
    // Get the smooth mesh B
    // Store the result to V_res

    int N = V_original.rows();
    Eigen::MatrixXd Vc;
    Eigen::VectorXi a(3);
    a << 0, 1, 2;
    igl::slice(V_original, handle_vertices.array(), a.array(), Vc);


    Eigen::MatrixXd b = -1 * Afc * Vc;
    Eigen::MatrixXd Vf = solver.solve(b);


    V_res = Eigen::MatrixXd(V_original);
    int countFree = 0;
    int countHandle = 0;
    for (int i = 0; i < N; i++) {
        if (handle_id[i] == -1) {
            V_res.row(i) = Vf.row(countFree);
            countFree++;
        }
        else {
            V_res.row(i) = Vc.row(countHandle);
            countHandle++;
        }
    }

    if (solver.info() != Eigen::Success) {
        // solving failed
        std::cerr << "Solving failed!" << std::endl;
    }


}

void Deformation::get_deformed_smooth_mesh(const Eigen::MatrixXd &handle_vertex_positions, Eigen::MatrixXd &V_res) {
    // Given the handle vertex positions, get the deformed smooth mesh B'
    // Store the result to V_res

    int n = V_original.rows();
    Eigen::MatrixXd Vc = handle_vertex_positions;


    Eigen::MatrixXd b = -1 * Afc * Vc;
    Eigen::MatrixXd Vf = solver.solve(b);


    V_res = Eigen::MatrixXd(V_original);
    int countFree = 0;
    int countHandle = 0;
    for (int i = 0; i < n; i++) {
        if (handle_id[i] == -1) {
            V_res.row(i) = Vf.row(countFree);
            countFree++;
        }
        else {
            V_res.row(i) = Vc.row(countHandle);
            countHandle++;
        }
    }

    if (solver.info() != Eigen::Success) {
        // solving failed
        std::cerr << "Solving failed!" << std::endl;
    }




}
// TODO change back 
void Deformation::get_deformed_mesh(const Eigen::MatrixXd& handle_vertex_positions, Eigen::MatrixXd& V_res) {
    //void Deformation::get_deformed_mesh(const Eigen::MatrixXd &handle_vertex_positions, Eigen::MatrixXd &V_res) {
    // Given the handle vertex positions, get the deformed mesh with details S'
    // Store the result to V_res

    int n = V_original.rows();
    Eigen::MatrixXd B1;
    get_smooth_mesh(B1);

    Eigen::MatrixXd N1;
    igl::per_vertex_normals(B1, F, N1);

    // project all outgoing edges to the per vertex normal
    std::vector<std::vector<int>> AA;
    igl::adjacency_list(F, AA);

    Eigen::MatrixXd Displ = V_original - B1;
    Eigen::MatrixXd D(B1);
    Eigen::MatrixXd Displ_proj_B1(B1);
    Eigen::VectorXi outgoing_e(n);




    for (int i = 0; i < n; i++) {
        Eigen::Vector3d n0 = N1.row(i).normalized();
        Eigen::Vector3d v0 = B1.row(i);
        Eigen::Vector3d e0, e1, e2;
        e2 = n0; 

        // project all neighbors
        // find it's longest edge
        double max_dis = -1; int max_nb = 0;
        for (int j = 0; j < AA[i].size(); j++) {
            int nb = AA[i][j];
            Eigen::Vector3d nb_vec = B1.row(nb);

            Eigen::Vector3d edge = (nb_vec - v0);
            Eigen::Vector3d p_edge = (edge - edge.dot(n0) * n0);
            double dis = p_edge.norm();
            //std::cout <<nb << " has dis" << dis << std::endl;
            if (max_dis < dis) {
                max_dis = dis;
                max_nb = nb;
                e0 = p_edge.normalized();
                //std::cout << "chosen " << nb << std::endl;
                outgoing_e[i] = nb;

            }

        }

        e1 = e2.cross(e0).normalized();

        D(i, 0) = e0.dot(Displ.row(i));
        D(i, 1) = e1.dot(Displ.row(i));
        D(i, 2) = e2.dot(Displ.row(i));

        Displ_proj_B1.row(i) = D(i, 0) * e0 + D(i, 1) * e1 + D(i, 2) * e2;


    }





    Eigen::MatrixXd B2;
    get_deformed_smooth_mesh(handle_vertex_positions, B2);
    Eigen::MatrixXd N2;
    igl::per_vertex_normals(B2, F, N2);
    Eigen::MatrixXd Displ_proj_B2(B2);

    // show projected displacement vectors
    for (int i = 0; i < n; i++) {
        Eigen::Vector3d v0 = B2.row(i).normalized();
        Eigen::Vector3d n0 = N2.row(i).normalized();
        int nb = outgoing_e[i];
        Eigen::Vector3d nb_vec = B2.row(nb);


        Eigen::Vector3d e0, e1, e2;
        e2 = n0;

        Eigen::Vector3d edge = (nb_vec - v0).transpose();
        Eigen::Vector3d proj_edge = edge - edge.dot(n0) * n0;

        e0 = proj_edge.normalized();
        e1 = e2.cross(e0);


        D(i, 0) = e0.dot(Displ_proj_B1.row(i));
        D(i, 1) = e1.dot(Displ_proj_B1.row(i));
        D(i, 2) = e2.dot(Displ_proj_B1.row(i));


        // project Displacements
        Displ_proj_B2.row(i) = D(i, 0) * e0 + D(i, 1) * e1 + D(i, 2) * e2;
        //Displ_proj_B2.row(i) = Displ_proj_B1(i, 0) * e0 + Displ_proj_B1(i, 1) * e1 + Displ_proj_B1(i, 2) * e2;
    }

    V_res = B2 + Displ_proj_B2;


}

void Deformation::get_deformed_mesh_deformation_transfer(const Eigen::MatrixXd &handle_vertex_positions, Eigen::MatrixXd &V_res)
{
    // Implement deformation transfer here.
    // Store the result to V_res

    int n = V_original.rows();
    Eigen::MatrixXd B1;
    get_smooth_mesh(B1);

    Eigen::MatrixXd N1;
    igl::per_face_normals(B1, F, N1);

    Eigen::MatrixXd B2;
    get_deformed_smooth_mesh(handle_vertex_positions, B2);
    Eigen::MatrixXd N2;
    igl::per_face_normals(B2, F, N2);

    int m = F.rows();

    //get derivatives of faces
    Eigen::SparseMatrix<double> G(3 * m, n);
    //use L instead of grad??
    igl::grad(V_original, F, G);

    // get areas of triangles
    Eigen::SparseMatrix<double> D(3 * m, 3 * m);
    Eigen::VectorXd double_areas;
    igl::doublearea(V_original, F, double_areas);
    Eigen::VectorXd areas = double_areas / 2.0;

    for (int i = 0; i < m; i++) {
        D.coeffRef(3*i + 0, 3*i + 0) = areas[i];
        D.coeffRef(3*i + 1, 3*i + 1) = areas[i];
        D.coeffRef(3*i + 2, 3*i + 2) = areas[i];
    }


    Eigen::MatrixXd S(3 * m, 3);

    for (int i = 0; i < m; i++) {
        int i0 = F(i, 0);
        int i1 = F(i, 1);
        int i2 = F(i, 2);
        Eigen::Matrix3d Sj1, Sj2, Sj;
        Sj1 << (B1.row(i0) - B1.row(i2)).normalized(), (B1.row(i1) - B1.row(i2)).normalized(), N1.row(i).normalized();
        Sj2 << (B2.row(i0) - B2.row(i2)).normalized(), (B2.row(i1) - B2.row(i2)).normalized(), N2.row(i).normalized();
        Sj = Sj2 * Sj1.inverse();
        S.block(i * 3, 0, 3, 3) << Sj;
    }

    Eigen::SparseMatrix<double> A = G.transpose() * D * G;
    Eigen::MatrixXd b = G.transpose() * D * S;

    Deformationsolver.compute(A);
    if (Deformationsolver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
    }
    Eigen::MatrixXd Vf = Deformationsolver.solve(b);
    if (Deformationsolver.info() != Eigen::Success) {
        std::cerr << "Solving failed!" << std::endl;
    }
    V_res = Vf;

}
