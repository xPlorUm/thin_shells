#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>
#include <igl/readTGF.h>
#include <igl/readDMAT.h>
#include <igl/forward_kinematics.h>
#include <igl/directed_edge_parents.h>
#include <igl/rotate_by_quat.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/cat.h>
#include <igl/harmonic.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/grad.h>



#include "Lasso.h"
#include "Colors.h"

#include "Deformation.h"

//activate this for alternate UI (easier to debug but no interactive updates, turn this OFF for your report)
//#define UPDATE_ONLY_ON_UP

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// PATHS
std::string PATH = "C:\\Users\\cedri\\source\\repos\\thin_shells\\assignment6";
std::string rel_data_path = "\\data\\hand\\";
std::string handles_file = "handles.dmat";
std::string original_pose_file = "pose.dmat";
std::string target_pose_file = "pose_target.dmat";
std::string skeleton_file = "skeleton.tgf";
std::string obj_file = "rest.obj";


//Assignment 6
Eigen::MatrixXd orig_pose, target_pose;
Eigen::VectorXi handles;

//vertex array, #V x3
Eigen::MatrixXd V, V_original;
//face array, #F x3
Eigen::MatrixXi F;
//colors for mesh
Eigen::MatrixXd colors;
//Vertices of Cage
Eigen::MatrixXd C, orig_C;
//Edges of Cage
Eigen::MatrixXi BE;
//Parents of Cage
Eigen::VectorXi P;
//Harmonic skinning weights
Eigen::MatrixXd W, W_f;

//Rotations
typedef std::vector < Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> > R_list;
R_list target_pose_list; //target animation
std::vector<R_list> poses; //sequential animation
double anim_t = 0; // current step of animation
double step_size = 1;
int num_poses_target = 50; //how many steps are used for target animation
bool sequential_anim = true;
bool animation = true;
bool show_skeleton = false;
bool face_based_coloring = false; // coloring of mesh

enum SkinningMode {
    LBS, DQS, PFB
};
SkinningMode skinning_mode = LBS;



int num_rotations, num_frames, num_faces;

//handles
int num_vertices, num_bones;
Eigen::MatrixXd M; //LBS Matrix

double distanceEdgeToVertice(Eigen::VectorXd e0, Eigen::VectorXd e1, Eigen::VectorXd v);

void computeHarmonicWeights(Eigen::MatrixXd& W, Eigen::MatrixXd& W_f);
void setColumnsToZero(Eigen::SparseMatrix<double>& matrix, const std::vector<int>& columnsToZero);

bool callback_pre_draw(Viewer &viewer);
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);


bool load_mesh(string path, string obj_file, string handles_file, string original_pose_file, string target_pose_file, string skeleton_file) {
    igl::read_triangle_mesh(path + obj_file, V, F);
    //Sequential animation
    igl::readDMAT(path + original_pose_file, orig_pose);
    //Target animation
    igl::readDMAT(path + target_pose_file, target_pose);
    //Handles
    igl::readDMAT(path + handles_file, handles);
    igl::readTGF(path + skeleton_file, C, BE);
    orig_C = C;
    V_original = Eigen::MatrixXd(V);

    

    //geometric handle selection
    //assign the closest vertice to each bone
    num_vertices = V.rows();
    num_bones = BE.rows();
    handles = Eigen::VectorXi::Constant(num_vertices, -1);
    Eigen::VectorXd selected_dis(num_bones);
    Eigen::MatrixXd D(num_bones, num_vertices);
    for (int b = 0; b < num_bones; b++) {
        //compute min_distance
        double min_dis = 100000;
        for (int i = 0; i < num_vertices; i++) {
            Eigen::Vector3d e0 = C.row(BE(b, 0));
            Eigen::Vector3d e1 = C.row(BE(b, 1));
            Eigen::Vector3d v = V.row(i);
            double dis = distanceEdgeToVertice(e0, e1, v);
            D(b, i) = dis;
            if (dis < min_dis) {
                min_dis = dis;
            }
        }
        selected_dis(b) = min_dis;
    }

     //assign vertices to bones
    for (int b = 0; b < num_bones; b++) {
        for (int i = 0; i < num_vertices; i++) {
            if (D(b, i) < 1.5 * selected_dis(b)) {
                // if there are multiple assignments per vertex then take the closer bone
                if (handles(i) != -1) {
                    if (D(handles(i), i) > D(b, i)) {
                        handles(i) = b;
                    }
                }
                else {
                    handles(i) = b;
                }
            }
        }
    }


    //set the color of the vertices of the first 3 bones
    num_vertices = V.rows();
    num_bones = BE.rows();


    //set colors on face basis with W_f
    if (face_based_coloring) {
        colors = Eigen::MatrixXd::Ones(F.rows(), 3);
    }
    else {
        colors = Eigen::MatrixXd::Ones(V.rows(), 3);
    }
;

    // set the skeleton with blue edges and red points
    viewer.data().clear();
    if (show_skeleton) {
        viewer.data().set_points(C, Eigen::RowVector3d(255.0, 0.0, 0.0));
        viewer.data().set_edges(C, BE, Eigen::RowVector3d(0.0, 0.0, 255.0));
        viewer.data().show_overlay_depth = false;
        viewer.core().align_camera_center(C);
    }
    else {
        viewer.data().set_mesh(V, F);
        viewer.data().set_colors(colors);
        viewer.core().align_camera_center(V);
    }

    return true;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage assignment6 mesh.off>" << endl;
        load_mesh(PATH + rel_data_path, obj_file, handles_file, original_pose_file, target_pose_file, skeleton_file);
    }
    //assignment 6
    igl::directed_edge_parents(BE, P);


    num_rotations = BE.rows();
    num_faces = F.rows();


    //Initialize different Poses 
    //Sequential animation
    num_frames = orig_pose.rows() / num_rotations;
    poses = std::vector<R_list>(num_frames);
    for (int f = 0; f < num_frames; f++) {
        R_list pose = R_list(num_rotations);
        for (int r = 0; r < num_rotations; r++) {
            pose[r] = Eigen::Quaterniond(orig_pose(f*num_rotations + r, 3), orig_pose(f * num_rotations + r, 0), orig_pose(f * num_rotations + r, 1), orig_pose(f * num_rotations + r, 2));
        }
        poses[f] = pose;

    }

    //Initialize R_lists for target_animation
    target_pose_list = R_list(num_rotations);
    for (int r = 0; r < num_rotations; r++) {
        target_pose_list[r] = Eigen::Quaterniond(target_pose(r, 3), target_pose(r, 0), target_pose(r, 1), target_pose(r, 2));
    }



    //solve laplacian to get scalar fields
    computeHarmonicWeights(W, W_f);
    

    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);


    viewer.callback_key_down = callback_key_down;
    viewer.callback_pre_draw = callback_pre_draw;

    viewer.data().point_size = 10;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.launch();
}

bool callback_pre_draw(Viewer& viewer) {
    Eigen::MatrixXd CT(C);
    Eigen::MatrixXi BET(BE);
    Eigen::MatrixXd U(V);
    if (animation) {
        R_list vQ;
        std::vector<Vector3d> vT;

        //Find interval
        int num_poses;
        if (sequential_anim) {
            num_poses = poses.size();
        }
        else {
            num_poses = num_poses_target;
        }

        int begin_step = (int)(anim_t) % num_poses;
        int end_step = (int)(floor(anim_t) + 1) % num_poses;

        //Get relative Rotations in dQ
        R_list dQ(num_rotations);
        for (int i = 0; i < num_rotations; i++) {
            if (sequential_anim) {
                dQ[i] = poses[begin_step][i];
            }
            else {
                Eigen::Quaterniond q0 = Eigen::Quaterniond::Identity();
                double t = begin_step / (double)num_poses;
                dQ[i] = q0.slerp(t, target_pose_list[i]);

            }
        }

        //Propagate relative rotations to get absolute rotation by using FK
        igl::forward_kinematics(C, BE, P, dQ, vQ, vT);


        std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> vA;

        Eigen::MatrixXd C_res = Eigen::MatrixXd::Zero(C.rows(), C.cols());
        ////Compute Deformation
        for (int i = 0; i < num_rotations; i++) {
            // Affine Matrix
            Affine3d a = Affine3d::Identity();
            a.translate(vT[i]);
            a.rotate(vQ[i]);

            vA.push_back(a);

        }



        //deform skeleton
        igl::deform_skeleton(C, BE, vA, CT, BET);


        //deform mesh
        if (!show_skeleton) {

            if (skinning_mode == LBS)
            {

                std::vector<Eigen::MatrixXd> MeshesPerBone; //inefficient yet

                for (int i = 0; i < num_rotations; i++) {
                    // get Affine transformation
                    Affine3d a = Affine3d::Identity();
                    a.translate(vT[i]);
                    a.rotate(vQ[i]);
                    Eigen::MatrixXd A = a.matrix();

                    // transform Vertices
                    Eigen::MatrixXd Vk_hom(num_vertices, 4);
                    Vk_hom << V_original, Eigen::VectorXd::Ones(num_vertices, 1);
                    Eigen::MatrixXd Vk_transformed(num_vertices, 4);
                    Vk_transformed = (A * Vk_hom.transpose()).transpose();

                    //multiply with weights
                    Eigen::MatrixXd Vk(num_vertices, 3);
                    for (int j = 0; j < num_vertices; j++) {
                        Vk.row(j) = W(j, i) * Vk_transformed.row(j).head(3);
                    }
                    MeshesPerBone.push_back(Vk);
                }

                //get weighted sum
                U.setZero();
                for (int i = 0; i < num_bones; i++) {
                    U = U + MeshesPerBone[i];
                }
            }
            else if (skinning_mode == DQS) {
                //check for all Quaternions if they are on the same hemisphere as refQ
                Eigen::Quaterniond refQ = vQ[0];
                for (int i = 0; i < num_rotations; i++) {
                    if (vQ[i].dot(refQ) < 0.0) {
                        vQ[i].coeffs() = -vQ[i].coeffs();
                    }
                }

                igl::dqs(V, W, vQ, vT, U);

            }
            else if (skinning_mode == PFB) {
                R_list vQ_avg;


                //get Quaternions as Matrix
                Eigen::MatrixXd vQ_m(num_bones, 4);
                for (int i = 0; i < vQ.size(); i++) {
                    vQ_m.row(i) = vQ[i].coeffs().transpose();
                }

                for (int i = 0; i < num_faces; i++) {
                    Eigen::MatrixXd Q(num_bones, 4);
                    Eigen::VectorXd w = W_f.row(i);


                    for (int j = 0; j < num_bones; j++) {
                        Q.row(j) = w[j] * vQ_m.row(j);
                    }


                    // find max q by using Eigen_Decomposition of QT*Q
                    Eigen::MatrixXd QQ(4, 4);
                    QQ = Q.transpose() * Q;


                    Eigen::EigenSolver<Eigen::Matrix4d> es;
                    es.compute(QQ);
                    Eigen::Vector4d eigenvalues = es.eigenvalues().real();
                    Eigen::Matrix4d eigenvectors = es.eigenvectors().real();


                    int maxIndex;
                    eigenvalues.maxCoeff(&maxIndex);
                    Eigen::Vector4d q_res_vec = eigenvectors.col(maxIndex);
                    Eigen::Quaterniond q_res(q_res_vec(0), q_res_vec(1), q_res_vec(2), q_res_vec(3));
                    q_res = q_res.normalized();


                    vQ_avg.push_back(q_res);


                }


                //poisson stitching 
                
                //get handles for bone 0 as constraints
                std::vector<int> first_handle;
                for (int i = 0; i < num_vertices; i++) {
                    if (handles[i] == 0) {
                        first_handle.push_back(i);
                    }
                }

                //make LHS
                // LHS = [G; C]T * [G; C]
                Eigen::SparseMatrix<double> LHS(num_vertices, num_vertices);
                Eigen::SparseMatrix<double> A(3 * num_faces + first_handle.size(), num_vertices);
                Eigen::SparseMatrix<double> C(first_handle.size(), num_vertices);


                Eigen::SparseMatrix<double> G;
                igl::grad(V, F, G); // G has size (3*num_faces, num_vertices)


                //remove cols in G with fixed vertices
                setColumnsToZero(G, first_handle);
                

                //make RHS
                Eigen::MatrixXd RHS(num_vertices, 3);
                //set rotations
                Eigen::MatrixXd R(3 * num_faces + first_handle.size(), 3);
                for (int i = 0; i < num_faces; i++) {
                    R.block(i * 3, 0, 3, 3) = vQ_avg[i].toRotationMatrix().transpose();
                }


                // add constraints on both sides
                for (int i = 0; i < first_handle.size(); i++) {
                    C.insert(i, first_handle[i]) = 1.0;
                    R.row(3*num_faces + i) = V.row(first_handle[i]);
                }



                igl::cat(1, G, C, A);


                RHS = A.transpose() * R;

                LHS = A.transpose() * A;

                //use sparsecholesky only symmetric matrices
                Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;

                solver.compute(LHS);

                // Check if the decomposition was successful
                if (solver.info() != Eigen::Success) {
                    std::cerr << "Decomposition failed" << std::endl;
                    return -1;
                }


                U = solver.solve(RHS);

                // Check if the solution was successful
                if (solver.info() != Eigen::Success) {
                    std::cerr << "Solving failed" << std::endl;
                    return -1;
                }

                //check if first_handle bones stay the same
                Eigen::VectorXi first_handle_vec(first_handle.size());
                for (int i = 0; i < first_handle.size(); i++) {
                    first_handle_vec[i] = first_handle[i];
                }
                Eigen::MatrixXd U_handle, V_handle;
                igl::slice(U, first_handle_vec, U_handle);
                igl::slice(V, first_handle_vec, V_handle);


            }

        }

        anim_t += step_size;
    }
    //render mesh
    viewer.data().clear();
    if (show_skeleton) {
        viewer.data().set_points(CT, Eigen::RowVector3d(255.0, 0.0, 0.0));
        viewer.data().set_edges(CT, BET, Eigen::RowVector3d(0.0, 0.0, 255.0));
        viewer.data().show_overlay_depth = false;
        viewer.core().align_camera_center(C);
    }
    else {
        // assign colors
        if (face_based_coloring) {
            colors.col(0) = Eigen::VectorXd::Ones(F.rows()) - W_f.col(0);
            colors.col(1) = Eigen::VectorXd::Ones(F.rows()) - W_f.col(1);
            colors.col(2) = Eigen::VectorXd::Ones(F.rows()) - W_f.col(2);
        }
        else {
            colors.col(0) = Eigen::VectorXd::Ones(V.rows()) - W.col(0);
            colors.col(1) = Eigen::VectorXd::Ones(V.rows()) - W.col(1);
            colors.col(2) = Eigen::VectorXd::Ones(V.rows()) - W.col(2);

        }

        viewer.data().clear();
        viewer.data().set_mesh(U, F);
        viewer.data().set_colors(colors);
        viewer.core().align_camera_center(U);
    }

    return false;
}


// Function to set certain columns of a sparse matrix to zero
void setColumnsToZero(Eigen::SparseMatrix<double>& matrix,const std::vector<int>& columnsToZero) {
    for (int col : columnsToZero) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, col); it; ++it) {
            it.valueRef() = 0.0;
        }
        
    }
}

void computeHarmonicWeights(Eigen::MatrixXd &W, Eigen::MatrixXd& W_f) {
    //get bones per assigned vertex in handle_bones
    std::vector<int> handle_ids;
    std::vector<int> handle_bones;
    for (int i = 0; i < num_vertices; i++) {
        if (handles(i) != -1) {
            handle_ids.push_back(i);
            handle_bones.push_back(handles(i));
        }
    }

    int num_constraints = handle_ids.size();

    //make Constraints matrix
    Eigen::MatrixXd bc(num_constraints, num_bones); //#b x #W
    bc.setZero();

    //matrix bc contains at i,j = 1 only if handle_vertex_i is assigned to bone j, else the values are 0
    for (int i = 0; i < num_constraints; i++) {
        bc(i, handle_bones[i]) = 1.0;
    }

    //convert to Eigen
    Eigen::VectorXi b = Eigen::VectorXi::Map(handle_ids.data(), handle_ids.size());

    igl::harmonic(V, F, b, bc, 1, W); // W of dimension #V x #W

    //get weights for each face
    W_f = Eigen::MatrixXd(num_faces, W.cols());
    for (int i = 0; i < num_faces; i++) {
        W_f.row(i) = (W.row(F(i, 0)) + W.row(F(i, 1)) + W.row(F(i, 2))) / 3.0;
    }


}


// calculates squared distance from vertice v to edge (e0, e1)
double distanceEdgeToVertice(Eigen::VectorXd e0, Eigen::VectorXd e1, Eigen::VectorXd v) {
    Eigen::Vector3d edge = e1 - e0;
    Eigen::Vector3d e0_v = v - e0;

    double t = edge.dot(e0_v) / edge.squaredNorm(); //projection factor t

    // Clamp t to the segment [0, 1]
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;


    // Compute the projection point on the segment
    Eigen::Vector3d projection = e0 + t * edge;

    // Compute the squared distance between w and the projection
    return sqrt((v - projection).squaredNorm());

}


bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
    bool handled = false;
    if (key == 'Q') {
        skinning_mode = LBS;
        handled = true;
    }

    if (key == 'W') {
        skinning_mode = DQS;
        handled = true;
    }

    if (key == 'E') {
        skinning_mode = PFB;
        handled = true;
    }

    return handled;
}

