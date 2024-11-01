#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>

#include <memory>

#include "Lasso.h"
#include "Colors.h"

#include "Deformation.h"

//activate this for alternate UI (easier to debug but no interactive updates, turn this OFF for your report)
//#define UPDATE_ONLY_ON_UP

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

string PATH = "/home/hugues/Desktop/ETHZ/MA3/PBS/thin_shells/data/";

string file = "woody-hi.off";

//vertex array, #V x3
Eigen::MatrixXd V(0, 3), V_original(0, 3);
//face array, #F x3
Eigen::MatrixXi F(0, 3);

Deformation solution;

//mouse interaction
enum MouseMode {
    SELECT, TRANSLATE, ROTATE, NONE
};


//for selecting vertices
std::unique_ptr<Lasso> lasso;
//list of currently selected vertices
Eigen::VectorXi selected_v(0, 1);

//for saving constrained vertices
//vertex-to-handle index, #V x1 (-1 if vertex is free)
Eigen::VectorXi handle_id(0, 1);
//list of all vertices belonging to handles, #HV x1
Eigen::VectorXi handle_vertices(0, 1);
//centroids of handle regions, #H x1
Eigen::MatrixXd handle_centroids(0, 3);
//updated positions of handle vertices, #HV x3
Eigen::MatrixXd handle_vertex_positions(0, 3);
//index of handle being moved
int moving_handle = -1;
//rotation and translation for the handle being moved
Eigen::Vector3f translation(0, 0, 0);
Eigen::Vector4f rotation(0, 0, 0, 1.);
typedef Eigen::Triplet<double> T;
//per vertex color array, #V x3
Eigen::MatrixXd vertex_colors;

//When false, use standard displacement vectors for details, when true use Deformation Transfer from part 2
// if true, solve will be called in the next pre-draw call
bool needs_solve = false;


void compute_handle_centroids();


bool callback_pre_draw(Viewer &viewer);


bool load_mesh(string filename) {
    igl::read_triangle_mesh(filename, V, F);
    viewer.data().clear();
    viewer.data().set_mesh(V, F);

    viewer.core().align_camera_center(V);
    V_original = V;
    handle_id.setConstant(V.rows(), 1, -1);
    // Initialize selector
    lasso = std::make_unique<Lasso>(V, F, viewer);

    selected_v.resize(0, 1);

    solution.set_initial_mesh(V, F);

    return true;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage thin_shells mesh.off>" << endl;
        load_mesh(PATH + file);
    } else {
        load_mesh(argv[1]);
    }

    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);


    //define the User Interface
    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content
        // menu.draw_viewer_menu();
    };

    //define callback functions of keys
    viewer.callback_pre_draw = callback_pre_draw;

    viewer.data().point_size = 10;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.launch();
}

// TODO : not sure how useful that is, so I keep it here for now
bool callback_pre_draw(Viewer &viewer) {
#ifndef UPDATE_ONLY_ON_UP
    if (needs_solve) {
        solve(viewer);
    }
#endif
    // initialize vertex colors
    vertex_colors = Eigen::MatrixXd::Constant(V.rows(), 3, .9);

    // first, color constraints
    for (int i = 0; i < V.rows(); ++i)
        if (handle_id[i] != -1) {
            int r = handle_id[i] % MAXNUMREGIONS;
            vertex_colors.row(i) << regionColors[r][0], regionColors[r][1], regionColors[r][2];
        }
    // then, color selection
    for (int i = 0; i < selected_v.size(); ++i)
        vertex_colors.row(selected_v[i]) << 131. / 255, 131. / 255, 131. / 255.;

    viewer.data().set_colors(vertex_colors);

    //clear points and lines
    viewer.data().set_points(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));
    viewer.data().set_edges(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXi::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));

    //draw the stroke of the selection
    for (unsigned int i = 0; i < lasso->strokePoints.size(); ++i) {
        viewer.data().add_points(lasso->strokePoints[i], Eigen::RowVector3d(0.4, 0.4, 0.4));
        if (i > 1)
            viewer.data().add_edges(lasso->strokePoints[i - 1], lasso->strokePoints[i],
                                    Eigen::RowVector3d(0.7, 0.7, 0.7));
    }


    //update the vertex position all the time
    viewer.data().set_mesh(V, F);

#ifdef UPDATE_ONLY_ON_UP
    //draw only the moving parts with a white line
    if (moving_handle>=0)
    {
      Eigen::MatrixXd edges(3*F.rows(),6);
      int num_edges = 0;
      for (int fi = 0; fi<F.rows(); ++fi)
      {
        int firstPickedVertex = -1;
        for(int vi = 0; vi<3 ; ++vi)
          if (handle_id[F(fi,vi)] == moving_handle)
          {
            firstPickedVertex = vi;
            break;
          }
        if(firstPickedVertex==-1)
          continue;


        Eigen::Matrix3d points;
        for(int vi = 0; vi<3; ++vi)
        {
          int vertex_id = F(fi,vi);
          if (handle_id[vertex_id] == moving_handle)
          {
            int index = -1;
            // if face is already constrained, find index in the constraints
            (handle_vertices.array()-vertex_id).cwiseAbs().minCoeff(&index);
            points.row(vi) = handle_vertex_positions.row(index);
          }
          else
            points.row(vi) =  V.row(vertex_id);

        }
        edges.row(num_edges++) << points.row(0), points.row(1);
        edges.row(num_edges++) << points.row(1), points.row(2);
        edges.row(num_edges++) << points.row(2), points.row(0);
      }
      edges.conservativeResize(num_edges, Eigen::NoChange);
      viewer.data().add_edges(edges.leftCols(3), edges.rightCols(3), Eigen::RowVector3d(0.9,0.9,0.9));
    }
#endif
    return false;

}

