#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/rotate_by_quat.h>
#include "DiscreteShell.h"
#include <iostream>


#include <memory>
#include <utility>

//activate this for alternate UI (easier to debug but no interactive updates, turn this OFF for your report)
//#define UPDATE_ONLY_ON_UP

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

string PATH = "/home/hugues/Desktop/ETHZ/MA3/PBS/thin_shells/data/";
// string PATH = "C:\\Users\\cedri\\source\\repos\\thin_shells\\data\\";

string file = "old_assignment/woody-hi.off";

//vertex array, #V x3
Eigen::MatrixXd V(0, 3), V_original(0, 3); //V: Vertices which are to be changed temporarily during animation
//face array, #F x3
Eigen::MatrixXi F(0, 3);

//mouse interaction
enum MouseMode {
    SELECT, TRANSLATE, ROTATE, NONE
};


//for selecting vertices
//list of currently selected vertices
Eigen::VectorXi selected_v(0, 1);

//for saving constrained vertices
//vertex-to-handle index, #V x1 (-1 if vertex is free)
Eigen::VectorXi handle_id(0, 1);
//rotation and translation for the handle being moved
Eigen::Vector3f translation(0, 0, 0);
Eigen::Vector4f rotation(0, 0, 0, 1.);
typedef Eigen::Triplet<double> T;


// animation variables
double anim_t = 0; // current step of animation
double step_size = 1; 
int total_steps = 50; // how many steps are used for animation
// False by default, need to start that in the GUI
bool animation = false;


// Discrete Shell
DiscreteShell ds;

// functions for libigl
bool callback_pre_draw(Viewer& viewer);


bool load_mesh(string filename) {
    igl::read_triangle_mesh(std::move(filename), V, F);
    std::cout << "v zero  : " << V.isZero() << std::endl;
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.core().align_camera_center(V);
    V_original = V;
    handle_id.setConstant(V.rows(), 1, -1);
    selected_v.resize(0, 1);
    return true;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Usage : thin_shells <path-to-mesh>.off>" << endl;
        return 1;
    }
    ds.initializeFromFile(argv[1]);
    // Initialize Viewer
    viewer.data().clear();
    auto V = ds.getPositions();
    auto F = ds.getFaces();
    viewer.data().set_mesh(*V, *F);
    viewer.core().align_camera_center(*V);

    viewer.core().is_animating = true; // Enables continuous rendering
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);


    //define the User Interface
    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content (optional, can be commented out)
        // menu.draw_viewer_menu();

        // Add a button to start/stop the animation
        if (ImGui::Button(animation ? "Stop Animation" : "Start Animation")) {
            animation = !animation; // Toggle animation state
        }
        // Optionally, add a slider to control animation speed
        ImGui::SliderInt("Total Steps", &total_steps, 10, 200);
    };
    viewer.callback_pre_draw = callback_pre_draw;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.launch();
}

// callback_pre_draw runs continuously as long as the mouse hovers over the window
bool callback_pre_draw(Viewer& viewer) {
    //clear points and lines
    viewer.data().set_points(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));
    viewer.data().set_edges(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXi::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));

    //continue animation
    if (animation) {
        int begin_step = (int)(anim_t) % total_steps; // starts again whenever anim_t > num_poses
        int end_step = (int)(floor(anim_t) + 1) % total_steps;
        //TODO advance time step of discrete shell
        ds.advanceOneStep(begin_step);
        //TODO change/add function in discrete shell class which returns the resulting V:Vertices Matrix
    }
    //add time_step
    anim_t += step_size;
    // set value of Mesh to V
    viewer.data().clear();
    // Draw the discrete shell
    auto V = ds.getPositions();
    auto F = ds.getFaces();
    viewer.data().set_mesh(*V, *F);
    return false;
}
