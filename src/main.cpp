#include <igl/read_triangle_mesh.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/rotate_by_quat.h>
#include <igl/unproject_onto_mesh.h>
#include "DiscreteShell.h"
#include <iostream>
#include <filesystem>

#include <memory>
#include <utility>
#include "constants.h"


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"


#include <string.h>

std::string PATH = "../data/";
std::string file = "rectangle_upsampled_4.off";

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;
//vertex array, #V x3
Eigen::MatrixXd V(0, 3); //V: Vertices which are to be changed temporarily during animation
//face array, #F x3
Eigen::MatrixXi F(0, 3);

//mouse interaction
enum MouseMode {
    SELECT, TRANSLATE, ROTATE, NONE
};

int frame_count;


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
bool callback_pre_draw(Viewer &viewer);

bool callback_mouse_down(Viewer &viewer, int button, int modifier);

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers);


bool load_mesh(string filename) {
    igl::read_triangle_mesh(std::move(filename), V, F);
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.core().align_camera_center(V);
    handle_id.setConstant(V.rows(), 1, -1);
    selected_v.resize(0, 1);
    return true;
}


int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage : thin_shells <path-to-mesh>.off>" << endl;
        ds.initializeFromFile(PATH + file);
    } else {
        ds.initializeFromFile(argv[1]);
    }
#if SAVE_FRAMES == 1
    // Create the frames directory if it doesn't exist
    std::string frames_dir = "frames"; // Relative path; adjust as needed
    if (!std::filesystem::exists(frames_dir)) {
        if (!std::filesystem::create_directory(frames_dir)) {
            cerr << "Error: Could not create frames directory!" << endl;
            return EXIT_FAILURE;
        }
    }
#endif


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
    anim_t += step_size;


    //define the User Interface
    menu.callback_draw_viewer_menu = [&]() {

        // Add a button to start/stop the animation
        if (ImGui::Button(animation ? "Stop Animation" : "Start Animation")) {
            animation = !animation; // Toggle animation state
        }
        // Optionally, add a slider to control animation speed
        ImGui::SliderInt("Total Steps", &total_steps, 10, 200);
    };
    viewer.callback_pre_draw = callback_pre_draw;
    viewer.callback_mouse_down = callback_mouse_down;
    viewer.callback_key_down = callback_key_down;

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.launch();
}

// callback_pre_draw runs continuously as long as the mouse hovers over the window
bool callback_pre_draw(Viewer &viewer) {
    //clear points and lines
    viewer.data().set_points(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));
    viewer.data().set_edges(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXi::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));

    //continue animation
    if (animation) {
        int begin_step = (int) (anim_t) % total_steps; // starts again whenever anim_t > num_poses
        int end_step = (int) (floor(anim_t) + 1) % total_steps;
        ds.advanceOneStep(begin_step);
        anim_t += step_size;
    }
    viewer.data().clear();
    // Draw the discrete shell
    auto V = ds.getPositions();
    auto F = ds.getFaces();
    viewer.data().set_mesh(*V, *F);
    viewer.core().align_camera_center(*V);

    // Save frame if requested
#if SAVE_FRAMES == 1
    if (!animation)
        return false;
// Allocate temporary buffers for 1280x800 image
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

// Draw the scene in the buffers
    viewer.core().draw_buffer(viewer.data(),false,R,G,B,A);

// Save it to a PNG
    int height = R.rows();
    int width  = R.cols();

    std::vector<unsigned char> rgba(width * height * 4);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int idx = (y * width + x) * 4;
            rgba[idx + 0] = R(y, x); // Red channel
            rgba[idx + 1] = G(y, x); // Green channel
            rgba[idx + 2] = B(y, x); // Blue channel
            rgba[idx + 3] = A(y, x); // Alpha channel
        }
    }

    char filename[256];
    sprintf(filename, "frames/frame_%04d.png", frame_count++);

// The last parameter (stride_in_bytes) is width * number_of_channels
// Here we have 4 channels (RGBA).
    if (!stbi_write_png(filename, width, height, 4, rgba.data(), width*4)) {
        std::cerr << "Error: Failed to write " << filename << std::endl;
    } else {
        std::cout << "Saved frame " << filename << std::endl;
    }
#endif
    return false;
}

// Display debug information on debug
bool callback_mouse_down(Viewer &viewer, int button, int modifier) {
    int vi = -1;
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view /* viewer.data().model*/,
                                 viewer.core().proj, viewer.core().viewport, *ds.getPositions(), *ds.getFaces(), fid,
                                 bc)) {
        // paint hit red
        bc.maxCoeff(&vi);

        vi = ds.getFaces()->row(fid)(vi);
    }
    if (vi == -1) {
        return false;
    }
    ds.debug_vertex(vi);
    return false;
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
    // If space, stop the animation
    if (key == ' ') {
        animation = !animation;
    }
    return false;
}
