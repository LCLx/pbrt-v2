// Tao Du
// taodu@csail.mit.edu
// Jan 9, 2017
#include <string>
#include "hex_mesh.h"
#include "quad_mesh.h"

// pbrt_file root_folder hex_mesh_name rho_mesh_name depth normalize
// pbrt_file root_folder t
// pbrt_file root_folder lattice_file displacement_file material_file lag_inf_point_file
// sing_point_file fine_intf_flag_file f_point_file psi_D_file density_file plot_surrounding_cells
// threshold normalize v0 v1 v2 v3
int main(int argc, char* argv[]) {
  if (argc == 6 || argc == 7) {
    int argument_index = 1;
    const std::string pbrt_file(argv[argument_index++]);
    const std::string root_folder(argv[argument_index++]);
    const std::string hex_mesh_name(argv[argument_index++]);
    const std::string rho_mesh_name(argv[argument_index++]);
    const double depth = atof(argv[argument_index++]);
    bool normalize = true;
    if (argument_index < argc) {
      normalize = atoi(argv[argument_index++]) != 0;
    }
    QuadMesh quad_mesh(root_folder + "\\" + hex_mesh_name,
      root_folder + "\\" + rho_mesh_name);
    if (normalize) quad_mesh.Normalize();
    quad_mesh.ToPBRT(pbrt_file, depth);
    return 0;
  } else if (argc == 4) {
    // Hard-code almost everything...
    int argument_index = 1;
    const std::string pbrt_file(argv[argument_index++]);
    const std::string root_folder(argv[argument_index++]);
    const double t = atof(argv[argument_index++]);
    HexMesh mesh0(root_folder + "\\0\\lattice",
      root_folder + "\\0\\ale_dis",
      root_folder + "\\0\\material", "NULL", "NULL",
      root_folder + "\\0\\fine_intf_flags");
    mesh0.Normalize();
    HexMesh mesh1(root_folder + "\\1\\lattice",
      root_folder + "\\1\\ale_dis",
      root_folder + "\\1\\material", "NULL", "NULL",
      root_folder + "\\0\\fine_intf_flags");
    mesh1.Normalize();
    (mesh0 * (1 - t) + mesh1 * t).ToPBRT(pbrt_file);
    return 0;
  } else {
    int argument_index = 1;
    const std::string pbrt_file(argv[argument_index++]);
    const std::string root_folder(argv[argument_index++]);

    std::string lattice_file(argv[argument_index++]);
    if (lattice_file != "NULL") lattice_file = root_folder + "\\" + lattice_file;

    std::string displacement_file(argv[argument_index++]);
    if (displacement_file != "NULL") displacement_file = root_folder + "\\" + displacement_file;

    std::string material_file(argv[argument_index++]);
    if (material_file != "NULL") material_file = root_folder + "\\" + material_file;

    std::string lag_inf_point_file(argv[argument_index++]);
    if (lag_inf_point_file != "NULL") lag_inf_point_file = root_folder + "\\" + lag_inf_point_file;

    std::string sing_point_file(argv[argument_index++]);
    if (sing_point_file != "NULL") sing_point_file = root_folder + "\\" + sing_point_file;

    std::string fine_intf_flag_file(argv[argument_index++]);
    if (fine_intf_flag_file != "NULL") fine_intf_flag_file = root_folder + "\\" + fine_intf_flag_file;

    std::string f_point_file = "NULL", psi_D_file = "NULL", density_file = "NULL";
    if (argument_index < argc) {
      f_point_file = std::string(argv[argument_index++]);
      if (f_point_file != "NULL") f_point_file = root_folder + "\\" + f_point_file;
    }
    if (argument_index < argc) {
      psi_D_file = std::string(argv[argument_index++]);
      if (psi_D_file != "NULL") psi_D_file = root_folder + "\\" + psi_D_file;
    }
    if (argument_index < argc) {
      density_file = std::string(argv[argument_index++]);
      if (density_file != "NULL") density_file = root_folder + "\\" + density_file;
    }
    bool plot_surrounding_cells = false;
    if (argument_index < argc) {
      plot_surrounding_cells = atoi(argv[argument_index++]) != 0;
    }
    double threshold = 0.5;
    if (argument_index < argc) {
      threshold = atof(argv[argument_index++]);
    }
    bool normalize = true;
    if (argument_index < argc) {
      normalize = atoi(argv[argument_index++]) != 0;
    }
    std::string v0_file = "NULL", v1_file = "NULL", v2_file = "NULL", v3_file = "NULL";
    if (argument_index < argc) {
      v0_file = std::string(argv[argument_index++]);
      if (v0_file != "NULL") v0_file = root_folder + "\\" + v0_file;
    }
    if (argument_index < argc) {
      v1_file = std::string(argv[argument_index++]);
      if (v1_file != "NULL") v1_file = root_folder + "\\" + v1_file;
    }
    if (argument_index < argc) {
      v2_file = std::string(argv[argument_index++]);
      if (v2_file != "NULL") v2_file = root_folder + "\\" + v2_file;
    }
    if (argument_index < argc) {
      v3_file = std::string(argv[argument_index++]);
      if (v3_file != "NULL") v3_file = root_folder + "\\" + v3_file;
    }

    HexMesh hex_mesh(lattice_file, displacement_file, material_file, lag_inf_point_file,
      sing_point_file, fine_intf_flag_file, f_point_file, psi_D_file, density_file,
      v0_file, v1_file, v2_file, v3_file);
    if (normalize) hex_mesh.Normalize();
    hex_mesh.ToPBRT(pbrt_file, plot_surrounding_cells, threshold);
    return 0;
  }
}