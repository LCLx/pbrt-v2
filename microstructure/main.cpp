// Tao Du
// taodu@csail.mit.edu
// Jan 9, 2017
#include <string>
#include "hex_mesh.h"
#include "quad_converter.h"

// pbrt_file root_folder t
// pbrt_file root_folder lattice_file displacement_file material_file lag_inf_point_file sing_point_file fine_intf_flag_file
int main(int argc, char* argv[]) {
  if (argc == 4) {
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
  } else if (argc == 9) {
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

    HexMesh hex_mesh(lattice_file, displacement_file, material_file, lag_inf_point_file,
      sing_point_file, fine_intf_flag_file);
    hex_mesh.Normalize();
    hex_mesh.ToPBRT(pbrt_file);
    return 0;
  }
}