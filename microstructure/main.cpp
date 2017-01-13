// Tao Du
// taodu@csail.mit.edu
// Jan 9, 2017
#include <string>
#include "hex_mesh.h"
#include "quad_converter.h"

// pbrt_file root_folder lattice_file dispalcement_file material_file lag_inf_point_file sing_point_file fine_intf_flag_file
int main(int argc, char* argv[]) {
  if (argc < 9) return 0;
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