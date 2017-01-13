// Tao Du
// taodu@csail.mit.edu
// Jan 9, 2017
#include <string>
#include "hex_mesh.h"
#include "quad_converter.h"

// pbrt_file lattice_file dispalcement_file material_file sing_point_file fine_intf_flag_file
int main(int argc, char* argv[]) {
  if (argc < 7) return 0;
  const std::string pbrt_file(argv[1]);
  const std::string lattice_file(argv[2]);
  const std::string displacement_file(argv[3]);
  const std::string material_file(argv[4]);
  const std::string sing_point_file(argv[5]);
  const std::string fine_intf_flag_file(argv[6]);
  HexMesh hex_mesh(lattice_file, displacement_file, material_file, sing_point_file, fine_intf_flag_file);
  hex_mesh.Normalize();
  hex_mesh.ToPBRT(pbrt_file);
  return 0;
}