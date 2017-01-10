// Tao Du
// taodu@csail.mit.edu
// Jan 9, 2017
#include <string>
#include "hex_converter.h"
#include "quad_converter.h"

// -q quad_file.txt rho_file.txt
// -h hex_file.txt rho_file.txt
int main(int argc, char* argv[]) {
  if (argc < 4) return 0;
  const std::string mesh_file(argv[2]);
  const std::string rho_file(argv[3]);
  std::string pbrt_file;
  if (argc == 4) {
    pbrt_file = mesh_file.substr(0, mesh_file.rfind(".")) + ".pbrt";
  } else {
    pbrt_file = std::string(argv[4]);
  }
  if (argv[1][1] == 'q')
    ConvertQuadToPbrt(mesh_file, rho_file, pbrt_file);
  else
    ConvertHexToPbrt(mesh_file, rho_file, pbrt_file);
  return 0;
}