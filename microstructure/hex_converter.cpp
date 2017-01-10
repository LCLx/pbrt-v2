// Tao Du
// taodu@csail.mit.edu
// Jan 9, 2017
#include <fstream>
#include <vector>
#include "hex_mesh_element.h"

void ConvertHexToPbrt(const std::string& hex_file, const std::string& rho_file, const std::string& pbrt_file) {
  const double radius = 5e-4;
  const double pi = 3.14159265358979323846264335;
  int hex_element_num;
  std::ifstream hex_input, rho_input;
  hex_input.open(hex_file);
  rho_input.open(rho_file);
  std::ofstream pbrt_output;
  pbrt_output.open(pbrt_file);
  hex_input >> hex_element_num;
  rho_input >> hex_element_num;
  const std::vector<std::array<int, 2>> default_edges = {
    { 0, 1 },
    { 1, 3 },
    { 3, 2 },
    { 2, 0 },
    { 4, 5 },
    { 5, 7 },
    { 7, 6 },
    { 6, 4 },
    { 1, 5 },
    { 3, 7 },
    { 2, 6 },
    { 0, 4 },
  };
  for (int i = 0; i < hex_element_num; ++i) {
    std::array<double, 24> point;
    for (int j = 0; j < 24; ++j)
      hex_input >> point[j];
    const HexMeshElement element(point);
    double rho;
    rho_input >> rho;
    // Write data to pbrt.
    pbrt_output << "AttributeBegin" << std::endl;
    const double color = 1 - rho * rho;
    pbrt_output << "Material \"glass\" \"rgb Kr\" [0.0 0.0 0.0] \"rgb Kt\" [" << color << " " << color << " " << color << "]"
      << " \"float index\" [1.0]" << std::endl;
    pbrt_output << "Shape \"trianglemesh\" \"point P\" [" << std::endl;
    for (int j = 0; j < element.NumOfPoints(); ++j) {
      pbrt_output << element.PointX(j) << " " << element.PointY(j) << " " << element.PointZ(j) << std::endl;
    }
    pbrt_output << "]" << std::endl;
    pbrt_output << "\"integer indices\" [" << std::endl
      << "0 3 1" << std::endl
      << "3 0 2" << std::endl
      << "4 5 7" << std::endl
      << "4 7 6" << std::endl
      << "0 1 5" << std::endl
      << "0 5 4" << std::endl
      << "3 2 7" << std::endl
      << "7 2 6" << std::endl
      << "1 3 7" << std::endl
      << "1 7 5" << std::endl
      << "2 0 6" << std::endl
      << "6 0 4" << std::endl
      << "]" << std::endl;
    pbrt_output << "AttributeEnd" << std::endl;

    // Write edges.
    pbrt_output << "AttributeBegin" << std::endl;
    pbrt_output << "Material \"plastic\" \"rgb Kd\" [8 8 8] \"rgb Ks\" [10 10 10]" << std::endl;
    for (int j = 0; j < 12; ++j) {
      const int index0 = default_edges[j][0], index1 = default_edges[j][1];
      const double x0 = element.PointX(index0), x1 = element.PointX(index1);
      const double y0 = element.PointY(index0), y1 = element.PointY(index1);
      const double z0 = element.PointZ(index0), z1 = element.PointZ(index1);
      const double x01 = x1 - x0, y01 = y1 - y0, z01 = z1 - z0;
      const double length = std::sqrt(x01 * x01 + y01 * y01 + z01 * z01);
      pbrt_output << "TransformBegin" << std::endl;
      pbrt_output << "Translate " << x0 << " " << y0 << " " << z0 << std::endl;
      const double angle = std::acos(z01 / length);
      if (std::abs(angle - pi) < 1e-6) {
        pbrt_output << "Rotate 180 1 0 0" << std::endl;
      } else if (std::abs(angle) > 1e-6) {
        const double n_length = std::sqrt(x01 * x01 + y01 * y01);
        const double nx = -y01 / n_length, ny = x01 / n_length;
        pbrt_output << "Rotate " << angle / pi * 180.0 << " " << nx << " " << ny << " 0" << std::endl;
      }
      pbrt_output << "Scale 1 1 " << length << std::endl;
      pbrt_output << "Shape \"cylinder\" \"float radius\" [" << radius << "] \"float zmin\" [0.0] \"float zmax\" [1.0]" << std::endl;
      pbrt_output << "TransformEnd" << std::endl;
    }
    pbrt_output << "AttributeEnd" << std::endl;
  }
  hex_input.close();
  rho_input.close();
  pbrt_output.close();
}