// Tao Du
// taodu@csail.mit.edu
// Jan 9, 2017
#include <assert.h>
#include <fstream>
#include <vector>
#include "hex_mesh_element.h"

void ConvertHexToPbrt(const std::string& hex_file, const std::string& rho_file, const std::string& pbrt_file) {
  const double radius = 1e-3;
  const double pi = 3.14159265358979323846264335;
  int hex_element_num, hex_double_num;
  std::ifstream hex_input, rho_input;
  hex_input.open(hex_file, std::ios::binary);
  rho_input.open(rho_file, std::ios::binary);
  std::ofstream pbrt_output;
  pbrt_output.open(pbrt_file);
  hex_input.read(reinterpret_cast<char*>(&hex_double_num), sizeof(int));
  rho_input.read(reinterpret_cast<char*>(&hex_element_num), sizeof(int));
  assert(hex_element_num * 24 == hex_double_num);
  double* point_data = new double[hex_double_num];
  double* rho_data = new double[hex_element_num];
  hex_input.read(reinterpret_cast<char*>(point_data), sizeof(double) * hex_double_num);
  rho_input.read(reinterpret_cast<char*>(rho_data), sizeof(double) * hex_element_num);
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
      point[j] = point_data[24 * i + j];
    const HexMeshElement element(point);
    const double rho = rho_data[i];
    // Write data to pbrt.
    pbrt_output << "AttributeBegin" << std::endl;
    const double color = 1 - rho * rho;
    pbrt_output << "Material \"glass\" \"rgb Kr\" [0.5 0.5 0.5] \"rgb Kt\" [" << color << " " << color << " " << color << "]"
      << " \"float index\" [1.0]" << std::endl;
    pbrt_output << "Shape \"trianglemesh\" \"point P\" [" << std::endl;
    for (int j = 0; j < element.NumOfPoints(); ++j) {
      pbrt_output << element.PointX(j) << " " << element.PointY(j) << " " << element.PointZ(j) << std::endl;
    }
    pbrt_output << "]" << std::endl;
    pbrt_output << "\"integer indices\" [" << std::endl
      << "1 3 7" << std::endl
      << "1 7 5" << std::endl
      << "2 0 6" << std::endl
      << "6 0 4" << std::endl
      << "]" << std::endl;
    pbrt_output << "AttributeEnd" << std::endl;

    // Write edges.
    pbrt_output << "AttributeBegin" << std::endl;
    pbrt_output << "Material \"metal\"" << std::endl;
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

  delete[] point_data;
  delete[] rho_data;

  hex_input.close();
  rho_input.close();
  pbrt_output.close();
}