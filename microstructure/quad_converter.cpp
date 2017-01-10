// Tao Du
// taodu@csail.mit.edu
// Jan 9, 2017
#include <fstream>
#include <vector>
#include "quad_mesh_element.h"

void ConvertQuadToPbrt(const std::string& quad_file, const std::string& rho_file, const std::string& pbrt_file) {
  const double radius = 1e-3;
  const double pi = 3.14159265358979323846264335;
  int quad_element_num;
  std::ifstream quad_input, rho_input;
  quad_input.open(quad_file);
  rho_input.open(rho_file);
  std::ofstream pbrt_output;
  pbrt_output.open(pbrt_file);
  quad_input >> quad_element_num;
  rho_input >> quad_element_num;
  const std::vector<std::array<int, 2>> default_edges = {
    { 0, 1 },
    { 0, 2 },
    { 1, 3 },
    { 2, 3 }
  };
  for (int i = 0; i < quad_element_num; ++i) {
    std::array<double, 8> point;
    for (int j = 0; j < 8; ++j)
      quad_input >> point[j];
    const QuadMeshElement element(point);
    double rho;
    rho_input >> rho;
    // Write data to pbrt.
    pbrt_output << "AttributeBegin" << std::endl;
    pbrt_output << "Material \"glass\" \"rgb Kr\" [0.0 0.0 0.0] \"rgb Kt\" [" << 1 - rho * rho * rho << " 0.0 0.0]"
      << " \"float index\" [1.0]" << std::endl;
    pbrt_output << "Shape \"trianglemesh\" \"point P\" [" << std::endl;
    for (int j = 0; j < element.NumOfPoints(); ++j) {
      pbrt_output << element.PointX(j) << " " << element.PointY(j) << " 0.0" << std::endl;
    }
    pbrt_output << "]" << std::endl;
    pbrt_output << "\"integer indices\" [0 1 3 0 3 2]" << std::endl;
    pbrt_output << "AttributeEnd" << std::endl;

    // Write edges.
    pbrt_output << "AttributeBegin" << std::endl;
    pbrt_output << "Material \"plastic\" \"rgb Kd\" [0.1 0.8 0.1] \"rgb Ks\" [0.3 1.0 0.2]" << std::endl;
    for (int j = 0; j < 4; ++j) {
      const int index0 = default_edges[j][0], index1 = default_edges[j][1];
      const double x0 = element.PointX(index0), x1 = element.PointX(index1);
      const double y0 = element.PointY(index0), y1 = element.PointY(index1);
      const double x01 = x1 - x0, y01 = y1 - y0;
      const double length = std::sqrt(x01 * x01 + y01 * y01);
      const double angle = std::acos(x01 / length) * (y01 > 0.0 ? 1.0 : -1.0);
      pbrt_output << "TransformBegin" << std::endl;
      pbrt_output << "Translate " << x0 << " " << y0 << " 0.0" << std::endl;
      pbrt_output << "Rotate " << angle / pi * 180.0 << " 0 0 1" << std::endl;
      pbrt_output << "Rotate 90 0 1 0" << std::endl;
      pbrt_output << "Scale 1 1 " << length << std::endl;
      pbrt_output << "Shape \"cylinder\" \"float radius\" [" << radius << "] \"float zmin\" [0.0] \"float zmax\" [1.0]" << std::endl;
      pbrt_output << "TransformEnd" << std::endl;
    }
    pbrt_output << "AttributeEnd" << std::endl;
  }
  quad_input.close();
  rho_input.close();
  pbrt_output.close();
}