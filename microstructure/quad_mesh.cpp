// Tao Du
// taodu@csail.mit.edu
// Jan 14, 2017
#include <fstream>
#include "quad_mesh.h"

QuadMesh::QuadMesh(const std::string& mesh_file_name, const std::string& rho_file_name) {
  std::ifstream mesh_file, rho_file;
  mesh_file.open(mesh_file_name, std::ios::binary);
  rho_file.open(rho_file_name, std::ios::binary);
  int point_num, cell_num;
  mesh_file.read(reinterpret_cast<char*>(&point_num), sizeof(int));
  rho_file.read(reinterpret_cast<char*>(&cell_num), sizeof(int));
  assert(point_num = cell_num * 4);
  // Read data.
  Eigen::Vector2d* point_data = new Eigen::Vector2d[point_num];
  double* rho_data = new double[cell_num];
  mesh_file.read(reinterpret_cast<char*>(point_data), sizeof(Eigen::Vector2d) * point_num);
  rho_file.read(reinterpret_cast<char*>(rho_data), sizeof(double) * cell_num);
  point_ = Eigen::Matrix2Xd::Zero(2, point_num);
  for (int i = 0; i < point_num; ++i)
    point_.col(i) = point_data[i];
  rho_ = Eigen::VectorXi::Zero(cell_num);
  for (int i = 0; i < cell_num; ++i)
    rho_(i) = static_cast<int>(rho_data[i]);

  delete[] point_data;
  delete[] rho_data;
  mesh_file.close();
  rho_file.close();
}

const int QuadMesh::CellType(const int i) const {
  assert(i >= 0 && i < NumOfCell());
  return rho_(i);
}

const Eigen::Matrix<double, 3, 8> QuadMesh::ExtrudedQuadElement(const int i, const double depth) const {
  Eigen::Matrix<double, 3, 8> point_3d = Eigen::Matrix<double, 3, 8>::Zero();
  for (int j = 0; j < 4; ++j) {
    const Eigen::Vector2d point = point_.col(4 * i + j);
    point_3d.col(j) = Eigen::Vector3d(point.x(), point.y(), 0.0);
    point_3d.col(j + 4) = Eigen::Vector3d(point.x(), point.y(), depth);
  }
  return point_3d;
}

void QuadMesh::Translate(const Eigen::Vector2d& translate_vector) {
  point_.colwise() += translate_vector;
}

void QuadMesh::Scale(const double scale_factor) {
  point_ *= scale_factor;
}

const Eigen::Vector2d QuadMesh::BoundingBoxMin() const {
  return point_.rowwise().minCoeff();
}

const Eigen::Vector2d QuadMesh::BoundingBoxMax() const {
  return point_.rowwise().maxCoeff();
}

void QuadMesh::Normalize() {
  const Eigen::Vector2d diagonal = BoundingBoxMax() - BoundingBoxMin();
  Scale(1.0 / diagonal.maxCoeff());
  Translate(-BoundingBoxMin());
}

void QuadMesh::ToPBRT(const std::string& pbrt_file_name, const double extruded_depth) const {
  const int cell_num = NumOfCell();
  std::ofstream pbrt_file;
  pbrt_file.open(pbrt_file_name);
  pbrt_file << "AttributeBegin" << std::endl;
  pbrt_file << "Material \"matte\" \"rgb Kd\" [0.75 0.75 0.75]" << std::endl;
  for (int i = 0; i < cell_num; ++i) {
    if (!CellType(i)) continue;
    // Draw this cube as white.
    const Eigen::Matrix<double, 3, 8> points = ExtrudedQuadElement(i, extruded_depth);
    pbrt_file << "Shape \"trianglemesh\"" << std::endl;
    pbrt_file << "\"integer indices\" [" << std::endl
      << "4 6 7" << std::endl
      << "4 7 5" << std::endl
      << "0 3 2" << std::endl
      << "0 1 3" << std::endl
      << "1 3 7" << std::endl
      << "1 7 5" << std::endl
      << "2 0 6" << std::endl
      << "6 0 4" << std::endl
      << "2 3 7" << std::endl
      << "2 7 6" << std::endl
      << "0 5 1" << std::endl
      << "0 4 5" << std::endl
      << "]" << std::endl;
    pbrt_file << "\"point P\" [" << std::endl;
    for (int j = 0; j < 8; ++j) {
      pbrt_file << points(0, j) << " " << points(1, j) << " " << points(2, j) << std::endl;
    }
    pbrt_file << "]" << std::endl;
  }
  pbrt_file << "AttributeEnd" << std::endl;
  pbrt_file.close();
}