// Tao Du
// taodu@csail.mit.edu
// Jan 12, 2017
#include <fstream>
#include "hex_mesh.h"

HexMesh::HexMesh(const std::string& lattice_file, const std::string& displacement_file,
  const std::string& material_file, const std::string& sing_point_file) {
  // Read cell_counts, dx and domain_min from lattice file.
  std::ifstream lattice;
  lattice.open(lattice_file, std::ios::binary);
  Eigen::Vector3i cell_count;
  lattice.read(reinterpret_cast<char*>(&cell_count), sizeof(Eigen::Vector3i));
  node_count_ = cell_count.array() + 1;
  lattice.read(reinterpret_cast<char*>(&dx_), sizeof(double));
  lattice.read(reinterpret_cast<char*>(&domain_min_), sizeof(Eigen::Vector3d));
  lattice.close();

  // Read displacement file.
  std::ifstream displacement;
  displacement.open(displacement_file, std::ios::binary);
  Eigen::Vector3i displacement_count;
  displacement.read(reinterpret_cast<char*>(&displacement_count), sizeof(Eigen::Vector3i));
  for (int i = 0; i < 3; ++i)
    assert(displacement_count(i) == node_count_(i));
  const int total_node_num = node_count_.prod();
  Eigen::Vector3d* displacement_data = new Eigen::Vector3d[total_node_num];
  displacement.read(reinterpret_cast<char*>(displacement_data), sizeof(Eigen::Vector3d) * total_node_num);
  displacement_ = Eigen::Matrix3Xd::Zero(3, total_node_num);
  for (int i = 0; i < total_node_num; ++i) {
    displacement_.col(i) = displacement_data[i];
  }
  delete[] displacement_data;
  displacement.close();

  // Read material.
  std::ifstream material;
  material.open(material_file, std::ios::binary);
  Eigen::Vector3i material_count;
  material.read(reinterpret_cast<char*>(&material_count), sizeof(Eigen::Vector3i));
  for (int i = 0; i < 3; ++i)
    assert(material_count(i) + 1 == node_count_(i));
  const int total_cell_num = material_count.prod();
  int* material_data = new int[total_cell_num];
  material.read(reinterpret_cast<char*>(material_data), sizeof(int) * total_cell_num);
  material_ = Eigen::VectorXi::Zero(total_cell_num);
  for (int i = 0; i < total_cell_num; ++i) {
    material_(i) = material_data[i];
  }
  delete[] material_data;
  material.close();

  // Read sing point.
  std::ifstream sing_point;
  sing_point.open(sing_point_file, std::ios::binary);
  int sing_point_count;
  sing_point.read(reinterpret_cast<char*>(&sing_point_count), sizeof(int));
  Eigen::Vector3d* sing_point_data = new Eigen::Vector3d[sing_point_count];
  sing_point.read(reinterpret_cast<char*>(sing_point_data), sizeof(Eigen::Vector3d) * sing_point_count);
  sing_point_ = Eigen::Matrix3Xd::Zero(3, sing_point_count);
  for (int i = 0; i < sing_point_count; ++i) {
    sing_point_.col(i) = sing_point_data[i];
  }
  delete[] sing_point_data;
  sing_point.close();
}

void HexMesh::Translate(const Eigen::Vector3d& translate_vector) {
  domain_min_ += translate_vector;
  sing_point_.colwise() += translate_vector;
}

void HexMesh::Scale(const double scale_factor) {
  domain_min_ *= scale_factor;
  dx_ *= scale_factor;
  displacement_ *= scale_factor;
  sing_point_ *= scale_factor;
}

const Eigen::Vector3d HexMesh::BoundingBoxMin() const {
  return domain_min_;
}

const Eigen::Vector3d HexMesh::BoundingBoxMax() const {
  return domain_min_.array() + (node_count_.cast<double>().array() - 1) * dx_;
}

void HexMesh::Normalize() {
  const Eigen::Vector3d diagonal = BoundingBoxMax() - BoundingBoxMin();
  Scale(1.0 / diagonal.maxCoeff());
  Translate(-BoundingBoxMin());
}

void HexMesh::ToPBRT(const std::string& pbrt_file) const {
  std::ofstream pbrt_output;
  pbrt_output.open(pbrt_file);

  const double radius = 1e-3;
  const double point_radius = 5e-3;
  const Eigen::Matrix<int, 2, 12> default_edges = (Eigen::Matrix<int, 2, 12>()
    << 0, 1, 3, 2, 4, 5, 7, 6, 1, 3, 2, 0,
       1, 3, 2, 0, 5, 7, 6, 4, 5, 7, 6, 4).finished();
  const int cell_x_num = NumOfCellX();
  const int cell_y_num = NumOfCellY();
  const int cell_z_num = NumOfCellZ();
  for (int i = 0; i < cell_x_num; ++i) {
    for (int j = 0; j < cell_y_num; ++j) {
      for (int k = 0; k < cell_z_num; ++k) {
        if (MaterialType(i, j, k) < 0) continue;
        const Eigen::Matrix<double, 3, 8> element = HexElement(i, j, k);
        // Write data to pbrt.
        pbrt_output << "AttributeBegin" << std::endl;
        const double color = 0.75;  // TODO: modify this when rho is available.
        pbrt_output << "Material \"glass\" \"rgb Kr\" [0.5 0.5 0.5] \"rgb Kt\" [" << color << " " << color << " " << color << "]"
          << " \"float index\" [1.0]" << std::endl;
        pbrt_output << "Shape \"trianglemesh\" \"point P\" [" << std::endl;
        for (int j = 0; j < 8; ++j) {
          pbrt_output << element(0, j) << " " << element(1, j) << " " << element(2, j) << std::endl;
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
          const int index0 = default_edges(0, j), index1 = default_edges(1, j);
          const Eigen::Vector3d p0 = element.col(index0), p1 = element.col(index1);
          const Eigen::Vector3d p01 = p1 - p0;
          pbrt_output << "TransformBegin" << std::endl;
          pbrt_output << "Translate " << p0.x() << " " << p0.y() << " " << p0.z() << std::endl;
          const Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), p01);
          const Eigen::AngleAxisd angle_axis(q);
          const Eigen::Vector3d axis = angle_axis.axis();
          pbrt_output << "Rotate " << angle_axis.angle() / M_PI * 180.0 << " " << axis.x() << " " << axis.y() << " " << axis.z() << std::endl;
          pbrt_output << "Scale 1 1 " << p01.norm() << std::endl;
          pbrt_output << "Shape \"cylinder\" \"float radius\" [" << radius << "] \"float zmin\" [0.0] \"float zmax\" [1.0]" << std::endl;
          pbrt_output << "TransformEnd" << std::endl;
        }
        pbrt_output << "AttributeEnd" << std::endl;
      }
    }
  }

  // Draw sing point.
  const int sing_point_num = static_cast<int>(sing_point_.cols());
  pbrt_output << "AttributeBegin" << std::endl;
  pbrt_output << "Material \"matte\" \"rgb Kd\" [0.75 0.05 0.05]" << std::endl;
  for (int i = 0; i < sing_point_num; ++i) {
    const Eigen::Vector3d center = sing_point_.col(i);
    pbrt_output << "TransformBegin" << std::endl;
    pbrt_output << "Translate " << center.x() << " " << center.y() << " " << center.z() << std::endl;
    pbrt_output << "Shape \"sphere\" \"float radius\" [" << point_radius << "]" << std::endl;
    pbrt_output << "TransformEnd" << std::endl;
  }
  pbrt_output << "AttributeEnd" << std::endl;
  pbrt_output.close();
}

const int HexMesh::CellSubToIdx(const int i, const int j, const int k) const {
  assert(i >= 0 && i < NumOfCellX());
  assert(j >= 0 && j < NumOfCellY());
  assert(k >= 0 && k < NumOfCellZ());
  return i * NumOfCellY() * NumOfCellZ() + j * NumOfCellZ() + k;
}

const int HexMesh::NodeSubToIdx(const int i, const int j, const int k) const {
  assert(i >= 0 && i < NumOfNodeX());
  assert(j >= 0 && j < NumOfNodeY());
  assert(k >= 0 && k < NumOfNodeZ());
  return i * NumOfNodeY() * NumOfNodeZ() + j * NumOfNodeZ() + k;
}

const Eigen::Matrix<double, 3, 8> HexMesh::HexElement(const int i, const int j, const int k) const {
  const Eigen::Vector3d corner_min = domain_min_ + Eigen::Vector3d(i, j, k) * dx_;
  Eigen::Matrix<double, 3, 8> hex_element;
  hex_element.col(0) = displacement_.col(NodeSubToIdx(i, j, k));
  hex_element.col(1) = Eigen::Vector3d(0, 0, 1) * dx_ + displacement_.col(NodeSubToIdx(i, j, k + 1));
  hex_element.col(2) = Eigen::Vector3d(0, 1, 0) * dx_ + displacement_.col(NodeSubToIdx(i, j + 1, k));
  hex_element.col(3) = Eigen::Vector3d(0, 1, 1) * dx_ + displacement_.col(NodeSubToIdx(i, j + 1, k + 1));
  hex_element.col(4) = Eigen::Vector3d(1, 0, 0) * dx_ + displacement_.col(NodeSubToIdx(i + 1, j, k));
  hex_element.col(5) = Eigen::Vector3d(1, 0, 1) * dx_ + displacement_.col(NodeSubToIdx(i + 1, j, k + 1));
  hex_element.col(6) = Eigen::Vector3d(1, 1, 0) * dx_ + displacement_.col(NodeSubToIdx(i + 1, j + 1, k));
  hex_element.col(7) = Eigen::Vector3d(1, 1, 1) * dx_ + displacement_.col(NodeSubToIdx(i + 1, j + 1, k + 1));
  hex_element.colwise() += corner_min;
  return hex_element;
}

const int HexMesh::MaterialType(const int i, const int j, const int k) const {
  return material_(CellSubToIdx(i, j, k));
}