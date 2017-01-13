// Tao Du
// taodu@csail.mit.edu
// Jan 12, 2017
#include <fstream>
#include "hex_mesh.h"

Eigen::Matrix<double, 3, 10> HexMesh::fine_intf_flag_colors_ = (Eigen::Matrix<double, 3, 10>() <<
  0.8147, 0.2, 0.2785, 0.3, 0.9572, 0.1419, 0.85, 0.0357, 0.3922, 0.6787,
  0.9058, 0.5, 0.5469, 0.9, 0.4854, 0.4218, 0.20, 0.8491, 0.6555, 0.7577,
  0.1270, 0.9, 0.9575, 0.2, 0.8003, 0.9157, 0.30, 0.9340, 0.1712, 0.7431).finished();

HexMesh::HexMesh(const std::string& lattice_file, const std::string& displacement_file,
  const std::string& material_file, const std::string& lag_inf_point_file,
  const std::string& sing_point_file, const std::string& fine_intf_flag_file) {
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

  // Read lag inf point.
  std::ifstream lag_inf_point;
  lag_inf_point_ = Eigen::Matrix3Xd::Zero(3, 0);
  if (lag_inf_point_file != "NULL") {
    lag_inf_point.open(lag_inf_point_file, std::ios::binary);
    int lag_inf_point_count;
    lag_inf_point.read(reinterpret_cast<char*>(&lag_inf_point_count), sizeof(int));
    Eigen::Vector3d* lag_inf_point_data = new Eigen::Vector3d[lag_inf_point_count];
    lag_inf_point.read(reinterpret_cast<char*>(lag_inf_point_data), sizeof(Eigen::Vector3d) * lag_inf_point_count);
    lag_inf_point_ = Eigen::Matrix3Xd::Zero(3, lag_inf_point_count);
    for (int i = 0; i < lag_inf_point_count; ++i) {
      lag_inf_point_.col(i) = lag_inf_point_data[i];
    }
    delete[] lag_inf_point_data;
    lag_inf_point.close();
  }

  // Read sing point.
  std::ifstream sing_point;
  sing_point_ = Eigen::Matrix3Xd::Zero(3, 0);
  if (sing_point_file != "NULL") {
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

  // Read fine_intf_flags.
  fine_intf_flags_ = Eigen::VectorXi::Ones(total_cell_num) * (-1);
  if (fine_intf_flag_file == "NULL") return;
  std::ifstream fine_intf_flag;
  fine_intf_flag.open(fine_intf_flag_file, std::ios::binary);
  Eigen::Vector3i fine_intf_flag_count;
  fine_intf_flag.read(reinterpret_cast<char*>(&fine_intf_flag_count), sizeof(Eigen::Vector3i));
  for (int i = 0; i < 3; ++i)
    assert(fine_intf_flag_count(i) + 1 == node_count_(i));
  double* fine_intf_flag_data = new double[total_cell_num];
  fine_intf_flag.read(reinterpret_cast<char*>(fine_intf_flag_data), sizeof(double) * total_cell_num);
  for (int i = 0; i < total_cell_num; ++i) {
    fine_intf_flags_(i) = static_cast<int>(fine_intf_flag_data[i]);
  }
  delete[] fine_intf_flag_data;
  fine_intf_flag.close();
}

void HexMesh::Translate(const Eigen::Vector3d& translate_vector) {
  domain_min_ += translate_vector;
  lag_inf_point_.colwise() += translate_vector;
  sing_point_.colwise() += translate_vector;
}

void HexMesh::Scale(const double scale_factor) {
  domain_min_ *= scale_factor;
  dx_ *= scale_factor;
  displacement_ *= scale_factor;
  lag_inf_point_ *= scale_factor;
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

  const double radius = 9e-4;
  const double point_radius = 5e-3;
  const Eigen::Matrix<int, 2, 12> default_edges = (Eigen::Matrix<int, 2, 12>()
    << 0, 1, 3, 2, 4, 5, 7, 6, 1, 3, 2, 0,
       1, 3, 2, 0, 5, 7, 6, 4, 5, 7, 6, 4).finished();
  const int cell_x_num = NumOfCellX();
  const int cell_y_num = NumOfCellY();
  const int cell_z_num = NumOfCellZ();
  const bool has_fine_intf_flags = fine_intf_flags_(0) != -1;
  for (int i = 0; i < cell_x_num; ++i) {
    for (int j = 0; j < cell_y_num; ++j) {
      for (int k = 0; k < cell_z_num; ++k) {
        if (MaterialType(i, j, k) < 0) continue;
        const Eigen::Matrix<double, 3, 8> element = HexElement(i, j, k);
        // Write data to pbrt.
        pbrt_output << "AttributeBegin" << std::endl;
        if (!has_fine_intf_flags) {
          const double color = 0.75;  // TODO: modify this when rho is available.
          pbrt_output << "Material \"glass\" \"rgb Kr\" [0.5 0.5 0.5] \"rgb Kt\" [" << color << " " << color << " " << color << "]"
            << " \"float index\" [1.0]" << std::endl;
          pbrt_output << "Shape \"trianglemesh\"" << std::endl;
          pbrt_output << "\"integer indices\" [" << std::endl
            << "1 3 7" << std::endl
            << "1 7 5" << std::endl
            << "2 0 6" << std::endl
            << "6 0 4" << std::endl
            << "]" << std::endl;
        } else {
          const Eigen::Vector3d color = fine_intf_flag_colors_.col(fine_intf_flags_(CellSubToIdx(i, j, k)));
          pbrt_output << "Material \"glass\" \"rgb Kr\" [0.5 0.5 0.5] \"rgb Kt\" [" << color.x() << " " << color.y() << " " << color.z() << "]"
            << " \"float index\" [1.0]" << std::endl;
          pbrt_output << "Shape \"trianglemesh\"" << std::endl;
          pbrt_output << "\"integer indices\" [" << std::endl
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
        }
        pbrt_output << "\"point P\" [" << std::endl;
        for (int j = 0; j < 8; ++j) {
          pbrt_output << element(0, j) << " " << element(1, j) << " " << element(2, j) << std::endl;
        }
        pbrt_output << "]" << std::endl;

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

  // Draw lag inf point.
  const int lag_inf_point_num = static_cast<int>(lag_inf_point_.cols());
  pbrt_output << "AttributeBegin" << std::endl;
  pbrt_output << "Material \"shinymetal\" \"rgb Ks\" [.1 .1 .1] \"rgb Kr\" [1 0.2 0.2]" << std::endl;
  for (int i = 0; i < lag_inf_point_num; ++i) {
    const Eigen::Vector3d center = lag_inf_point_.col(i);
    pbrt_output << "TransformBegin" << std::endl;
    pbrt_output << "Translate " << center.x() << " " << center.y() << " " << center.z() << std::endl;
    pbrt_output << "Shape \"sphere\" \"float radius\" [" << point_radius << "]" << std::endl;
    pbrt_output << "TransformEnd" << std::endl;
  }
  pbrt_output << "AttributeEnd" << std::endl;

  // Draw sing point.
  const int sing_point_num = static_cast<int>(sing_point_.cols());
  pbrt_output << "AttributeBegin" << std::endl;
  pbrt_output << "Material \"shinymetal\" \"rgb Ks\" [.1 .1 .1] \"rgb Kr\" [0.2 1 0.2]" << std::endl;
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