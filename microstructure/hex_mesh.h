// Tao Du
// taodu@csail.mit.edu
// Jan 12, 2017
#ifndef _HEX_MESH_H_
#define _HEX_MESH_H_

#include "Eigen/Dense"

class HexMesh {
public:
  explicit HexMesh(const std::string& lattice_file = "NULL",
    const std::string& displacement_file = "NULL",
    const std::string& material_file = "NULL",
    const std::string& lag_inf_point_file = "NULL",
    const std::string& sing_point_file = "NULL",
    const std::string& fine_intf_flag_file = "NULL",
    const std::string& f_point_file = "NULL",
    const std::string& psi_D_file = "NULL",
    const std::string& density_file = "NULL",
    const std::string& v0_file = "NULL",
    const std::string& v1_file = "NULL",
    const std::string& v2_file = "NULL",
    const std::string& v3_file = "NULL"
  );
  HexMesh(const HexMesh& other);
  ~HexMesh() {}

  HexMesh& operator*(const double t);
  HexMesh& operator+(const HexMesh& other);
  HexMesh& operator=(const HexMesh& other);

  void Translate(const Eigen::Vector3d& translate_vector);
  void Scale(const double scale_factor);
  const Eigen::Vector3d BoundingBoxMin() const;
  const Eigen::Vector3d BoundingBoxMax() const;
  void Normalize();
  void ToPBRT(const std::string& pbrt_file, const bool plot_surrounding_cells = true, const double threshold = 0.5) const;

private:
  const int NumOfNodeX() const { return node_count_.x(); }
  const int NumOfNodeY() const { return node_count_.y(); }
  const int NumOfNodeZ() const { return node_count_.z(); }

  const int NumOfCellX() const { return node_count_.x() - 1; }
  const int NumOfCellY() const { return node_count_.y() - 1; }
  const int NumOfCellZ() const { return (node_count_.z() == 1 ? 1 : node_count_.z() - 1); }

  const int CellSubToIdx(const int i, const int j, const int k) const;
  const int NodeSubToIdx(const int i, const int j, const int k) const;

  const Eigen::Matrix3Xd HexElement(const int i, const int j, const int k) const;
  const int MaterialType(const int i, const int j, const int k) const;

  static Eigen::Matrix<double, 3, 10> fine_intf_flag_colors_;

  Eigen::Vector3d domain_min_;
  Eigen::Vector3i node_count_;
  double dx_;
  Eigen::Matrix3Xd displacement_;
  Eigen::VectorXi material_;
  Eigen::Matrix3Xd lag_inf_point_;
  Eigen::Matrix3Xd sing_point_;
  Eigen::VectorXi fine_intf_flags_;
  Eigen::Matrix3Xd f_point_;
  Eigen::Matrix3Xd psi_D_;
  Eigen::Matrix3Xd density_;
};

#endif