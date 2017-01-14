// Tao Du
// taodu@csail.mit.edu
// Jan 14, 2017
#ifndef _QUAD_MESH_H_
#define _QUAD_MESH_H_

#include "Eigen/Dense"

class QuadMesh {
public:
  QuadMesh(const std::string& mesh_file_name, const std::string& rho_file_name);
  ~QuadMesh() {}

  const int NumOfCell() const { return static_cast<int>(rho_.size()); }
  const int CellType(const int i) const;
  // Assume point_ is defined in XOY plane and depth is the extruded length along Z.
  const Eigen::Matrix<double, 3, 8> ExtrudedQuadElement(const int i, const double depth) const;

  void Translate(const Eigen::Vector2d& translate_vector);
  void Scale(const double scale_factor);
  const Eigen::Vector2d BoundingBoxMin() const;
  const Eigen::Vector2d BoundingBoxMax() const;
  void Normalize();

  void ToPBRT(const std::string& pbrt_file_name, const double extruded_depth) const;

private:
  Eigen::Matrix2Xd point_;
  Eigen::VectorXi rho_;
};


#endif