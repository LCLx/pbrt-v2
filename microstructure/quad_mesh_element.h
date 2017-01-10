// Tao Du
// taodu@csail.mit.edu
// Jan 9, 2017
#ifndef _QUAD_MESH_ELEMENT_H_
#define _QUAD_MESH_ELEMENT_H_

#include <array>

class QuadMeshElement {
public:
  QuadMeshElement(const std::array<double, 8>& points)
    : points_(points) {}

  const int NumOfPoints() const { return 4; }
  const std::array<double, 2> Point(const int i) const {
    return std::array<double, 2>{points_[2 * i], points_[2 * i + 1]};
  }
  const double PointX(const int i) const {
    return points_[2 * i];
  }
  const double PointY(const int i) const {
    return points_[2 * i + 1];
  }

private:
  std::array<double, 8> points_;
};

#endif