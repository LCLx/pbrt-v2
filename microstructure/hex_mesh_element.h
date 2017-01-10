// Tao Du
// taodu@csail.mit.edu
// Jan 9, 2017
#ifndef _HEX_MESH_ELEMENT_H_
#define _HEX_MESH_ELEMENT_H_

#include <array>

class HexMeshElement {
public:
  HexMeshElement(const std::array<double, 24>& points)
    : points_(points) {}

  const int NumOfPoints() const { return 8; }
  const std::array<double, 3> Point(const int i) const {
    return std::array<double, 3>{points_[3 * i], points_[3 * i + 1], points_[3 * i + 2]};
  }
  const double PointX(const int i) const {
    return points_[3 * i];
  }
  const double PointY(const int i) const {
    return points_[3 * i + 1];
  }
  const double PointZ(const int i) const {
    return points_[3 * i + 2];
  }

private:
  std::array<double, 24> points_;
};


#endif