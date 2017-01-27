// Tao Du
// taodu@csail.mit.edu
// Jan 27, 2017
#include "helper.h"

// Usage: density_file_name pbrt_file_name threshold
int main(int argc, char* argv[]) {
  const std::string density_file_name(argv[1]);
  const std::string pbrt_file_name(argv[2]);
  const double threshold = std::atof(argv[3]);
  std::vector<double> density;
  ReadDensity(density_file_name, density);
  for (double& value : density) {
    if (value < threshold) value = 0.0;
  }
  WriteDensityToPbrt(pbrt_file_name, density);
  return 0;
}