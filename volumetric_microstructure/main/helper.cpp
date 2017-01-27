// Tao Du
// taodu@csail.mit.edu
// Jan 27, 2017
#include <assert.h>
#include "helper.h"

void ReadDensity(std::ifstream& fin, std::vector<double>& density) {
  int nx, ny, nz;
  fin.read(reinterpret_cast<char*>(&nx), sizeof(int));
  fin.read(reinterpret_cast<char*>(&ny), sizeof(int));
  fin.read(reinterpret_cast<char*>(&nz), sizeof(int));
  assert(nx == 64 && ny == 64 && nz == 64);
  density = std::vector<double>(nx * ny * nz, 0.0);
  fin.read(reinterpret_cast<char*>(density.data()), sizeof(double) * nx * ny * nz);
}

void WriteDensityToSolidPbrt(std::ofstream& fout, const std::vector<double>& density, const double threshold) {
  const int n = 64;
  assert(static_cast<int>(density.size()) == n * n * n);
  const double dx = 1.0 / n;
  const std::vector<std::vector<double>> points{
    {0, 0, 0},
    {0, 0, dx},
    {0, dx, 0},
    {0, dx, dx},
    {dx, 0, 0},
    {dx, 0, dx},
    {dx, dx, 0},
    {dx, dx, dx}
  };
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k) {
        const double corner_x = i * dx, corner_y = j * dx, corner_z = k * dx;
        const double d = density[i * n * n + j * n + k] * 0.9;
        if (d <= threshold) continue;
        fout << "AttributeBegin" << std::endl;
        fout << "  Material \"matte\" \"rgb Kd\" [" << d << " " << d << " " << d << "]" << std::endl
          << "Shape \"trianglemesh\"" << std::endl
          << "\"integer indices\" [" << std::endl
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
        fout << "\"point P\" [" << std::endl;
        for (int s = 0; s < 8; ++s) {
          fout << corner_x + points[s][0] << " " << corner_y + points[s][1] << " " << corner_z + points[s][2] << std::endl;
        }
        fout << "]" << std::endl;
        fout << "AttributeEnd" << std::endl;
      }
}

void WriteDensityToVolumePbrt(std::ofstream& fout, const std::vector<double>& density, const double threshold) {
  const int n = 64;
  assert(static_cast<int>(density.size()) == n * n * n);
  fout << "AttributeBegin" << std::endl;
  fout << "  Volume \"volumegrid\" \"integer nx\" 64 \"integer ny\" 64 \"integer nz\" 64" << std::endl;
  fout << "  \"float density\" [\n";
  for (const double value : density) {
    fout << (value < threshold ? value : 0.0) << std::endl;
  }
  fout << "]\n";
  fout << "\"color sigma_a\" [80 80 80]" << std::endl;
  fout << "\"color Le\" [5 5 5]" << std::endl;
  fout << "AttributeEnd" << std::endl;
}