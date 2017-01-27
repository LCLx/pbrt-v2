// Tao Du
// taodu@csail.mit.edu
// Jan 27, 2017
#include <iostream>
#include <limits>
#include "helper.h"

// Usage: density_file_name pbrt_file_name threshold
int main(int argc, char* argv[]) {
  if (argc == 1) {
    std::cout << "Usage: main.exe -m {solid|volume|mixed} -d <density_file_name> -p <pbrt_file_name> [-t <threshold>]" << std::endl;
    return 0;
  }
  int current_argument_index = 1;
  enum Mode {kSolid = 0, kVolume, kMixed, kInvalid};
  Mode mode = kInvalid;
  std::string density_file_name(""), pbrt_file_name("");
  double threshold = std::numeric_limits<double>::infinity();
  while (current_argument_index < argc) {
    const std::string option(argv[current_argument_index++]);
    if (option == "-m") {
      if (mode != kInvalid) {
        std::cout << "Warning: duplicated -m." << std::endl;
      }
      if (current_argument_index >= argc) {
        std::cout << "Error: missing arguments after -m." << std::endl;
        return 0;
      }
      const std::string mode_string(argv[current_argument_index++]);
      if (mode_string == "solid") mode = kSolid;
      else if (mode_string == "volume") mode = kVolume;
      else if (mode_string == "mixed") mode = kMixed;
      else {
        std::cout << "Error: invalid mode arguments." << std::endl;
        return 0;
      }
    } else if (option == "-d") {
      if (!density_file_name.empty()) {
        std::cout << "Warning: duplicated -d." << std::endl;
      }
      if (current_argument_index >= argc) {
        std::cout << "Error: missing arguments after -d." << std::endl;
        return 0;
      }
      density_file_name = std::string(argv[current_argument_index++]);
    } else if (option == "-p") {
      if (!pbrt_file_name.empty()) {
        std::cout << "Warning: duplicated -p." << std::endl;
      }
      if (current_argument_index >= argc) {
        std::cout << "Error: missing arguments after -p." << std::endl;
        return 0;
      }
      pbrt_file_name = std::string(argv[current_argument_index++]);
    } else if (option == "-t") {
      if (threshold != std::numeric_limits<double>::infinity()) {
        std::cout << "Warning: duplicated -t." << std::endl;
      }
      if (current_argument_index >= argc) {
        std::cout << "Error: missing arguments after -t." << std::endl;
        return 0;
      }
      threshold = std::atof(argv[current_argument_index++]);
    } else {
      std::cout << "Error: invalid arguments." << std::endl;
      return 0;
    }
  }
  // By default we use 0.0 as the threshold.
  if (threshold == std::numeric_limits<double>::infinity()) threshold = 0.0;
  if (mode == kInvalid) {
    std::cout << "Error: missing -m." << std::endl;
    return 0;
  }
  if (density_file_name.empty()) {
    std::cout << "Error: missing -d." << std::endl;
    return 0;
  }
  if (pbrt_file_name.empty()) {
    std::cout << "Error: missing -p." << std::endl;
    return 0;
  }
  std::vector<double> density;
  std::ifstream density_file(density_file_name, std::ios::binary);
  ReadDensity(density_file, density);
  std::ofstream pbrt_file(pbrt_file_name);
  if (mode == kSolid) {
    WriteDensityToSolidPbrt(pbrt_file, density, threshold);
  } else if (mode == kVolume) {
    // Volume.
    WriteDensityToVolumePbrt(pbrt_file, density, threshold);
  } else {
    // Mixed.
    WriteDensityToSolidPbrt(pbrt_file, density, threshold);
    WriteDensityToVolumePbrt(pbrt_file, density, threshold);
  }
  return 0;
}