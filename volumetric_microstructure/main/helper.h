// Tao Du
// taodu@csail.mit.edu
// Jan 27, 2017
#ifndef _HELPER_H_
#define _HELPER_H_

#include <fstream>
#include <string>
#include <vector>

void ReadDensity(std::ifstream& fin, std::vector<double>& density);
void WriteDensityToSolidPbrt(std::ofstream& fout, const std::vector<double>& density, const double threshold);
void WriteDensityToVolumePbrt(std::ofstream& fout, const std::vector<double>& density, const double threshold);

#endif