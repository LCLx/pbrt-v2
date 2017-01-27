// Tao Du
// taodu@csail.mit.edu
// Jan 27, 2017
#ifndef _HELPER_H_
#define _HELPER_H_

#include <string>
#include <vector>

void ReadDensity(const std::string& file_name, std::vector<double>& density);
void WriteDensityToPbrt(const std::string& file_name, const std::vector<double>& density);

#endif