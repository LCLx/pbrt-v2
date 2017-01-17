#!/usr/bin/python

import os
import sys

n_start = int(sys.argv[1])
n_end = int(sys.argv[2])
scale_factor = float(sys.argv[3])
n_step = 1
if len(sys.argv) >= 5:
  n_step = int(sys.argv[4])

# Now generate pbrt scripts.
root_folder = '/home/ubuntu/external/j15_bridge3_ale_64_3_defo/rendering'

if not os.path.exists(root_folder):
  os.system('mkdir -p %s' % root_folder)

for i in range(n_start, n_end + 1, n_step):
  # Rendering.
  os.system('../../../src/bin/pbrt --ncores 36 --outfile %s/j15_bridge3_ale_64_3_defo_%03d.exr %d/j15_bridge3_ale_64_3_defo.pbrt' % (root_folder, i, i))
  os.system('../../../src/bin/exrtotiff -scale %f %s/j15_bridge3_ale_64_3_defo_%03d.exr %s/j15_bridge3_ale_64_3_defo_%03d.tiff' % (scale_factor, root_folder, i, root_folder, i))
