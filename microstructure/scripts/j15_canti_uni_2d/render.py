#!/usr/bin/python

import os
import sys

n_start = int(sys.argv[1])
n_end = int(sys.argv[2])
scale_factor = float(sys.argv[3])

# Now generate pbrt scripts.
root_folder = '/home/ubuntu/external/j15_bridge3_uni_64/rendering'

example_name = 'j15_canti_uni_2d'

if not os.path.exists(root_folder):
  os.system('mkdir -p %s' % root_folder)

for i in range(n_start, n_end + 1):
  # Rendering.
  os.system('../../../src/bin/pbrt --ncores 36 --outfile %s/%s_%03d.exr %d/%s.pbrt' % (root_folder, example_name, i, i, example_name))
  os.system('../../../src/bin/exrtotiff -scale %f %s/%s_%03d.exr %s/%s_%03d.tiff' % (scale_factor, root_folder, example_name, i, root_folder, example_name, i))
