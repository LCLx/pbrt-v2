#!/usr/bin/python

import os
import sys

n_start = int(sys.argv[1])
n_end = int(sys.argv[2])
scale_factor = float(sys.argv[3])

# Now generate pbrt scripts.
root_folder = '/home/ubuntu/external/j14_sin_ani_2/linear_interp/rendering'

if not os.path.exists(root_folder):
  os.system('mkdir -p %s' % root_folder)

for i in range(n_start, n_end + 1):
  # Rendering.
  os.system('../../../../src/bin/pbrt --ncores 36 --outfile %s/j14_sin_ani_2_%03d.exr %d/jj14_sin_ani_2.pbrt' % (root_folder, i, i))
  os.system('../../../../src/bin/exrtotiff -scale %f %s/j14_sin_ani_2_%03d.exr %s/j14_sin_ani_2_%03d.tiff' % (scale_factor, root_folder, i, root_folder, i))
