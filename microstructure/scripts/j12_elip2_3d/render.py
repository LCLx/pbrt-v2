#!/usr/bin/python

import os
import sys

# Arguments: Number of frames Sampling number
n_frame = int(sys.argv[1])
scale_factor = float(sys.argv[2])

# Now generate pbrt scripts.
if not os.path.exists('rendering'):
  os.makedirs('rendering')

for i in range(n_frame):
  # Rendering.
  os.system('../../../src/bin/pbrt --ncores 36 --outfile %d/j12_elip2_3d_%03d.exr %d/j12_elip2_3d.pbrt' % (i, i, i))
  os.system('../../../src/bin/exrtotiff -scale %f %d/j12_elip2_3d_%03d.exr rendering/j12_elip2_3d_%03d.tiff' % (scale_factor, i, i, i))
