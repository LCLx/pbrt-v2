#!/usr/bin/python

import os
import sys

n_start = int(sys.argv[1])
n_end = int(sys.argv[2])
scale_factor = float(sys.argv[3])
n_step = 1
if len(sys.argv) >= 5:
  n_step = int(sys.argv[4])
print n_step

# Now generate pbrt scripts.
name = 'j15_canti_ale_2d_4_defo_final'
root_folder = '/home/ubuntu/external/' + name + '/rendering'

if not os.path.exists(root_folder):
  os.system('mkdir -p %s' % root_folder)

for i in range(n_start, n_end + 1, n_step):
  # Rendering.
  os.system('../../../src/bin/pbrt --ncores 36 --outfile %s/%s_%03d.exr %d/%s.pbrt' % (root_folder, name, i, i, name))
  os.system('../../../src/bin/exrtotiff -scale %f %s/%s_%03d.exr %s/%s_%03d.tiff' % (scale_factor, root_folder, name, i, root_folder, name, i))

