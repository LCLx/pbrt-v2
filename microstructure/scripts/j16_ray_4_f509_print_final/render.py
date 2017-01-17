#!/usr/bin/python

import math
import os
import sys

# Arguments: Number of frames Sampling number
n_frame = int(sys.argv[1])
delta = 360.0 / n_frame
n_sample = int(sys.argv[2])
scale_factor = float(sys.argv[3])

name = 'j16_ray_4_f509_print_final'

# Now generate pbrt scripts.
root_folder = '/home/ubuntu/external/' + name + '/rendering'

if not os.path.exists(root_folder):
  os.system('mkdir -p %s' % root_folder)

for i in range(0, n_frame):
  f = open(name + '_camera.pbrt', 'w')
  f.write('Film "image" "integer xresolution" [1600] "integer yresolution" [1200]\n')
  f.write('LookAt 0 4 -10 0.0 0.0 0.0 0.0 1.0 0.0\n')
  f.write('Camera "perspective" "float fov" [15]\n')
  f.write('Sampler "lowdiscrepancy" "integer pixelsamples" [%d]\n' % n_sample)
  f.write('SurfaceIntegrator "directlighting" "integer maxdepth" [3]\n')

  f.write('WorldBegin\n')
  f.write('AttributeBegin\n')
  f.write('  LightSource "infinite" "rgb L" [.01 .01 .01]\n')
  f.write('AttributeEnd\n')


  f.write('AttributeBegin\n')
  f.write('TransformBegin\n')
  f.write('  AreaLightSource "diffuse" "rgb L" [5 5 5]\n')
  f.write('  Translate -3 1 -3\n')
  f.write('  Shape "disk" "float radius" [1.0]\n')
  f.write('TransformEnd\n')
  f.write('TransformBegin\n')
  f.write('  AreaLightSource "diffuse" "rgb L" [2.5 2.5 2.5]\n')
  f.write('  Translate -3 4.5 -3\n')
  f.write('  Rotate 90 1 0 0\n')
  f.write('  Shape "disk" "float radius" [4.0]\n')
  f.write('TransformEnd\n')
  f.write('AttributeEnd\n')

  f.write('AttributeBegin\n')
  f.write('Material "matte" "rgb Kd" [0.28 0.3 0.34]\n')
  f.write('TransformBegin\n')
  f.write('  Translate 0 -0.2 0\n')
  f.write('  Scale 100 1 100\n')
  f.write('  Shape "trianglemesh" "integer indices" [0 1 2 2 1 3]\n')
  f.write('  "point P" [-1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1]\n')
  f.write('TransformEnd\n')
  f.write('AttributeEnd\n')
  f.write('TransformBegin\n')
  f.write('  Rotate %f 0 1 0\n' % (i * delta))
  f.write('  Translate 0 -1.1 0\n')
  f.write('  Scale 1.2 1.2 1.2\n')
  f.write('  Rotate 45 0 0 1\n')
  f.write('  Include "mesh.pbrt"\n')
  f.write('TransformEnd\n')
  f.write('WorldEnd\n')
  f.close()

  # Rendering.
  os.system('../../../src/bin/pbrt --ncores 36 --outfile %s/%s_%03d.exr %s_camera.pbrt' % (root_folder, name, i, name))
  os.system('../../../src/bin/exrtotiff -scale %f %s/%s_%03d.exr %s/%s_%03d.tiff' % (scale_factor, root_folder, name, i, root_folder, name, i))
