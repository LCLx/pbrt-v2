#!/usr/bin/python

import math
import os
import sys

# Arguments: Number of frames Sampling number
n_frame = int(sys.argv[1])
delta = 360.0 / n_frame
n_sample = int(sys.argv[2])

# Now generate pbrt scripts.
if not os.path.exists('rendering'):
  os.makedirs('rendering')

for i in range(n_frame):
  f = open('hex.pbrt', 'w')
  f.write('Film "image" "integer xresolution" [1600] "integer yresolution" [1200]\n')
  f.write('LookAt 3.5 4 -10 0.0 0.0 0.0 0.0 1.0 0.0\n')
  f.write('Camera "perspective" "float fov" [15]\n')
  f.write('Sampler "lowdiscrepancy" "integer pixelsamples" [%d]\n' % n_sample)

  f.write('WorldBegin\n')
  f.write('AttributeBegin\n')
  f.write('  LightSource "infinite" "rgb L" [.01 .01 .01]\n')
  f.write('AttributeEnd\n')


  f.write('AttributeBegin\n')
  f.write('  AreaLightSource "diffuse" "rgb L" [2 2 2]\n')
  f.write('  Translate 0 5 0\n')
  f.write('  Rotate 90 1 0 0\n')
  f.write('  Shape "disk" "float radius" [6.0]\n')
  f.write('AttributeEnd\n')

  f.write('AttributeBegin\n')
  f.write('Material "matte" "rgb Kd" [1.4 1.5 1.7]\n')
  f.write('TransformBegin\n')
  f.write('  Scale 100 1 100\n')
  f.write('  Shape "trianglemesh" "integer indices" [0 1 2 2 1 3]\n')
  f.write('  "point P" [-1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1]\n')
  f.write('TransformEnd\n')
  f.write('AttributeEnd\n')
  f.write('TransformBegin\n')
  f.write('  Rotate %f 0 1 0\n' % (i * delta))
  f.write('  Scale 3 3 3\n')
  f.write('  Translate -0.6 -0.25 0\n')
  f.write('  Include "hex_mesh.pbrt"\n')
  f.write('TransformEnd\n')
  f.write('WorldEnd\n')
  f.close()

  # Rendering.
  os.system('../../../src/bin/pbrt --ncores 36 --outfile rendering/hex_%03d.exr hex.pbrt' % i)
  os.system('../../../src/bin/exrtotiff rendering/hex_%03d.exr rendering/hex_%03d.tiff' % (i, i))
