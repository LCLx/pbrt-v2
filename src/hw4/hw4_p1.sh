#!/bin/bash

../bin/exrtotiff lines.exr lines.tiff
../bin/exrtotiff manykilleroos_ref.exr manykilleroos_ref.tiff
../bin/exrtotiff occluder_ref.exr occluder_ref.tiff
../bin/exrtotiff occluder_shadow_ref.exr occluder_shadow_ref.tiff

../bin/pbrt occluder_ld_1.pbrt
../bin/pbrt occluder_ld_2.pbrt
../bin/pbrt occluder_ld_4.pbrt
../bin/pbrt occluder_ld_8.pbrt
../bin/pbrt occluder_ld_16.pbrt
../bin/pbrt occluder_ld_32.pbrt
../bin/pbrt occluder_ld_64.pbrt
../bin/pbrt occluder_rd_1.pbrt
../bin/pbrt occluder_rd_2.pbrt
../bin/pbrt occluder_rd_4.pbrt
../bin/pbrt occluder_rd_8.pbrt
../bin/pbrt occluder_rd_16.pbrt
../bin/pbrt occluder_rd_32.pbrt
../bin/pbrt occluder_rd_64.pbrt

../bin/exrtotiff occluder_ld_1.exr occluder_ld_1.tiff
../bin/exrtotiff occluder_ld_2.exr occluder_ld_2.tiff
../bin/exrtotiff occluder_ld_4.exr occluder_ld_4.tiff
../bin/exrtotiff occluder_ld_8.exr occluder_ld_8.tiff
../bin/exrtotiff occluder_ld_16.exr occluder_ld_16.tiff
../bin/exrtotiff occluder_ld_32.exr occluder_ld_32.tiff
../bin/exrtotiff occluder_ld_64.exr occluder_ld_64.tiff
../bin/exrtotiff occluder_rd_1.exr occluder_rd_1.tiff
../bin/exrtotiff occluder_rd_2.exr occluder_rd_2.tiff
../bin/exrtotiff occluder_rd_4.exr occluder_rd_4.tiff
../bin/exrtotiff occluder_rd_8.exr occluder_rd_8.tiff
../bin/exrtotiff occluder_rd_16.exr occluder_rd_16.tiff
../bin/exrtotiff occluder_rd_32.exr occluder_rd_32.tiff
../bin/exrtotiff occluder_rd_64.exr occluder_rd_64.tiff

#	difference
rm hw4_data.txt
../bin/exrdiff occluder_ld_1.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_ld_2.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_ld_4.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_ld_8.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_ld_16.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_ld_32.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_ld_64.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_rd_1.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_rd_2.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_rd_4.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_rd_8.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_rd_16.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_rd_32.exr occluder_ref.exr >> hw4_data.txt
../bin/exrdiff occluder_rd_64.exr occluder_ref.exr >> hw4_data.txt
