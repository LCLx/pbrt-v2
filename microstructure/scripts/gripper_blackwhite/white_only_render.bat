..\..\..\bin\pbrt.exe --ncores 12 --outfile white_only.exr gripper_blackwhite.pbrt
..\..\..\bin\exrtotiff -scale 2 white_only.exr white_only.tiff