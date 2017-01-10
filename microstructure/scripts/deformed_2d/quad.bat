..\..\x64\Release\microstructure.exe -q hex_mesh.txt rho.txt quad_mesh.pbrt
..\..\..\bin\pbrt.exe --ncores 12 --outfile quad.exr quad.pbrt
..\..\..\bin\exrtotiff.exe -scale 1 quad.exr quad.tiff