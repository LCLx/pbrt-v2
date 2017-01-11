..\..\x64\Release\microstructure.exe -q hex_mesh rho quad_mesh.pbrt
..\..\..\bin\pbrt.exe --ncores 12 --outfile quad.exr quad.pbrt
..\..\..\bin\exrtotiff.exe -scale 1 quad.exr quad.tiff