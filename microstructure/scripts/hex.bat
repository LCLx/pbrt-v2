..\x64\Release\microstructure.exe -h hex_mesh.txt rho.txt hex_mesh.pbrt
..\..\bin\pbrt.exe --ncores 12 --outfile hex.exr hex.pbrt
..\..\bin\exrtotiff.exe -scale 1 hex.exr hex.tiff