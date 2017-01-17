..\..\..\bin\pbrt.exe --ncores 12 --outfile deformed_grid.exr render_grid.pbrt
..\..\..\bin\exrtotiff -scale 2 deformed_grid.exr deformed_grid.tiff