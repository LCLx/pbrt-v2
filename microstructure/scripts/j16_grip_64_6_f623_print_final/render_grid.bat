..\..\..\bin\pbrt.exe --ncores 12 --outfile grid.exr render_grid.pbrt
..\..\..\bin\exrtotiff -scale 2 grid.exr grid.tiff