set name=j16_stretch_bar_1_f270_print

..\..\..\bin\pbrt.exe --ncores 12 --outfile %name%.exr %name%.pbrt
..\..\..\bin\exrtotiff -scale 3 %name%.exr %name%.tiff