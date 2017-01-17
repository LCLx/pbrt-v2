set name=j16_ray_4_f509_print_final

..\..\..\bin\pbrt.exe --ncores 12 --outfile %name%.exr %name%.pbrt
..\..\..\bin\exrtotiff -scale 2 %name%.exr %name%.tiff