set name=j16_curved_grip_2d_64_f470_print_6_final

..\..\..\bin\pbrt.exe --ncores 12 --outfile %name%.exr %name%.pbrt
..\..\..\bin\exrtotiff -scale 2 %name%.exr %name%.tiff