..\..\x64\Release\microstructure.exe \\SCALAR\share_topo\rendering2\j12_elip2_3d\0\lattice \\SCALAR\share_topo\rendering2\j12_elip2_3d\0\ale_dis \\SCALAR\share_topo\rendering2\j12_elip2_3d\0\material hex_mesh.pbrt
..\..\..\bin\pbrt.exe --ncores 12 --outfile hex.exr hex.pbrt
..\..\..\bin\exrtotiff.exe -scale 1 hex.exr hex.tiff