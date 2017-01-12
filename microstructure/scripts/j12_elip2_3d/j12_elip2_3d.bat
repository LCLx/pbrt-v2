set frame_index=%1
set root_folder=\\SCALAR\share_topo\rendering2\j12_elip2_3d\%frame_index%
..\..\x64\Release\microstructure.exe %root_folder%\lattice %root_folder%\ale_dis %root_folder%\material %root_folder%\sing_points j12_elip2_3d_mesh.pbrt
..\..\..\bin\pbrt.exe --ncores 12 --outfile j12_elip2_3d_%frame_index%.exr j12_elip2_3d.pbrt
..\..\..\bin\exrtotiff.exe -scale 4 j12_elip2_3d_%frame_index%.exr j12_elip2_3d_%frame_index%.tiff