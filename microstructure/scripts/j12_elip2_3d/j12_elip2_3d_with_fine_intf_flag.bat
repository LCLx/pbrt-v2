set frame_index=%1
set root_folder=\\SCALAR\share_topo\rendering2\j12_elip2_3d\%frame_index%
..\..\x64\Release\microstructure.exe j12_elip2_3d_mesh.pbrt %root_folder% lattice ale_dis material NULL sing_points fine_intf_flags
..\..\..\bin\pbrt.exe --ncores 12 --outfile j12_elip2_3d_00%frame_index%_with_fine_intf_flags.exr j12_elip2_3d.pbrt
..\..\..\bin\exrtotiff -scale 2 j12_elip2_3d_00%frame_index%_with_fine_intf_flags.exr j12_elip2_3d_00%frame_index%_with_fine_intf_flags.tiff