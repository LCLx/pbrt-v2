set frame_index=%1
set root_folder=\\SCALAR\share_topo\rendering2\j14_sin_ani_2\%frame_index%
..\..\x64\Release\microstructure.exe j14_sin_ani_2_mesh.pbrt %root_folder% lattice ale_dis material NULL NULL fine_intf_flags
..\..\..\bin\pbrt.exe --ncores 12 --outfile j14_sin_ani_2_00%frame_index%_with_fine_intf_flags.exr j14_sin_ani_2.pbrt
..\..\..\bin\exrtotiff -scale 2 j14_sin_ani_2_00%frame_index%_with_fine_intf_flags.exr j14_sin_ani_2_00%frame_index%_with_fine_intf_flags.tiff