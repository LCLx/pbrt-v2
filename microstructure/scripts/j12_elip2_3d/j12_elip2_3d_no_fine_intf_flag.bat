set frame_index=%1
set root_folder=\\SCALAR\share_topo\rendering2\j12_elip2_3d\%frame_index%
..\..\x64\Release\microstructure.exe %frame_index%\j12_elip2_3d_mesh.pbrt %root_folder% lattice ale_dis material lag_intf_points NULL NULL