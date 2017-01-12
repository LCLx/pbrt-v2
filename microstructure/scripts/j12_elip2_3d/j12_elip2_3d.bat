set frame_index=%1
set root_folder=\\SCALAR\share_topo\rendering2\j12_elip2_3d\%frame_index%
..\..\x64\Release\microstructure.exe %root_folder%\lattice %root_folder%\ale_dis %root_folder%\material %root_folder%\lag_intf_points %frame_index%\j12_elip2_3d_mesh.pbrt