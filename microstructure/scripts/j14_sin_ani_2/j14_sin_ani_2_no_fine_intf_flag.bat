set frame_index=%1
set root_folder=\\SCALAR\share_topo\rendering2\j14_sin_ani_2\%frame_index%
..\..\x64\Release\microstructure.exe %frame_index%\j14_sin_ani_2_mesh.pbrt %root_folder% lattice ale_dis material lag_intf_points NULL NULL