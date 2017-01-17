set frame_index=%1
set name=j14_field_ani
set root_folder=\\SCALAR\share_topo\rendering2\%name%\%frame_index%
..\..\x64\Release\microstructure.exe %frame_index%\%name%_mesh.pbrt %root_folder% lattice ale_dis material lag_intf_points NULL NULL