set frame_index=%1
set root_folder=\\SCALAR\share_topo\rendering2\j15_bridge3_ale_64_3_defo\%frame_index%
..\..\x64\Release\microstructure.exe %frame_index%\mesh.pbrt %root_folder% lattice ale_dis material NULL sing_points fine_intf_flags