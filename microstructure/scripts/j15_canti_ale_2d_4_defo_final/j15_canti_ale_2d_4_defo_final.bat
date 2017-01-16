set frame_index=%1
set root_folder=\\SCALAR\share_topo\rendering2\j15_canti_ale_2d_4_defo_final\%frame_index%
..\..\x64\Release\microstructure.exe %frame_index%\mesh.pbrt %root_folder% lattice ale_dis material NULL sing_points fine_intf_flags