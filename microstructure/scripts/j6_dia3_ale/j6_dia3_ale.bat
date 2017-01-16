set frame_index=%1
set root_folder=\\SCALAR\share_topo\rendering2\j6_dia3_ale\%frame_index%
..\..\x64\Release\microstructure.exe %frame_index%\mesh.pbrt %root_folder% lattice ale_dis material NULL NULL NULL f_points psi_Ds density