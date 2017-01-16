set frame_index=%1
set root_folder=\\SCALAR\share_topo\rendering2\j6_dia3_ale\%frame_index%
..\..\x64\Release\microstructure.exe mesh.pbrt %root_folder% lattice ale_dis material NULL NULL fine_intf_flags NULL NULL density 0
..\..\..\bin\pbrt.exe --ncores 12 --outfile %1_no_surrounding_cells.exr j6_dia3_ale.pbrt
..\..\..\bin\exrtotiff %1_no_surrounding_cells.exr %1_no_surrounding_cells.tiff