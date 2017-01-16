set frame_index=%1
set threshold=%2
set root_folder=\\SCALAR\share_topo\rendering2\j15_canti_ale_2d_4\%frame_index%
..\..\x64\Release\microstructure.exe mesh.pbrt %root_folder% lattice ale_dis material NULL NULL fine_intf_flags NULL NULL density 0 %threshold%
..\..\..\bin\pbrt.exe --ncores 12 --outfile %1_no_surrounding_cells.exr j15_canti_ale_2d_4.pbrt
..\..\..\bin\exrtotiff %1_no_surrounding_cells.exr %1_no_surrounding_cells.tiff