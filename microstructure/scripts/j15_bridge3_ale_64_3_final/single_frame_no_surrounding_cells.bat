set frame_index=%1
set threshold=%2
set name=j15_bridge3_ale_64_3_final
set root_folder=\\SCALAR\share_topo\rendering2\%name%\%frame_index%
..\..\x64\Release\microstructure.exe mesh.pbrt %root_folder% lattice ale_dis material NULL NULL fine_intf_flags NULL NULL density 0 %threshold%
..\..\..\bin\pbrt.exe --ncores 12 --outfile %1_no_surrounding_cells.exr %name%.pbrt
..\..\..\bin\exrtotiff %1_no_surrounding_cells.exr %1_no_surrounding_cells.tiff