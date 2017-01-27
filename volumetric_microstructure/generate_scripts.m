function [ ] = generate_scripts( density_folder, pbrt_folder )
% Tao Du
% taodu@csail.mit.edu
% Jan 27, 2017
%
% Read density files from density folder, and write the corresponding pbrt
% scripts into pbrt folder.
% density_folder: folder that contains binary files. The structure should
% look like this:
% density_folder/
%   0/density
%   1/density
%   2/density
% ...
% pbrt_folder: a folder that contains the pbrt scripts we are going to
% generate. The structure looks like this:
% pbrt_folder/
%   0/geometry.pbrt
%   1/geometry.pbrt
%   2/geometry.pbrt
% ...
idx = 0;
converter_location = 'D:\library\pbrt-v2\volumetric_microstructure\main\x64\Release\main.exe';

while true
  density_file_name = fullfile(density_folder, num2str(idx), 'density');
  if ~exist(density_file_name, 'file')
    break;
  end
  pbrt_idx_folder = fullfile(pbrt_folder, num2str(idx));
  if ~exist(pbrt_idx_folder, 'dir')
    mkdir(pbrt_idx_folder);
  end
  pbrt_file_name = fullfile(pbrt_idx_folder, 'geometry.pbrt');
  system(sprintf('%s -m mixed -d %s -p %s -t 0.5', converter_location, density_file_name, pbrt_file_name));
  copyfile('render.pbrt', pbrt_idx_folder);
  idx = idx + 1;
end

end
