function [] = render(pbrt_script_folder)
% Tao Du
% taodu@csail.mit.edu
% Jan 27, 2017
%
% Generate images from the given pbrt scripts.
% pbrt_script_folder: the folder that contains the pbrt scripts. The
% structure should be as follows:
% pbrt_script_folder/
%   0/render.pbrt
%   1/render.pbrt
%   2/render.pbrt
% ...
% We will create a new folder named 'render' under pbrt_script_folder, and
% put images under that folder.
idx = 0;
pbrt_exe_location = 'D:\library\pbrt-v2\bin\';

render_folder = fullfile(pbrt_script_folder, 'render');
if ~exist(render_folder, 'dir')
  mkdir(render_folder);
end

while true
  pbrt_idx_folder = fullfile(pbrt_script_folder, num2str(idx));
  if ~exist(pbrt_idx_folder, 'dir')
    break;
  end
  pbrt_file_name = fullfile(pbrt_idx_folder, 'render.pbrt');
  exr_file_name = fullfile(render_folder, [num2str(idx, '%03d'), '.exr']);
  tiff_file_name = fullfile(render_folder, [num2str(idx, '%03d'), '.tiff']);
  system(sprintf('%s --ncores 12 --outfile %s %s\n', ...
    fullfile(pbrt_exe_location, 'pbrt.exe'), exr_file_name, pbrt_file_name));
  system(sprintf('%s %s %s\n', ...
    fullfile(pbrt_exe_location, 'exrtotiff.exe'), exr_file_name, tiff_file_name));
  idx = idx + 1;
end

end