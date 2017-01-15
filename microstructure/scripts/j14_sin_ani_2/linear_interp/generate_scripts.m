% Tao Du
% taodu@csail.mit.edu
% Jan 13, 2017
root_folder = '\\SCALAR\share_topo\rendering2\j14_sin_ani_2\';
n = 60; % 2s.
for i = 0 : n
    t = i / n;
    mkdir(num2str(i));
    system(['j14_sin_ani_2_linear_interp.bat ', num2str(i), ' ', num2str(t)]);
    copyfile('..\j14_sin_ani_2.pbrt', num2str(i));
%    system(['..\..\..\..\bin\pbrt.exe --ncores 12 --outfile rendering\', num2str(i, '%03d'), '.exr ', ...
%        num2str(i), '\j12_elip2_3d.pbrt']);
%    system(['..\..\..\..\bin\exrtotiff.exe -scale 2 rendering\', num2str(i, '%03d'), '.exr rendering\', ...
%        num2str(i, '%03d'), '.tiff']);
end