% Tao Du
% taodu@csail.mit.edu
% Jan 12, 2017
root_folder = '\\SCALAR\share_topo\rendering2\j12_elip2_3d\';
n = numel(dir(root_folder)) - 2;

for i = 0 : (n - 1)
    mkdir(num2str(i));
    system(['j12_elip2_3d.bat ', num2str(i)]);
    copyfile('j12_elip2_3d.pbrt', num2str(i));
end