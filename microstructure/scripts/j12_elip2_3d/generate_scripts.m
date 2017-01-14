% Tao Du
% taodu@csail.mit.edu
% Jan 12, 2017
root_folder = '\\SCALAR\share_topo\rendering2\j12_elip2_3d\';
n = 200;

for i = 0 : 200
    mkdir(num2str(i));
    system(['j12_elip2_3d_no_fine_intf_flag.bat ', num2str(i)]);
    copyfile('j12_elip2_3d.pbrt', num2str(i));
end