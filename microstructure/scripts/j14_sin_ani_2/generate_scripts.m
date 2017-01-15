% Tao Du
% taodu@csail.mit.edu
% Jan 12, 2017
root_folder = '\\SCALAR\share_topo\rendering2\j14_sin_ani_2\';
n = 165;

for i = 0 : n
    mkdir(num2str(i));
    system(['j14_sin_ani_2_no_fine_intf_flag.bat ', num2str(i)]);
    copyfile('j14_sin_ani_2.pbrt', num2str(i));
end