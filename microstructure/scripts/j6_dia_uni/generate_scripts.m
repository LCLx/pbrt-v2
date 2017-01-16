% Tao Du
% taodu@csail.mit.edu
% Jan 12, 2017
root_folder = '\\SCALAR\share_topo\rendering2\j6_dia_uni\';
n = 16;

for i = 0 : n
    mkdir(num2str(i));
    system(['j6_dia_uni.bat ', num2str(i)]);
    copyfile('j6_dia_uni.pbrt', num2str(i));
end