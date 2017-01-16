% Tao Du
% taodu@csail.mit.edu
% Jan 12, 2017
root_folder = '\\SCALAR\share_topo\rendering2\j6_dia3_ale\';
n = 75;

for i = 0 : n
    mkdir(num2str(i));
    system(['j6_dia3_ale.bat ', num2str(i)]);
    copyfile('j6_dia3_ale.pbrt', num2str(i));
end