% Tao Du
% taodu@csail.mit.edu
% Jan 12, 2017
root_folder = '\\SCALAR\share_topo\rendering2\j15_bridge3_uni_64\';
n = 100;

for i = 0 : n
    mkdir(num2str(i));
    system(['j15_bridge3_uni_64.bat ', num2str(i)]);
    copyfile('j15_bridge3_uni_64.pbrt', num2str(i));
end