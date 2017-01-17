% Tao Du
% taodu@csail.mit.edu
% Jan 12, 2017
name = 'j15_bridge3_ale_64_3_final';

root_folder = ['\\SCALAR\share_topo\rendering2\', name, '\'];
n = 100;

for i = 0 : n
    mkdir(num2str(i));
    system([name, '.bat ', num2str(i)]);
    copyfile([name, '.pbrt'], num2str(i));
end