% Tao Du
% taodu@csail.mit.edu
% Jan 12, 2017
name = 'j6_dia3_ale_defo';
root_folder = ['\\SCALAR\share_topo\rendering2\', name, '\'];
n = 111;

for i = 0 : n
    mkdir(num2str(i));
    system([name, '.bat ', num2str(i)]);
    copyfile([name, '.pbrt'], num2str(i));
end