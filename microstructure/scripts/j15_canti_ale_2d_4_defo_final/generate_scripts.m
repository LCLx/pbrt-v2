% Tao Du
% taodu@csail.mit.edu
% Jan 12, 2017
name = 'j15_canti_ale_2d_4_defo_final';
root_folder = ['\\SCALAR\share_topo\rendering2\', name, '\'];
n = 78;

for i = 78 : 136
    mkdir(num2str(i));
    system([name, '.bat ', num2str(i)]);
    copyfile([name, '.pbrt'], num2str(i));
end