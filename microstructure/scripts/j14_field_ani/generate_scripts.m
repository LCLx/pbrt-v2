% Tao Du
% taodu@csail.mit.edu
% Jan 12, 2017
name = 'j14_field_ani';
root_folder = ['\\SCALAR\share_topo\rendering2\', name, '\'];
n = 100;

for i = 0 : n
    mkdir(num2str(i));
    system([name, '_no_fine_intf_flag.bat ', num2str(i)]);
    copyfile([name, '.pbrt'], num2str(i));
end