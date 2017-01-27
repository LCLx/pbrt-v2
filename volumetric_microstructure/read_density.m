function [ density ] = read_density( file_name )
% Tao Du
% taodu@csail.mit.edu
% Jan 27, 2017
%
% Read density data from a binary file.
% file_name: name of the binary file. The file contains 64 x 64 x 64 double
% values, in the following order:
% id_x id_y id_z
% 0 0 0
% 0 0 1
% 0 0 2
% ...
% 0 0 63
% 0 1 0
% 0 1 1
% 0 1 2
% ...
% 0 1 63
% 0 2 0
% ...
% ...
% 1 0 0
% 1 0 1
% ...
% density: a 64 x 64 x 64 3-D tensor. density(id_x, id_y, id_z) gives the
% corresponding values at cell (id_x, id_y, id_z).

file_id = fopen(file_name);
nx = fread(file_id, 1, 'int');
ny = fread(file_id, 1, 'int');
nz = fread(file_id, 1, 'int');
density = fread(file_id, inf, 'double');
density = permute(reshape(density, [nz, ny, nx]), [3, 2, 1]);
fclose(file_id);

end

