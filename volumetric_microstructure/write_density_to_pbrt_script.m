function [ ] = write_density_to_pbrt_script( density, pbrt_file_name )
% Tao Du
% taodu@csail.mit.edu
% Jan 27, 2017
%
% Given a 64 x 64 x 64 density cube, write it to a pbrt script that defines
% a triangle mesh, which contains 64 x 64 x 64 cells, and the rgb values of
% the material is specified by the density.
% density: a 64 x 64 x 64 3D tensor, output of read_density.
% pbrt_file_name: the name of the pbrt file.

file_id = fopen(pbrt_file_name, 'w');

fprintf(file_id, 'AttributeBegin\n');
% Loop over all the elements.
dx = 1.0 / 64;
points = [0, 0, 0
  0, 0, dx
  0, dx, 0
  0, dx, dx
  dx, 0, 0
  dx, 0, dx
  dx, dx, 0
  dx, dx, dx];
for i = 1 : 64
  for j = 1 : 64
    for k = 1 : 64
      d = density(i, j, k);
      if d == 0
        continue;
      end
      string_format = ['Material "matte" "rgb Kd" [%d %d %d]\n', ...
        'Shape "trianglemesh"\n', ...
        '"integer indices" [\n', ...
        '4 6 7 4 7 5 0 3 2 0 1 3 1 3 7 1 7 5 2 0 6 6 0 4 2 3 7 2 7 6 0 5 1 0 4 5\n', ...
        ']\n'];
      fprintf(file_id, string_format, d, d, d);
      corner = ([i, j, k] - 1) * dx;
      fprintf(file_id, '"point P" [\n');
      fprintf(file_id, '%d %d %d\n', points + corner);
      fprintf(file_id, ']\n');
    end
  end
end
fprintf(file_id, 'AttributeEnd\n');

fclose(file_id);

end