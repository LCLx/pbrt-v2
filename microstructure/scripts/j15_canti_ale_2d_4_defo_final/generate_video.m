name = 'j15_canti_ale_2d_4_defo_final';

writerObj = VideoWriter([name, '.avi']);
writerObj.FrameRate = 30;
% open the video writer
open(writerObj);
% write the frames to the video
img_num = (length(dir('rendering')) - 2) / 2;
step = 3;
for u = 0 : (img_num - 1)
 % convert the image to a frame
 img = imread(['rendering/', name, '_', num2str(step * u, '%03d'), '.tiff']);
 frame = im2frame(img(:, :, 1 : 3));
 writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);