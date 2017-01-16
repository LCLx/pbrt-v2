writerObj = VideoWriter('j12_elip2_3d_linear_interp.avi');
writerObj.FrameRate = 30;
% open the video writer
open(writerObj);
% write the frames to the video
img_num = 12;
for u = 0 : img_num
    % convert the image to a frame
    img = imread(['rendering/j12_elip2_3d_', num2str(u, '%03d'), '.tiff']);
    frame = im2frame(img(:, :, 1 : 3));
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);