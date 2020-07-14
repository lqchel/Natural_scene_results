% this function generate grey masks that cuts a square image into
% equal-size patches. img_size = width of image in pixels, start_x and start_y = upper
% left corner of the area of the image you want to cut patch from (for entire
% image, start_x = 1, start_y = 1). patch_size = width of patches.
% num_patches = number of patches you want to get from the image; has to be square number.


function [mask] = Draw_Mask(img_size, start_x,start_y,patch_size, num_patch)

grey_mask = uint8(128+zeros(img_size,img_size,3));
grey_mask = rgb2gray(grey_mask);

for i = 1:sqrt(num_patch)

    x_begin = start_x +(i-1).*patch_size;
    x_end = i.*patch_size;
    
    for a = 1:sqrt(num_patch)
        tp = zeros(img_size,img_size) + 1;
        y_begin = start_y + (a-1).*patch_size;
        y_end = a.*patch_size;
        tp(x_begin:x_end,y_begin:y_end) = 0;
        save_path = ['mask_' num2str(a+(i-1).*sqrt(num_patch)) '.png'];
        imwrite(grey_mask,save_path,'Alpha',tp);
        clear tp
    end
end
end

