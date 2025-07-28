function plot_image_patchgrid(Glom, sizePatches, downsampleFactor, color)
%plot_image_patchgrid produces a plot with a grid overlay representing distinct patches

%% Resize mask and image
Mask  = imresize(Glom.Mask, round(size(Glom.Mask)/sqrt(downsampleFactor)));
Image = im2double(imresize(Glom.Image, round(size(Glom.Image, 1:2)/sqrt(downsampleFactor))));
%% Calculate rows and columns
[rows, cols, ~] = size(Image);
%% Produce plot
imshow(Image .* Mask)
hold on
for row = 0:sizePatches:rows
    line([1, cols], [row, row], 'Color', color, 'LineWidth', 2);
end
for col = 0:sizePatches:cols
    line([col, col], [1, rows], 'Color', color, 'LineWidth', 2);
end

end

