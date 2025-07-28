function [DataStream, numPatches, indices] = extract_patches(Glom, sizePatches, downsampleFactor, method, colorspace, fill_status, include_segments)
%extract_patches Given a Glom object, extracts patches (in DataStream
%format) of given size.
%   Inputs:
%       Structure:          Struct representing a Glom
%       sizePatches:        Scalar, side length of square patch
%       downsampleFactor:   Downsample factor for patches
%       method:             'distinct' or 'sliding'
%       colorspace:         'RGB', 'HSV', or 'LAB'
%       fill_status:        'any', 'all'
%       include_segments:   Segment labels to include
%   Outputs:
%       DataStream:     2-D matrix
%             Note: Use data2patches to see actual patches!

%% Resize mask and image
Mask  = imresize(Glom.Mask, round(size(Glom.Mask)/downsampleFactor));
Image = im2double(imresize(Glom.Image, round(size(Glom.Image, 1:2)/downsampleFactor)));

%% Choose patches to ignore
%   If the patch's mask vector contains zeros, ignore
if ~isfield(Glom, 'Segmentation')
    mask_datastream = im2col(Mask, [sizePatches, sizePatches], method);
    numPatches = size(mask_datastream, 2);      % Equivalent to prod(size(Mask) - sizePatches + 1)
    switch fill_status
        case "any"
            indices = any(mask_datastream == 1);
        case "all"
            indices = all(mask_datastream == 1);
    end
else
    mask_datastream_sp = im2col(Glom.Segmentation.Clusters_sp, [sizePatches, sizePatches], method);
    mask_datastream_nn = im2col(Glom.Segmentation.Clusters_nn, [sizePatches, sizePatches], method);
    numPatches = size(mask_datastream_sp, 2);
    switch fill_status
        case "any"
            indices_sp = any(ismember(mask_datastream_sp, include_segments));
            indices_nn = any(ismember(mask_datastream_nn, include_segments));
        case "all"
            indices_sp = all(ismember(mask_datastream_sp, include_segments));
            indices_nn = all(ismember(mask_datastream_nn, include_segments));
    end
end

%% Convert colorspaces
switch colorspace
    case 'RGB'
    case 'HSV'
        Image = rgb2hsv(Image);
    case 'CIELAB'
        Image = rgb2lab(Image);
end
        
%% Convert the image to datastream
C1 = im2col(Image(:,:,1), [sizePatches, sizePatches], method);
if (colorspace ~= "CIELAB")
    C2 = im2col(Image(:,:,2), [sizePatches, sizePatches], method);
    C3 = im2col(Image(:,:,3), [sizePatches, sizePatches], method);
else
    C2 = zeros(size(C1));
    C3 = zeros(size(C2));
end

DataStream_temp = [C1; C2; C3];

%% Remove bad patches
if ~isfield(Glom, 'Segmentation')
    DataStream = DataStream_temp(:, indices);
else
    DataStream = cell(2,1);
    DataStream{1} = DataStream_temp(:, indices_sp);
    DataStream{2} = DataStream_temp(:, indices_nn);
    indices = {indices_sp; indices_nn};
end


