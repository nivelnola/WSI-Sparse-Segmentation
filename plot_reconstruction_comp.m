function plot_reconstruction_comp(DataStream, Reconstruction, dim, colorspace)
%plot_reconstruction_comp Shows a comparison of original data and
%reconstructions

% Prepare window arrangement
rows = dim(1);
%cols = dim(2);

% Obtain sizePatches
sizePatches = sqrt(size(DataStream, 1) / 3);

% Choose patches randomly
    patch_nums = randi(size(DataStream, 2), [rows, 1]); %dim);

    for patch_ticker = 1:rows %prod(dim)
        % Choose a patch number
        patch_number = patch_nums(patch_ticker);
        
        % Choose current patch
        switch colorspace
            case 'RGB'
                currPatch_orig  = rescale(data2patches(DataStream(:,patch_number), sizePatches));
                currPatch_recon = rescale(data2patches(Reconstruction(:,patch_number), sizePatches));
            case 'HSV'
                currPatch_orig  = hsv2rgb(rescale(data2patches(DataStream(:,patch_number), sizePatches)));
                currPatch_recon = hsv2rgb(rescale(data2patches(Reconstruction(:,patch_number), sizePatches)));
            case 'CIELAB'
                currPatch_orig  = lab2rgb(rescale(data2patches(DataStream(:,patch_number), sizePatches)));
                currPatch_recon = lab2rgb(rescale(data2patches(Reconstruction(:,patch_number), sizePatches)));
        end
        
        % Create the plot
        subplot(rows, 2, 2*patch_ticker-1)
        imshow(imshowpair(currPatch_orig, currPatch_recon, 'montage').CData)
        title('Original vs. Reconstruction')
        subplot(rows, 2, 2*patch_ticker)
        imshow(imshowpair(currPatch_orig, currPatch_recon).CData)
        title('Difference')
        axis equal;
    end

end

