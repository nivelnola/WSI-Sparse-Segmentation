function plot_rand_patches(D, dim, colorspace)
%plot_dictionary Given a dictionary of patches D (arranged
%[x,y,z,p_i]), plots a random selection of them in the form dim

% Prepare window arrangement
rows = dim(1);
cols = dim(2);

% Colored data
if(ndims(D) == 4)
    % Choose patches randomly
    patch_nums = randi(size(D, 4), dim);

    for patch_ticker = 1:prod(dim)
        % Choose current patch
        switch colorspace
            case 'RGB'
                currPatch = rescale(D(:,:,:,patch_nums(patch_ticker)));
            case 'HSV'
                currPatch = hsv2rgb(rescale(D(:,:,:,patch_nums(patch_ticker))));
            case 'CIELAB'
                currPatch = lab2rgb(rescale(D(:,:,:,patch_nums(patch_ticker))));
        end

        % Create the plot
        subplot(rows, cols, patch_ticker)
        imshow(currPatch)
        title("Patch " + patch_nums(patch_ticker));    
    end
end

% Grayscale data
if(ndims(D) == 3)
    % Choose patches randomly
    patch_nums = randi(size(D, 3), dim);

    for patch_ticker = 1:prod(dim)
        % Choose current patch
        currPatch = rescale(D(:,:,patch_nums(patch_ticker)));

        % Create the plot
        subplot(rows, cols, patch_ticker)
        image(currPatch)
        title("Patch " + patch_nums(patch_ticker));    
    end
end

end

