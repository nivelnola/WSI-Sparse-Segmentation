function plot_dictionary(D, colorspace)
%plot_dictionary Given a dictionary of RGB patches D (arranged
%[x,y,z,p_i]), displays all patches in a figure

% 3D data
if(ndims(D) == 4)
    % Patch-dimensions
    x_dim = size(D,2);
    y_dim = size(D,1);
    
    % Calculate subplot arrangement
    [dims, ~] = numSubplots(size(D, 4));
    rows = dims(1);
    cols = dims(2);
    all_blocks = nan(1+rows*(y_dim+1), 1+cols*(x_dim+1), 3);
    
    for patch_ticker = 1:size(D, 4)
        % Choose current patch
        switch colorspace
            case 'RGB'
                currPatch = rescale(D(:,:,:,patch_ticker));
            case 'HSV'
                currPatch = hsv2rgb(rescale(D(:,:,:,patch_ticker)));
            case 'CIELAB'
                currPatch = lab2rgb(rescale(D(:,:,:,patch_ticker)));
        end
        
        % Add patch to image
        start_row = (y_dim+1)*(mod(patch_ticker-1, rows)) + 2;
        start_col = (x_dim+1)*(floor((patch_ticker-1) / rows)) + 2;
        all_blocks(start_row:start_row+y_dim-1, start_col:start_col+x_dim-1, :) = currPatch;
    end
    
end

% 2D data
if(ndims(D) == 3)
    % Patch-dimensions
    x_dim = size(D,2);
    y_dim = size(D,1);
    
    % Calculate subplot arrangement
    [dims, ~] = numSubplots(size(D, 4));
    rows = dims(1);
    cols = dims(2);
    all_blocks = nan(1+rows*(y_dim+1), 1+cols*(x_dim+1), 3);

    for patch_ticker = 1:size(D, 3)
        % Choose current patch
        currPatch = rescale(D(:,:,patch_ticker));
        
        % Add patch to image
        start_row = (y_dim+1)*(mod(patch_ticker-1, rows)) + 2;
        start_col = (x_dim+1)*(floor((patch_ticker-1) / rows)) + 2;
        all_blocks(start_row:start_row+y_dim-1, start_col:start_col+x_dim-1, :) = currPatch;
    end
end

%% Create the plot
imagesc(all_blocks);
hold on
for row = 1:(y_dim+1):(1+rows*(y_dim+1))
    line([0, 1+cols*(x_dim+1)], [row, row], 'Color', 'g', 'LineWidth', 2);
end
for col = 1:(x_dim+1):(1+cols*(x_dim+1))
    line([col, col], [0, 1+rows*(y_dim+1)], 'Color', 'g', 'LineWidth', 2);
end