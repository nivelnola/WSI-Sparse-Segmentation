%% Glomerulus Segmentation - Stage 3
%   Compute a neural network-based segmentation, and compare.
%   Alon S. Levin

%% Prepare the environment
clc, clear, nnet.guis.closeAllViews(), %close all
format compact
rmpath(genpath('D:\Program Files\MATLAB\Custom Packages\sporco-m0.0.9'))

fprintf('#######################################\n')
fprintf('## Glomerulus Segmentation - Stage 3 ##\n')
fprintf('#######################################\n')
fprintf('\tAuthor:\t\tAlon S. Levin\n')
fprintf('\tDate:\t\tAugust 20, 2021\n')
fprintf('\tVersion:\t2.5\n')
fprintf('===============================================\n')

%% Restore Data
fprintf('Restoring data... ')
[file,path] = uigetfile('F:\PAS_Dataset\Segmentation_Thesis\Glom_2_*.mat', 'Select a file.');
load([path, file])
clear path file
fprintf('Complete!\n') 
fprintf('===============================================\n')

%% Settings
% Data settings
numGloms     = length(Gloms);                     % Number of glom images to use for clustering
sizePatches  = Gloms{1}.Segmentation.sizePatches; % Size of patches, side length
downsample   = Gloms{1}.Segmentation.downsample;  % Downsample factor
                        % NOTE: True number of pixels in a patch is
                        % downsample * sizePatches^2
numFeatures = 3*sizePatches^2;

% Model settings
numDicts    = 3;                        % Number of dictionaries
numAtomsTot = numFeatures*(numDicts+0); % Total number of atoms across all dictionaries
colorspace  = "HSV";                    % Color space to use ['RGB', 'HSV', 'CIELAB']
colormodel  = "Concatenation";          % Color model to use ['Concatenation', 'Quaternion']
removeMeans = false;                    % Remove means, boolean

% SOM Settings
dimensions = [1, numDicts];
coverSteps = 100;
initNeighbor = numDicts;
topologyFcn = 'hextop';
distanceFcn = 'linkdist';
SoM_full = selforgmap(dimensions,coverSteps,initNeighbor,topologyFcn,distanceFcn);

%% Obtain DataStream
[DataStream, means] = obtain_datastream(Gloms, sizePatches, downsample, colorspace, numFeatures, removeMeans, [1,2,3]);
DataStream = DataStream{1};
means = means{1};
%% Train SoM

fprintf('Training neural net... ')
SoM_full = train(SoM_full,DataStream);
fprintf('Complete!\n')
fprintf('===============================================\n')

%% Compare to previous results
% Rebuild
for glom_ticker = 1:length(Gloms)
    fprintf('Segmenting Glom %i...\n', glom_ticker)
    Glom_test = Gloms{glom_ticker};
    fprintf('Complete!\n')

    fprintf('\tConverting glom to datastream... ')
    [DataStream_test, numPatches_test, indices_test] = extract_patches(Glom_test, sizePatches, downsample, 'sliding', colorspace, 'any', [1,2,3]);
    DataStream_test = DataStream_test{1};
    indices_test = indices_test{1};
    means_test = zeros(size(DataStream_test));
    if removeMeans
        means_test(0*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream_test(0*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
        means_test(1*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream_test(1*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
        means_test(2*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream_test(2*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);

        DataStream_test = DataStream_test - means_test;
    end
    fprintf('Complete!\n')

    fprintf('\tAssigning classes to patches... \n')
    y = SoM_full(DataStream_test);
    classes = vec2ind(y);
    
    temp_assignments = zeros(1, numPatches_test);
    temp_assignments(indices_test) = classes;
    ClusterAssignments_test = temp_assignments;
    tabulate(ClusterAssignments_test);
    clear temp_assignments ic counts

    fprintf('\tAssigning classes to pixels...\n')
    Mask  = imresize(Glom_test.Mask, round(size(Glom_test.Mask)/downsample));
    pixel_grid = reshape(1:numel(Mask), size(Mask));
    pixel_grid_patches = im2col(pixel_grid, [sizePatches, sizePatches], 'sliding');
    clusters_cols = repmat(ClusterAssignments_test, sizePatches.^2, 1);
    pixel_clusters = zeros(size(pixel_grid));
    for pixel_ticker = 1:numel(pixel_grid)
        if(Mask(pixel_ticker))
            pixel_clusters(pixel_ticker) = mode(clusters_cols(pixel_grid_patches == pixel_ticker));
        end
    end
    fprintf('Complete!\n')
    
    fprintf('\tNeural net segment sizes:\n')
    t = tabulate(pixel_clusters(:));
    t(2:numDicts+1, :) = sortrows(t(2:numDicts+1, :), 3);
    pixel_clusters = pixel_clusters + 100;
    for ticker = 1:numDicts
        pixel_clusters(pixel_clusters == t(ticker+1, 1)+100) = ticker;
    end
    pixel_clusters(pixel_clusters == 100) = 0;
    tabulate(pixel_clusters(:));
    Glom_test.Segmentation.Clusters_som = pixel_clusters;

    fprintf('\tResizing segmentations to original size...')
    Glom_test.Segmentation.Clusters_sp_full = ...
        imresize(Glom_test.Segmentation.Clusters_sp, size(Glom_test.Image, 1:2), 'method', 'nearest');
    Glom_test.Segmentation.Clusters_nn_full = ...
        imresize(Glom_test.Segmentation.Clusters_nn, size(Glom_test.Image, 1:2), 'method', 'nearest');
    Glom_test.Segmentation.Clusters_som_full = ...
        imresize(Glom_test.Segmentation.Clusters_som, size(Glom_test.Image, 1:2), 'method', 'nearest');
    
    fprintf('\tStoring segmentations...')
    Mask = Glom_test.Mask;
    Glom_test.comparison_sp = (sum(Glom_test.Segmentation.Clusters_som_full == Glom_test.Segmentation.Clusters_sp_full, 'all') - sum(~Mask, 'all'))/sum(Mask, 'all');
    Glom_test.comparison_nn = (sum(Glom_test.Segmentation.Clusters_som_full == Glom_test.Segmentation.Clusters_nn_full, 'all') - sum(~Mask, 'all'))/sum(Mask, 'all');
    Gloms{glom_ticker} = Glom_test;
    fprintf('Complete!\n')
    
    fprintf('\tVisualizing segmentations... ')
    Image = im2double(Glom_test.Image);
    Image_masked = Image .* Mask;    
    figure
        subplot(3,3,1)
            image(Image_masked)
            title('Original Image')
        subplot(3,3,2)
            image(labeloverlay(Image_masked, Glom_test.Segmentation.Clusters_sp_full))
            title('Spectral Clustering (Overlay)')
        subplot(3,3,3)
            imagesc(Glom_test.Segmentation.Clusters_sp_full)
            cmap = colormap(parula(length(Glom_test.Segmentation.Dictionaries)+1));
            cbh = colorbar ;
            cbh.Ticks = linspace(0, length(Glom_test.Segmentation.Dictionaries)+1, length(Glom_test.Segmentation.Dictionaries)+2);
            cbh.TickLabels = num2cell(0:length(Glom_test.Segmentation.Dictionaries)+1);
            title('Spectral Clustering (Segmentation)')
        subplot(3,3,4)
            imshowpair(Glom_test.Segmentation.Clusters_sp_full, Glom_test.Segmentation.Clusters_nn_full);
            title('Spec. Clustering (Pink) vs. Hybrid Method (Green)')
        subplot(3,3,5)
            image(labeloverlay(Image_masked, Glom_test.Segmentation.Clusters_nn_full))
            title('Hybrid Method (Overlay)')
        subplot(3,3,6)
            imagesc(Glom_test.Segmentation.Clusters_nn_full)
            cmap = colormap(parula(length(Glom_test.Segmentation.Dictionaries)+2));
            cbh = colorbar ;
            cbh.Ticks = linspace(0, length(Glom_test.Segmentation.Dictionaries)+2, length(Glom_test.Segmentation.Dictionaries)+3);
            cbh.TickLabels = num2cell(0:length(Glom_test.Segmentation.Dictionaries)+2);
            title('Hybrid Method (Segmentation)')
        subplot(3,3,7)
            imshowpair(Glom_test.Segmentation.Clusters_som_full, Glom_test.Segmentation.Clusters_nn_full);
            title('Neural Net (Pink) vs. Hybrid Method (Green)')
        subplot(3,3,8)
            image(labeloverlay(Image_masked, Glom_test.Segmentation.Clusters_som_full))
            title('Neural Net (Overlay)')
        subplot(3,3,9)
            imagesc(Glom_test.Segmentation.Clusters_som_full)
            cmap = colormap(parula(length(Glom_test.Segmentation.Dictionaries)+2));
            cbh = colorbar ;
            cbh.Ticks = linspace(0, length(Glom_test.Segmentation.Dictionaries)+2, length(Glom_test.Segmentation.Dictionaries)+3);
            cbh.TickLabels = num2cell(0:length(Glom_test.Segmentation.Dictionaries)+2);
            title('Neural Net (Segmentation)')

    figure
        for ticker = 1:(numDicts)
            subplot(3, numDicts, ticker)
                imshowpair(Image_masked, Glom_test.Segmentation.Clusters_sp_full == ticker, 'falsecolor')
                title(['Spectral Clustering - Class ', num2str(ticker), ' (Pink)']);
            subplot(3, numDicts, ticker+3)
                imshowpair(Image_masked, Glom_test.Segmentation.Clusters_nn_full == ticker, 'falsecolor')
                title(['Hybrid Method - Class ', num2str(ticker), ' (Pink)']);
            subplot(3, numDicts, ticker+6)
                imshowpair(Image_masked, Glom_test.Segmentation.Clusters_som_full == ticker, 'falsecolor')
                title(['Neural Net - Class ', num2str(ticker), ' (Pink)']);
        end
        
    figure
        subplot(1,2,1)
            imshowpair(Glom_test.Segmentation.Clusters_som_full, Glom_test.Segmentation.Clusters_sp_full);
            title('Neural Net vs. Spectral Clustering Method')
        subplot(1,2,2)
            imshowpair(Glom_test.Segmentation.Clusters_som_full, Glom_test.Segmentation.Clusters_nn_full);
            title('Neural Net vs. SoM Clustering Method')
            
    fprintf('Complete!\n')
    fprintf('-----------------------------------------------\n')
    
end

confirmsave = input("Would you like to save these results? Y/N [Y]: ", 's');
    switch confirmsave
        case {'N', 'n'}
            warning('Results not saved manually. Do not clear before saving.')
        case {'', 'Y', 'y'}
            fprintf('Saving... ')
            timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss'));
            save(['F:\PAS_Dataset\Segmentation_Thesis\Glom_3_', timestamp, '.mat'], 'Gloms', '-v7.3')
            fprintf('Complete!\n')
    end
    
fprintf('===============================================\n')
fprintf('######### Stage 3 processing complete #########\n')
fprintf('===============================================\n\n')
