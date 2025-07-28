%% Glomerulus Segmentation - Stage 2
%   Using a partially segmented glomerulus, segment the out-of-class
%   segment.
%   Alon S. Levin

%% Prepare the environment
clc, clear, %close all, nnet.guis.closeAllViews()
format compact

fprintf('#######################################\n')
fprintf('## Glomerulus Segmentation - Stage 2 ##\n')
fprintf('#######################################\n')
fprintf('\tAuthor:\t\tAlon S. Levin\n')
fprintf('\tDate:\t\tAugust 27, 2021\n')
fprintf('\tVersion:\t3.2\n')
fprintf('===============================================\n')

%% Restore Data
fprintf('Restoring data... ')
[file,path] = uigetfile('F:\PAS_Dataset\Segmentation_Thesis\Glom_1_*.mat', 'Select a file.');
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
                        % (downsample * sizePatches)^2
numFeatures = 3*sizePatches^2;

% Model settings
numDicts    = 2;                        % Number of dictionaries
numAtomsTot = numFeatures*(numDicts); % Total number of atoms across all dictionaries
colorspace  = "HSV";                    % Color space to use ['RGB', 'HSV', 'CIELAB']
colormodel  = "Concatenation";          % Color model to use ['Concatenation', 'Quaternion']
removeMeans = false;                    % Remove means, boolean

% Learning settings
lambda       = 1.2 / sqrt(numFeatures);
numOptimIter = 5;

% Self-Organizing Map settings
rmpath(genpath('D:\Program Files\MATLAB\Custom Packages\sporco-m0.0.9'))
dimensions = [1, numDicts];
coverSteps = 100;
initNeighbor = numDicts;
topologyFcn = 'hextop';
distanceFcn = 'linkdist';
SoM = selforgmap(dimensions,coverSteps,initNeighbor,topologyFcn,distanceFcn);

% Report settings
fprintf('GOAL: Segment the internals of the glomeruli.\n')
fprintf('-----------------------------------------------\n')
fprintf('Data Settings:\n')
fprintf('\tnumGloms:\t\t%i\n', numGloms)
fprintf('\tsizePatches:\t%i\n', sizePatches)
fprintf('\tdownsample:\t\t%i\n', downsample)
fprintf('Model Settings:\n')
fprintf('\tnumDicts:\t\t%i\n', numDicts)
fprintf('\tnumAtomsTot:\t%i\n', numAtomsTot)
fprintf('\tcolorspace:\t\t%s\n', colorspace)
fprintf('\tcolormodel:\t\t%s\n', colormodel)
fprintf('\tremoveMeans:\t%s\n', mat2str(removeMeans))
fprintf('Learning Settings:\n')
fprintf('\tlambda:\t\t\t%f\n', lambda)
fprintf('\tnumOptimIter:\t%i\n', numOptimIter)
fprintf('For more settings, see set_opts.m\n')
fprintf('===============================================\n')

%% Obtain DataStream
[DataStream, means] = obtain_datastream(Gloms, sizePatches, downsample, colorspace, numFeatures, removeMeans, 2);

%% Dictionary Initialization
addpath(genpath('D:\Program Files\MATLAB\Custom Packages\sporco-m0.0.9'))
fprintf('Spectral-Clustered Dictionary\n')
D_0_sp = Initialize_Dictionary(DataStream{1}, numAtomsTot, colormodel, lambda);
fprintf('-----------------------------------------------\n')
fprintf('Neural Net-Clustered Dictionary\n')
D_0_nn = Initialize_Dictionary(DataStream{2}, numAtomsTot, colormodel, lambda);
fprintf('===============================================\n')

%% Clustering
fprintf('Spectral-Clustered Dictionary\n')
[Dictionary_sp, A_0_sp, ~, cluster_ids_sp, ~] = Cluster_Atoms(DataStream{1}, D_0_sp, numDicts, colormodel, lambda, SoM, "sp");
fprintf('-----------------------------------------------\n')
fprintf('Neural Net-Clustered Dictionary\n')
[Dictionary_nn, A_0_nn, SoM, ~, cluster_ids_nn] = Cluster_Atoms(DataStream{2}, D_0_nn, numDicts, colormodel, lambda, SoM, "nn");
fprintf('===============================================\n')

%% Check clustering
fprintf('Clustering distribution:\n')
fprintf('\tSpectral Clustering distribution:\n')
tabulate(cluster_ids_sp)
fprintf('\tNeural Net Clustering distribution:\n')
tabulate(cluster_ids_nn)
fprintf('===============================================\n')

%% Iterative Optimization
fprintf('Performing optimization...\n')
Dictionary_sp = Optimize_Dictionary(Dictionary_sp, DataStream{1}, lambda, colormodel, numOptimIter);
Dictionary_nn = Optimize_Dictionary(Dictionary_nn, DataStream{2}, lambda, colormodel, numOptimIter);
Dictionary = [Dictionary_sp; Dictionary_nn];
fprintf('Complete!\n')
fprintf('===============================================\n')

%% Segment Gloms, Store
fprintf('Applying Dictionaries to Gloms...\n')
fprintf('-----------------------------------------------\n')

for glom_ticker = 1:length(Gloms)
    fprintf('Segmenting Glom %i...\n', glom_ticker)
    Glom_test = Gloms{glom_ticker};
    fprintf('Complete!\n')

    fprintf('\tConverting glom to datastream... ')
    [DataStream_test, numPatches_test, indices_test] = extract_patches(Glom_test, sizePatches, downsample, 'sliding', colorspace, 'any', 2);
    DataStream_test_sp = DataStream_test{1};
    DataStream_test_nn = DataStream_test{2};
    clear DataStream_test
    
    means_test_sp = zeros(size(DataStream_test_sp));
    means_test_nn = zeros(size(DataStream_test_nn));
    if removeMeans
        means_test_sp(0*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream_test_sp(0*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
        means_test_sp(1*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream_test_sp(1*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
        means_test_sp(2*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream_test_sp(2*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
        DataStream_test_sp = DataStream_test_sp - means_test;
        
        means_test_nn(0*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream_test_nn(0*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
        means_test_nn(1*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream_test_nn(1*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
        means_test_nn(2*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream_test_nn(2*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
        DataStream_test_nn = DataStream_test_nn - means_test;
    end

    fprintf('Complete!\n')

    fprintf('\tReconstructing patches... \n')
    A_test = cell(size(Dictionary));
    for dict_ticker = 1:numDicts
        fprintf('\t\tDictionary %i... ', dict_ticker)
        A_test{1, dict_ticker} = sparse(bpdn(Dictionary{1, dict_ticker}, DataStream_test_sp, lambda, set_opts("test_code")));
        A_test{2, dict_ticker} = sparse(bpdn(Dictionary{2, dict_ticker}, DataStream_test_nn, lambda, set_opts("test_code")));
        fprintf('Complete!\n')
    end

    fprintf('\tCalculating classification metric... \n')
    R_hat_test_sp = nan(numDicts, size(DataStream_test_sp, 2));
    R_hat_test_nn = nan(numDicts, size(DataStream_test_nn, 2));
    for dict_ticker = 1:numDicts
        fprintf('\t\tDictionary %i... ', dict_ticker)
        R_hat_test_sp(dict_ticker, :) = calculate_reconstruction_metric(DataStream_test_sp, Dictionary{1, dict_ticker}, A_test{1, dict_ticker}, lambda);
        R_hat_test_nn(dict_ticker, :) = calculate_reconstruction_metric(DataStream_test_nn, Dictionary{2, dict_ticker}, A_test{2, dict_ticker}, lambda);
        fprintf('Complete!\n')
    end

    fprintf('\tAssigning classes to patches... \n')
    [~,ClusterAssignments_sp] = min(R_hat_test_sp);
    [~,ClusterAssignments_nn] = min(R_hat_test_nn);
    
    ClusterAssignments_sp = ClusterAssignments_sp + 1;
    ClusterAssignments_nn = ClusterAssignments_nn + 1;
    
    temp_assignments = zeros(1, numPatches_test);
    temp_assignments(indices_test{1}) = ClusterAssignments_sp;
    ClusterAssignments_sp = temp_assignments;
    
    temp_assignments = zeros(1, numPatches_test);
    temp_assignments(indices_test{2}) = ClusterAssignments_nn;
    ClusterAssignments_nn = temp_assignments;
    
    fprintf('\t\tSpectral clustering:\n')
    tabulate(ClusterAssignments_sp);
    fprintf('\t\tSoM clustering:\n')
    tabulate(ClusterAssignments_nn);
    clear temp_assignments ic counts

    fprintf('\tAssigning classes to pixels...')
    clusters_cols_sp = repmat(ClusterAssignments_sp, sizePatches.^2, 1);
    clusters_cols_nn = repmat(ClusterAssignments_nn, sizePatches.^2, 1);
    
    Mask_sp = (Glom_test.Segmentation.Clusters_sp == 2);
    Mask_nn = (Glom_test.Segmentation.Clusters_nn == 2);
    
    pixel_grid = reshape(1:numel(Mask_sp), size(Mask_sp));
    pixel_grid_patches = im2col(pixel_grid, [sizePatches, sizePatches], 'sliding');
    segmentation_sp = zeros(size(pixel_grid));
    segmentation_nn = zeros(size(pixel_grid));
    
    for pixel_ticker = 1:numel(pixel_grid)
        if(Mask_sp(pixel_ticker))
            segmentation_sp(pixel_ticker) = mode(clusters_cols_sp(pixel_grid_patches == pixel_ticker));
        end
        if(Mask_nn(pixel_ticker))
            segmentation_nn(pixel_ticker) = mode(clusters_cols_nn(pixel_grid_patches == pixel_ticker));
        end
    end
    fprintf('Complete!\n')
    
    fprintf('\tStoring segmentations... ')
    Glom_test.Segmentation.Clusters_sp = ...
        Glom_test.Segmentation.Clusters_sp.*(Glom_test.Segmentation.Clusters_sp <= 1) ...
      + segmentation_sp.*(Glom_test.Segmentation.Clusters_sp == 2);
    Glom_test.Segmentation.Clusters_nn = ...
        Glom_test.Segmentation.Clusters_nn.*(Glom_test.Segmentation.Clusters_nn <= 1) ...
      + segmentation_nn.*(Glom_test.Segmentation.Clusters_nn == 2);
    Glom_test.Segmentation.Dictionaries_2 = Dictionary;
    Gloms{glom_ticker}.Segmentation = Glom_test.Segmentation;
    fprintf('Complete!\n')
    
    fprintf('\tVisualizing segmentations... ')
    Image = im2double(imresize(Glom_test.Image, round(size(Glom_test.Image, 1:2)/downsample)));
    Image_masked = Image .* im2double(imresize(Glom_test.Mask, round(size(Glom_test.Mask, 1:2)/downsample)));    
    figure
        subplot(2,3,1)
            image(Image_masked)
            title('Original Image')
        subplot(2,3,2)
            image(labeloverlay(Image_masked, Glom_test.Segmentation.Clusters_sp))
            title('Spectral Clustering (Overlay)')
        subplot(2,3,3)
            imagesc(Glom_test.Segmentation.Clusters_sp)
            cmap = colormap(parula(length(Glom_test.Segmentation.Dictionaries)+2));
            cbh = colorbar ;
            cbh.Ticks = linspace(0, length(Glom_test.Segmentation.Dictionaries)+2, length(Glom_test.Segmentation.Dictionaries)+3);
            cbh.TickLabels = num2cell(0:length(Glom_test.Segmentation.Dictionaries)+2);
            title('Hybrid Method (Segmentation)')
            title('Spectral Clustering (Segmentation)')
        subplot(2,3,4)
            imshowpair(Glom_test.Segmentation.Clusters_sp, Glom_test.Segmentation.Clusters_nn);
            title('Spec. Clustering (Pink) vs. Hybrid Method (Green)')
        subplot(2,3,5)
            image(labeloverlay(Image_masked, Glom_test.Segmentation.Clusters_nn))
            title('Hybrid Method (Overlay)')
        subplot(2,3,6)
            imagesc(Glom_test.Segmentation.Clusters_nn)
            cmap = colormap(parula(length(Glom_test.Segmentation.Dictionaries)+2));
            cbh = colorbar ;
            cbh.Ticks = linspace(0, length(Glom_test.Segmentation.Dictionaries)+2, length(Glom_test.Segmentation.Dictionaries)+3);
            cbh.TickLabels = num2cell(0:length(Glom_test.Segmentation.Dictionaries)+2);
            title('Hybrid Method (Segmentation)')

    figure
        for ticker = 1:(numDicts+1)
            subplot(2, numDicts+1, ticker)
                imshowpair(Image_masked, Glom_test.Segmentation.Clusters_sp == ticker, 'falsecolor')
                title(['Spectral Clustering - Class ', num2str(ticker), ' (Pink)']);
            subplot(2, numDicts+1, ticker+3)
                imshowpair(Image_masked, Glom_test.Segmentation.Clusters_nn == ticker, 'falsecolor')
                title(['Hybrid Method - Class ', num2str(ticker), ' (Pink)']);
        end
    
    figure
        subplot(2,1,1)
            plot_dictionary(data2patches(D_0_sp, sizePatches), 'HSV')
            title('Original Dictionary (SP)')
        subplot(2,1,2)
            plot_dictionary(data2patches(D_0_nn, sizePatches), 'HSV')
            title('Original Dictionary (SoM)')

    figure
        subplot(2,2,1)
            plot_dictionary(data2patches(Dictionary_sp{1}, sizePatches), 'HSV')
            title('Class 1 (SP)')
        subplot(2,2,2)
            plot_dictionary(data2patches(Dictionary_sp{2}, sizePatches), 'HSV')
            title('Class 2 (SP)')
        subplot(2,2,3)
            plot_dictionary(data2patches(Dictionary_nn{1}, sizePatches), 'HSV')
            title('Class 1 (SoM)')
        subplot(2,2,4)
            plot_dictionary(data2patches(Dictionary_nn{2}, sizePatches), 'HSV')
            title('Class 2 (SoM)')
        
        
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
            save(['F:\PAS_Dataset\Segmentation_Thesis\Glom_2_', timestamp, '.mat'], 'Gloms', '-v7.3')
            fprintf('Complete!\n')
    end
    
fprintf('===============================================\n')
fprintf('######### Stage 2 processing complete #########\n')
fprintf('===============================================\n\n')