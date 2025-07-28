%% Glomerulus Segmentation - Stage 1
%   Segment a glomerulus into two classes.
%   Alon S. Levin

%% Prepare the environment
clc, nnet.guis.closeAllViews(), clear, close all
format compact

fprintf('#######################################\n')
fprintf('## Glomerulus Segmentation - Stage 1 ##\n')
fprintf('#######################################\n')
fprintf('\tAuthor:\t\tAlon S. Levin\n')
fprintf('\tDate:\t\tAugust 27, 2021\n')
fprintf('\tVersion:\t4.2\n')
fprintf('===============================================\n')

%% Settings
% Data settings
numGloms       = 1;     % Number of glom images to use for clustering
sizePatches    = 8;     % Size of patches, side length
downsample     = 2;     % Downsample factor
                        % NOTE: True number of pixels in a patch is
                        % (downsample * sizePatches)^2
numFeatures = 3*sizePatches^2;

% Model settings
numDicts    = 2;                        % Number of dictionaries
numAtomsTot = numFeatures*(numDicts+0); % Total number of atoms across all dictionaries
colorspace  = "HSV";                    % Color space to use ['RGB', 'HSV', 'CIELAB']
colormodel  = "Concatenation";          % Color model to use ['Concatenation', 'Quaternion']
removeMeans = false;                    % Remove means, boolean

% Dictionary Learning settings
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
fprintf('GOAL: Extract the lining of the glomeruli.\n')
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
fprintf('Dictionary Learning Settings:\n')
fprintf('\tlambda:\t\t\t%f\n', lambda)
fprintf('\tnumOptimIter:\t%i\n', numOptimIter)
fprintf('\tFor more settings, see set_opts.m\n')
fprintf('Self-Organizing Map Settings:\n')
fprintf('\tdimensions:\t\t[%i, %i]\n', dimensions)
fprintf('\tcoverSteps:\t\t%i\n', coverSteps)
fprintf('\tinitNeighbor:\t%i\n', initNeighbor)
fprintf('\ttopologyFcn:\t%s\n', topologyFcn)
fprintf('\tdistanceFcn:\t%s\n', distanceFcn)
fprintf('===============================================\n')


%% Obtain DataStream
[DataStream, Gloms, means] = Obtain_Data(numGloms, sizePatches, downsample, colorspace, numFeatures, removeMeans, 1);

%% Dictionary Initialization
addpath(genpath('D:\Program Files\MATLAB\Custom Packages\sporco-m0.0.9'))
D_0 = Initialize_Dictionary(DataStream, numAtomsTot, colormodel, lambda);
fprintf('===============================================\n')

% Visualize several atoms and full dictionary
fprintf('Visualizing dictionary... ')
figure('Name', 'Sample Atoms from General Dictionary')
    plot_rand_patches(data2patches(D_0, sizePatches), [5,5], colorspace)

figure('Name', 'Initial Dictionary Visualization')
visualize_dict_layers(D_0);

fprintf('Complete!\n')
fprintf('===============================================\n')


%% Clustering
[Dictionary, A_0, SoM, cluster_ids_sp, cluster_ids_nn] = ...
    Cluster_Atoms(DataStream, D_0, numDicts, colormodel, lambda, SoM, "both");
fprintf('===============================================\n')

% Visualize reconstructions
fprintf('Reconstructing patches... ')
Reconstruction = D_0 * A_0;
fprintf('Complete!\n')

figure('Name', 'Sample Reconstructions')
    plot_reconstruction_comp(DataStream+means, Reconstruction+means, [5 5], colorspace)
fprintf('===============================================\n')

%% Check clustering
fprintf('Spectral clustering distribution:\n')
tabulate(cluster_ids_sp)
fprintf('Neural net clustering distribution:\n')
tabulate(cluster_ids_nn)
fprintf('===============================================\n')

%% Iterative Optimization
fprintf('Performing optimization...\n')
fprintf('\tSpectral clustering...\n')
Dictionary_opt(1, :) = Optimize_Dictionary(Dictionary(1, :), DataStream, lambda, colormodel, numOptimIter);
fprintf('\tHybrid model...\n')
Dictionary_opt(2, :) = Optimize_Dictionary(Dictionary(2, :), DataStream, lambda, colormodel, numOptimIter);
fprintf('Complete!\n')
fprintf('===============================================\n')


%% Segment Gloms, Store
fprintf('Applying Dictionaries to Gloms...\n')
fprintf('-----------------------------------------------\n')

for glom_ticker = 1:length(Gloms)
    fprintf('Segmenting Glom %i... ', glom_ticker)
    Glom_test = Gloms{glom_ticker};
    fprintf('Complete!\n')

    fprintf('\tConverting glom to datastream... ')
    [DataStream_test, numPatches_test, indices_test] = extract_patches(Glom_test, sizePatches, downsample, 'sliding', colorspace, 'any', 1);
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

    figure('Name', 'Test Glom Patchified')
        plot_image_patchgrid(Glom_test, sizePatches, downsample, 'g');
        title('Test Glom')

    fprintf('Complete!\n')

    fprintf('\tReconstructing patches... \n')
    A_test = cell(size(Dictionary));
    for dict_ticker = 1:numDicts
        fprintf('\t\tDictionary %i... ', dict_ticker)
        A_test{1, dict_ticker} = sparse(bpdn(Dictionary{1, dict_ticker}, DataStream_test, lambda, set_opts("test_code")));
        A_test{2, dict_ticker} = sparse(bpdn(Dictionary{2, dict_ticker}, DataStream_test, lambda, set_opts("test_code")));
        fprintf('Complete!\n')
    end

    fprintf('\tCalculating classification metric... \n')
    R_hat_test = nan(numDicts, size(DataStream_test, 2), 2);
    for dict_ticker = 1:numDicts
        fprintf('\t\tDictionary %i... ', dict_ticker)
        R_hat_test(dict_ticker, :, 1) = calculate_reconstruction_metric(DataStream_test, Dictionary{1, dict_ticker}, A_test{1, dict_ticker}, lambda);
        R_hat_test(dict_ticker, :, 2) = calculate_reconstruction_metric(DataStream_test, Dictionary{2, dict_ticker}, A_test{2, dict_ticker}, lambda);
        fprintf('Complete!\n')
    end

    fprintf('\tAssigning classes to patches... \n')
    [~,ClusterAssignments_sp] = min(R_hat_test(:,:,1));
    [~,ClusterAssignments_nn] = min(R_hat_test(:,:,2));
    
    temp_assignments = zeros(1, numPatches_test);
    temp_assignments(indices_test) = ClusterAssignments_sp;
    ClusterAssignments_sp = temp_assignments;
    
    temp_assignments = zeros(1, numPatches_test);
    temp_assignments(indices_test) = ClusterAssignments_nn;
    ClusterAssignments_nn = temp_assignments;
    
    fprintf('\t\tSpectral clustering:\n')
    tabulate(ClusterAssignments_sp);
    fprintf('\t\tSoM clustering:\n')
    tabulate(ClusterAssignments_nn);
    clear temp_assignments ic counts

    fprintf('\tAssigning classes to pixels... ')
    Mask = imresize(Glom_test.Mask, round(size(Glom_test.Mask)/downsample));
    pixel_grid = reshape(1:numel(Mask), size(Mask));
    pixel_grid_patches = im2col(pixel_grid, [sizePatches, sizePatches], 'sliding');
    clusters_cols_sp = repmat(ClusterAssignments_sp, sizePatches.^2, 1);
    clusters_cols_nn = repmat(ClusterAssignments_nn, sizePatches.^2, 1);
    pixel_clusters_sp = zeros(size(pixel_grid));
    pixel_clusters_nn = zeros(size(pixel_grid));
    for pixel_ticker = 1:numel(pixel_grid)
        if(Mask(pixel_ticker))
            pixel_clusters_sp(pixel_ticker) = mode(clusters_cols_sp(pixel_grid_patches == pixel_ticker));
            pixel_clusters_nn(pixel_ticker) = mode(clusters_cols_nn(pixel_grid_patches == pixel_ticker));
        end
    end    
    fprintf('Complete!\n')
        
    fprintf('\tStoring segmentations...')
    Segmentation.downsample = downsample;
    Segmentation.sizePatches = sizePatches;
    Segmentation.Clusters_sp = pixel_clusters_sp;
    Segmentation.Clusters_nn = pixel_clusters_nn;
    Segmentation.Dictionaries = Dictionary;
    
    if(sum(Segmentation.Clusters_sp == 1, 'all') > sum(Segmentation.Clusters_sp == 2, 'all'))
        Segmentation = swap_segmentation_labels(Segmentation, 1, 2, "SP");
    end
    if(sum(Segmentation.Clusters_nn == 1, 'all') > sum(Segmentation.Clusters_nn == 2, 'all'))
        Segmentation = swap_segmentation_labels(Segmentation, 1, 2, "NN");
    end
    
    Gloms{glom_ticker}.Segmentation = Segmentation;
    Glom_test = Gloms{glom_ticker};
    fprintf('Complete!\n')

    fprintf('\tVisualizing segmentations... ')
    Image = im2double(imresize(Glom_test.Image, round(size(Glom_test.Image, 1:2)/downsample)));
    Image_masked = Image .* Mask;    
    figure
        subplot(2,3,1)
            image(Image_masked)
            title('Original Image')
        subplot(2,3,2)
            image(labeloverlay(Image_masked, Segmentation.Clusters_sp))
            title('Spectral Clustering (Overlay)')
        subplot(2,3,3)
            imagesc(Segmentation.Clusters_sp)
            cmap = colormap(parula(length(Segmentation.Dictionaries)+1));
            cbh = colorbar;
            cbh.Ticks = linspace(0, length(Segmentation.Dictionaries)+1, length(Segmentation.Dictionaries)+2);
            cbh.TickLabels = num2cell(0:length(Segmentation.Dictionaries)+1);
            title('Spectral Clustering (Segmentation)')
        subplot(2,3,4)
            imshowpair(Segmentation.Clusters_sp, Segmentation.Clusters_nn);
            title('Spec. Clustering (Pink) vs. Hybrid Method (Green)')
        subplot(2,3,5)
            image(labeloverlay(Image_masked, Segmentation.Clusters_nn))
            title('Hybrid Method (Overlay)')
        subplot(2,3,6)
            imagesc(Segmentation.Clusters_nn)
            cmap = colormap(parula(length(Segmentation.Dictionaries)+1));
            cbh = colorbar;
            cbh.Ticks = linspace(0, length(Segmentation.Dictionaries)+1, length(Segmentation.Dictionaries)+2);
            cbh.TickLabels = num2cell(0:length(Segmentation.Dictionaries)+1);
            title('Hybrid Method (Segmentation)')
    figure
        for ticker = 1:numDicts
            subplot(2, numDicts, ticker)
                imshowpair(Image_masked, Segmentation.Clusters_sp == ticker, 'falsecolor')
                title(['Spectral Clustering - Class ', num2str(ticker), ' (Pink)']);
            subplot(2, numDicts, ticker+2)
                imshowpair(Image_masked, Segmentation.Clusters_nn == ticker, 'falsecolor')
                title(['Hybrid Method - Class ', num2str(ticker), ' (Pink)']);
        end
        
    figure
        plot_dictionary(data2patches(D_0, sizePatches), 'HSV')
        title('Original Dictionary')

    figure
        subplot(2,2,1)
            plot_dictionary(data2patches(Segmentation.Dictionaries{1,1}, sizePatches), 'HSV')
            title('Class 1 (SP)')
        subplot(2,2,2)
            plot_dictionary(data2patches(Segmentation.Dictionaries{1,2}, sizePatches), 'HSV')
            title('Class 2 (SP)')
        subplot(2,2,3)
            plot_dictionary(data2patches(Segmentation.Dictionaries{2,1}, sizePatches), 'HSV')
            title('Class 1 (SoM)')
        subplot(2,2,4)
            plot_dictionary(data2patches(Segmentation.Dictionaries{2,2}, sizePatches), 'HSV')
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
            save(['F:\PAS_Dataset\Segmentation_Thesis\Glom_1_', timestamp, '.mat'], 'Gloms', '-v7.3')
            fprintf('Complete!\n')
    end
    
fprintf('===============================================\n')
fprintf('######### Stage 1 processing complete #########\n')
fprintf('===============================================\n\n')