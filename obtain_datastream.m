function [DataStream, means] = obtain_datastream(Gloms, sizePatches, downsample, colorspace, numFeatures, removeMeans, include_segments)
%obtain_datastream extracts all viable patches from glom images

%% Case 1: No initial segmentation
if ~isfield(Gloms{1}, 'Segmentation')
    %% Convert to DataStream of patches
    fprintf('Extracting patches...\n')
    DataStream = zeros(numFeatures, 0);
    numGloms = length(Gloms);
    for glom_ticker = 1:numGloms
        fprintf('\tGlom %i...\t', glom_ticker)
        DataStream = [DataStream, extract_patches(Gloms{glom_ticker}, sizePatches, downsample, 'sliding', colorspace, 'all', include_segments)];
        fprintf('Complete!\n')
    end
    numDataPts = size(DataStream, 2);
    fprintf('Number of data points: %i\n', numDataPts);
    fprintf('===============================================\n')

    %% Remove means
    means = zeros(size(DataStream));
    if removeMeans
        means(0*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream(0*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
        means(1*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream(1*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
        means(2*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
            repmat(mean(DataStream(2*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);

        DataStream = DataStream - means;
    end
    
else 
%% Case 2: Initial segmentation exists
    %% Convert to DataStream of patches
    fprintf('Extracting patches...\n')
    DataStream = cell(2,1);
    DataStream{1} = zeros(numFeatures, 0);
    DataStream{2} = zeros(numFeatures, 0);
    numGloms = length(Gloms);
    for glom_ticker = 1:numGloms
        fprintf('\tGlom %i...\t', glom_ticker)
        DataStream_temp = extract_patches(Gloms{glom_ticker}, sizePatches, downsample, 'sliding', colorspace, 'all', include_segments);
        DataStream{1} = [DataStream{1}, DataStream_temp{1}];
        DataStream{2} = [DataStream{2}, DataStream_temp{2}];
        fprintf('Complete!\n')
    end
    numDataPts = [size(DataStream{1}, 2), size(DataStream{2}, 2)];
    clear DataStream_temp
    fprintf('Number of data points: %i (Spectral method)\n', numDataPts(1));
    fprintf('Number of data points: %i (NeuralNet method)\n', numDataPts(2));
    fprintf('===============================================\n')

    %% Remove means
    means = cell(2,1);
    means{1} = zeros(size(DataStream{1}));
    means{2} = zeros(size(DataStream{2}));
    if removeMeans
        for ticker = 1:2
            currDataStream = DataStream{ticker};
            currMeans = means{ticker};
            
            currMeans(0*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
                repmat(mean(currDataStream(0*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
            currMeans(1*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
                repmat(mean(currDataStream(1*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
            currMeans(2*(numFeatures/3)+(1:(numFeatures/3)),:) = ...
                repmat(mean(currDataStream(2*(numFeatures/3)+(1:(numFeatures/3)),:)), numFeatures/3, 1);
            
            means{ticker} = currMeans;
            DataStream{ticker} = DataStream{ticker} - means{ticker};
        end
    end
end

end