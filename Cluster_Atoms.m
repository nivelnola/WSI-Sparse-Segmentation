function [Dictionary, A_0, SoM, cluster_ids_sp, cluster_ids_nn] = Cluster_Atoms(DataStream, D_0, numDicts, colormodel, lambda, SoM, cluster_type)
% Given a learned total dictionary, cluster the atoms into separate
% dictionaries

%% Step 1 - Calculate sparse representations
fprintf('Calculating reconstructions...\n')
tic
switch colormodel
    case "Concatenation"
        A_0 = bpdn(D_0, DataStream, lambda, set_opts("init_code_concat"));
        A_0_sparse = sparse(A_0);
    case "Quaternion"
        A_0 = QOMP(concat2quaternion(D_0), concat2quaternion(DataStream), set_opts("init_code_quater").L);
        A_0_sparse = sparse(quaternion2concat(A_0));
end
t = toc;
fprintf('\tComplete!\n')
fprintf('Runtime: %f sec\n', t)

fprintf('Sparsity: %2.2f%%\n', 100*(1-nnz(A_0)/numel(A_0)))
fprintf('Memory comparison:\n\n')
whos A*
A_0 = A_0_sparse;

%% Step 2 - Calculate similarity matrix
fprintf('Calculating similarity matrix... ')
S = abs(A_0 * A_0');
fprintf('Complete!\n')

%% Step 3 - Cluster on atoms
fprintf('Clustering atoms...\n')
switch cluster_type
    case "both"
        fprintf('\tSpectral clustering... ')
        cluster_ids_sp = spectralcluster(S, numDicts, 'Distance', 'precomputed', 'LaplacianNormalization', 'symmetric');
        fprintf('Complete!\n')
        
        fprintf('\tSoM Clustering... ')
        rmpath(genpath('D:\Program Files\MATLAB\Custom Packages\sporco-m0.0.9'))
        SoM = train(SoM, full(S));
        cluster_ids_nn = vec2ind(SoM(full(S)));
        addpath(genpath('D:\Program Files\MATLAB\Custom Packages\sporco-m0.0.9'))
        
    case "sp"
        fprintf('\tSpectral clustering... ')
        cluster_ids_sp = spectralcluster(S, numDicts, 'Distance', 'precomputed', 'LaplacianNormalization', 'symmetric');
        
    case "nn"
        fprintf('\tSoM Clustering... ')
        rmpath(genpath('D:\Program Files\MATLAB\Custom Packages\sporco-m0.0.9'))
        SoM = train(SoM, full(S));
        cluster_ids_nn = vec2ind(SoM(full(S)))';
        addpath(genpath('D:\Program Files\MATLAB\Custom Packages\sporco-m0.0.9'))
end
fprintf('Complete!\n')

%% Step 4 - Separate dictionaries
fprintf('Separating dictionaries... ')
switch cluster_type
    case "both"
        Dictionary = cell(2,numDicts);
        for dict_ticker = 1:numDicts
            Dictionary{1,dict_ticker} = D_0(:, cluster_ids_sp == dict_ticker, :);
            Dictionary{2,dict_ticker} = D_0(:, cluster_ids_nn == dict_ticker, :);
        end
    case "sp"
        Dictionary = cell(1,numDicts);
        for dict_ticker = 1:numDicts
            Dictionary{dict_ticker} = D_0(:, cluster_ids_sp == dict_ticker, :);
        end
    case "nn"
        Dictionary = cell(1,numDicts);
        for dict_ticker = 1:numDicts
            Dictionary{dict_ticker} = D_0(:, cluster_ids_nn == dict_ticker, :);
        end
end
fprintf('Complete!\n')

%% Step 6 - Visualizations
switch cluster_type
    case "sp"
        cluster_ids_nn = NaN;
        SoM = NaN;
    case "nn"
        cluster_ids_sp = NaN;
end

%% Visualization

cluster_nn_1 = A_0(cluster_ids_nn == 1, :);
cluster_nn_2 = A_0(cluster_ids_nn == 2, :);
A_clustered_nn = [cluster_nn_1; cluster_nn_2];
S_clustered_nn = abs(A_clustered_nn) * abs(A_clustered_nn)';

cluster_sp_1 = A_0(cluster_ids_sp == 1, :);
cluster_sp_2 = A_0(cluster_ids_sp == 2, :);
A_clustered_sp = [cluster_sp_1; cluster_sp_2];
S_clustered_sp = abs(A_clustered_sp) * abs(A_clustered_sp)';

figure('Name', 'Unclustered Similarity Matrix')
    imagesc(full(abs(A_0) * abs(A_0)'))
    colormap('turbo')
    title('Similarity Matrix, Unclustered')
    colorbar
    axis square;

figure('Name', 'Clustered Similarity Matrix')
switch cluster_type
    case "both"
        subplot(1,2,1)
            imagesc(S_clustered_sp)
            title('Similarity Matrix, Spectrally Clustered')
            colormap('turbo')
            colorbar   
            axis square;
        subplot(1,2,2)
            imagesc(S_clustered_nn)
            title('Similarity Matrix, SoM Clustered')
            colormap('turbo')
            colorbar
            axis square;
    case "sp"
        imagesc(S_clustered_sp)
        title('Similarity Matrix, Spectrally Clustered')
        colormap('turbo')
        colorbar   
        axis square;
    case "nn"
        imagesc(S_clustered_nn)
        title('Similarity Matrix, SoM Clustered')
        colormap('turbo')
        colorbar   
        axis square;
end

if cluster_type == "sp" || cluster_type == "both"
    graph_unclustered = graph(S);
    graph_cluster_1   = subgraph(graph_unclustered, find(cluster_ids_sp == 1));
    graph_cluster_2   = subgraph(graph_unclustered, find(cluster_ids_sp == 2));
    figure('Name', 'Spectral Clustering Graph')
        subplot(1,2,1)
            p1 = plot(graph_unclustered, 'NodeLabel',  {}, ...
                                         'NodeColor', 	'#006400', ...
                                         'EdgeColor', [0.4660 0.6740 0.1880], ...
                                         'Marker', '*', ...
                                         'MarkerSize', 7);
            title('Original Graph')
        subplot(1,2,2)
            hold on
            plot(graph_cluster_1, 'XData', p1.XData(find(cluster_ids_sp == 1)), ...
                                  'YData', p1.YData(find(cluster_ids_sp == 1)), ...
                                  'NodeLabel',  {}, ...
                                  'EdgeColor', [0.3010 0.7450 0.9330], ...
                                  'NodeColor', 'b', ...
                                  'Marker', '*', ...
                                  'MarkerSize', 7);
            plot(graph_cluster_2, 'XData', p1.XData(find(cluster_ids_sp == 2)), ...
                                  'YData', p1.YData(find(cluster_ids_sp == 2)), ...
                                  'NodeLabel',  {}, ...
                                  'EdgeColor', [0.9290 0.6940 0.1250], ...
                                  'NodeColor', 'r', ...
                                  'Marker', '*', ...
                                  'MarkerSize', 7);
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            title('Clustered Subgraphs')
            legend({'Cluster 1', 'Cluster 2'})
end