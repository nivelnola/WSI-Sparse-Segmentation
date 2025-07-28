function [Segmentation] = swap_segmentation_labels(Segmentation, L1, L2, method)
%swap_segmentation_labels Swaps the label of the segmentation, and changes
%the order of the corresponding dictionaries

%% Replace labels in clusters
switch method
    case 'SP'
        Segmentation.Clusters_sp(Segmentation.Clusters_sp == L1) = inf;
        Segmentation.Clusters_sp(Segmentation.Clusters_sp == L2) = L1;
        Segmentation.Clusters_sp(Segmentation.Clusters_sp == inf) = L2;
        
    case 'NN'
        Segmentation.Clusters_nn(Segmentation.Clusters_nn == L1) = inf;
        Segmentation.Clusters_nn(Segmentation.Clusters_nn == L2) = L1;
        Segmentation.Clusters_nn(Segmentation.Clusters_nn == inf) = L2;
        
    case 'SOM'
        Segmentation.Clusters_som(Segmentation.Clusters_som == L1) = inf;
        Segmentation.Clusters_som(Segmentation.Clusters_som == L2) = L1;
        Segmentation.Clusters_som(Segmentation.Clusters_som == inf) = L2;

end

end

