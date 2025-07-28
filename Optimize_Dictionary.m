function [Dictionary] = Optimize_Dictionary(Dictionary, DataStream, lambda, colormodel, numOptimIter)
% Given dictionary clusters, optimize them.
% Needs to be updated for quaternions!

numDicts = size(Dictionary, 2);
numDataPts = size(DataStream, 2);

fprintf('-----------------------------------------------\n')
for iteration = 1:numOptimIter
    % Report iteration number
    fprintf('Iteration %i/%i\n', iteration, numOptimIter)
    
    % Step 1 - Assignment
    fprintf('Step 1 - Assignment\n')
    
    fprintf('\tReconstructing patches... \n')
    A = cell(size(Dictionary));
    for dict_ticker = 1:numDicts
        fprintf('\t\tDictionary %i... ', dict_ticker)
        switch colormodel
            case "Concatenation"
                A{dict_ticker} = sparse(bpdn(Dictionary{dict_ticker}, DataStream, lambda, set_opts("optim_code_concat")));
            case "Quaternion"
                A{dict_ticker} = QOMP(concat2quaternion(Dictionary{dict_ticker}), concat2quaternion(DataStream), set_opts("optim_code_quater").L);
                A{dict_ticker} = sparse(quaternion2concat(A{dict_ticker}));
        end
        fprintf('Complete!\n')
    end
    
    fprintf('\tCalculating classification metric... \n')
    R_hat = nan(numDicts, numDataPts, 1);
    for dict_ticker = 1:numDicts
        fprintf('\t\tDictionary %i... ', dict_ticker)
        switch colormodel
            case "Concatenation"
                R_hat(dict_ticker, :) = calculate_reconstruction_metric(DataStream, Dictionary{dict_ticker}, A{dict_ticker}, lambda);
            case "Quaternion"
                R_hat(dict_ticker, :) = sqrt(sum((concat2quaternion(DataStream) - Qmult(concat2quaternion(D_0), concat2quaternion(A{dict_ticker}))).^2)) / size(DataStream, 2);
        end
        fprintf('Complete!\n')
    end
    
    fprintf('\tSeparating patches... \n')
    [~, classes] = max(R_hat);
    X = cell(size(Dictionary));
    for dict_ticker = 1:numDicts
        fprintf('\t\tDictionary %i... ', dict_ticker)
        X{dict_ticker} = DataStream(:, classes == dict_ticker);
        fprintf('Complete!\n')
    end
    
    tabulate(classes);

    % Step 2 - Update
    fprintf('\n\nStep 2 - Update\n')
    
    fprintf('\tUpdating dictionaries...\n')
    for dict_ticker = 1:numDicts
        fprintf('\t\tDictionary %i... ', dict_ticker)
        Dictionary{dict_ticker} = bpdndl(Dictionary{dict_ticker}, X{dict_ticker}, lambda, set_opts("optim_learn"));
        fprintf('Complete!\n')
    end
    fprintf('-----------------------------------------------\n')
end