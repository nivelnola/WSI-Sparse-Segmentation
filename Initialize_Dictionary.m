function D_0 = Initialize_Dictionary(DataStream, numAtoms, colormodel, lambda)
% Given input data, returns the learned dictionary

fprintf('Initializing dictionary...\n')
tic

% Step 1 - Unpartitioned dictionary containing random pre-existing patches
fprintf('\tSelecting random patches... ')
D_0 = DataStream(:, randperm(size(DataStream, 2), numAtoms));
if colormodel == "Quaternion"
    D_0 = concat2quaternion(D_0);
end
fprintf('Complete!\n')

% Step 2 - Train to reconstruct dataset
fprintf('\tLearning dictionary...\n')
switch colormodel
    case "Concatenation"
        opts = set_opts("init_learn_concat");
        D_0 = bpdndl(D_0, DataStream, lambda, opts);
    case "Quaternion"
        DataStream = concat2quaternion(DataStream);
        opts = set_opts("init_learn_quater");
        opts.K = numAtoms;
        opts.initialDictionary = D_0;
        [D_0, ~] = K_QSVD(DataStream, opts);
        D_0 = quaternion2concat(D_0);
end

t = toc;
fprintf('\tComplete!\n')
fprintf('Runtime: %f sec\n', t)

end

