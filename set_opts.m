function opts = set_opts(str)
%set_opts Simply imports options for dictionary learning and coding

switch str
    % Initialization - Learning
    case "init_learn_concat"
        opts.Verbose            = 1;
        opts.MaxMainIter        = 300;
        opts.AutoSigma          = 1;
        opts.AutoRho            = 1;
        opts.AutoSigmaScaling   = 1;
        opts.AutoRhoScaling     = 1;
    case "init_learn_quater"
        % Make sure to specify opts.K and opts.InitialDictionary!
        opts.displayProgress      = 1;
        opts.numIteration         = 5;
        opts.preserveDCAtom       = 0;
        opts.InitializationMethod = 'GivenMatrix';
        opts.errorFlag            = 0;
        opts.L                    = 5;
        
    % Initialization - Coding
    case "init_code_concat"
        opts.Verbose            = 1;
        opts.MaxMainIter        = 150;
    case "init_code_quater"
        opts.L                  = 5; 
        
    % Optimization - Learning
    case "optim_learn"
        opts.Verbose            = 0;
        opts.MaxMainIter        = 75;
        opts.AutoSigma          = 1;
        opts.AutoRho            = 1;
        opts.AutoSigmaScaling   = 1;
        opts.AutoRhoScaling     = 1; 

    % Optimization - Coding
    case "optim_code_concat"
        opts.Verbose            = 0;
        opts.MaxMainIter        = 25;
    case "optim_code_quater"
        opts.L                  = 5; 
        
    % Testing - coding
    case "test_code"
        opts.Verbose            = 0;
        opts.MaxMainIter        = 100;
end

end

