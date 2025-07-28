%% Summarize glom data
%   Load the saved glomeruli items, and summarize them in a csv
%
%   Alon S. Levin

%% Prepare the environment
%clear, clc, close all
format compact

%% Load AllGloms file
AllGloms_file = 'F:\PAS_Dataset\AllGloms.mat';
fprintf('Loading file:\n\t%s\n', AllGloms_file)
if ~exist('Glom_files', 'var')
    load(AllGloms_file)
    fprintf('Complete!\n')
else
    fprintf('File already loaded!\n')
end
fprintf('===============================================\n')

%% Prepare CSV file

% Open the file for writing
CSV_file = 'F:\PAS_Dataset\AllGloms.csv';
fprintf('Writing to:\n\t%s\n', CSV_file)
fid = fopen(CSV_file, 'w');

% Write a header
cellHeader = {'FileID', 'GlomID', 'NumVertices', 'mppX', 'mppY', 'AreaPixl', 'AreaTrue', 'Circularity'};
textHeader = strjoin(cellHeader, ',');
fprintf(fid,'%s\n', textHeader);
fprintf('Headers:\n\t%s\n', textHeader);
fprintf('===============================================\n')

%% Iterate through Glom_files
for file_ticker = 1:length(Glom_files)
    curr_Glom_file = Glom_files{file_ticker};
    
    fprintf('Reading file %i/%i... ', file_ticker, length(Glom_files));
    
    % Iterate through gloms in curr_Glom_file
    for Glom_ticker = 1:length(curr_Glom_file.Gloms)
        curr_Glom = curr_Glom_file.Gloms{Glom_ticker};
        
        % FileID, GlomID, NumVertices, mppX, mppY
        curr_Line = [num2str(file_ticker), ',', ...
                     num2str(Glom_ticker), ',', ...
                     num2str(size(curr_Glom.Vertices, 1)), ',' ...
                     num2str(curr_Glom_file.mppX), ',', ...
                     num2str(curr_Glom_file.mppX)];
                 
        % Calculate region props of interest
        area_pixl = regionprops(curr_Glom.Mask, 'area');
        area_pixl = sum([area_pixl(:).Area]);
        area_true = area_pixl * curr_Glom_file.mppX * curr_Glom_file.mppY;
        circularity = regionprops(curr_Glom.Mask, 'Circularity').Circularity;
                
        % Append region props of interest
        curr_Line = [curr_Line, ',', ...
                     num2str(area_pixl), ',', ...
                     num2str(area_true), ',', ...
                     num2str(circularity)];
        
        % Write to file
        fprintf(fid,'%s\n', curr_Line); 
    end
    fprintf('Complete!\n');
end
fprintf('===============================================\n')

%% Close the file
fclose(fid);
fprintf('Completed summary. It can be found at:\n\t%s\n', CSV_file);