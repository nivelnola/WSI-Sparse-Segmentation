function [DataStream, Gloms, means] = Obtain_Data(numGloms, sizePatches, downsample, colorspace, numFeatures, removeMeans, include_segments)

%% Obtain gloms
GoodGloms_file = 'F:\PAS_Dataset\GoodGloms.mat';
fprintf('Importing database:\n\t%s\n', GoodGloms_file)
if ~exist('Glom_files_clean', 'var') || ~exist('Gloms_Table', 'var')
    load(GoodGloms_file)
    fprintf('Complete!\n')
else
    fprintf('File already loaded!\n')
end
fprintf('===============================================\n')

fprintf('Choosing training gloms...\n')
ID_Table = Gloms_Table(randperm(size(Gloms_Table,1), numGloms), 1:2);
Gloms = cell(numGloms, 1);
for glom_ticker = 1:numGloms
    FileID = ID_Table.FileID(glom_ticker);
    GlomID = ID_Table.GlomID(glom_ticker);
    fprintf('\t%i.%i\t==>\tGlom %i...\t', FileID, GlomID, glom_ticker)
    Gloms{glom_ticker} = Glom_files_clean{FileID}.Gloms{GlomID};
    fprintf('Complete!\n')
end

fprintf('Unloading database... ')
clear Gloms_files_clean Gloms_Table
fprintf('Complete!\n')
fprintf('===============================================\n')

%% Obtain DataStream
[DataStream, means] = obtain_datastream(Gloms, sizePatches, downsample, colorspace, numFeatures, removeMeans, include_segments);
    
end