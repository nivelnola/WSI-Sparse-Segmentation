%% Extract All Glomeruli
%   Given a folder of WSIs and a corresponding folder of XMLs, extract
%   all the glomeruli and save them
%
%   Alon S. Levin

%% Prepare the environment
clear, clc, close all
format compact
profile on

%% Prepare openslide library
openslide_load_library
disp(['OpenSlide version: ',openslide_get_version()])
fprintf('===============================================\n')

%% Obtain relevant folders
%{
The library can read virtual slides in the following formats:
    Aperio (.svs, .tif)
    Hamamatsu (.vms, .vmu, .ndpi)
    Leica (.scn)
    MIRAX (.mrxs)
    Philips (.tiff)
    Sakura (.svslide)
    Trestle (.tif)
    Ventana (.bif, .tif)
    Generic tiled TIFF (.tif)
%}

% Obtain WSI file list
% IMG_folder = uigetdir('F:\PAS_Dataset', 'Select IMG Directory');
IMG_folder = 'F:\PAS_Dataset\IMG';
IMG_fileList = dir(IMG_folder);
IMG_fileList = IMG_fileList(3:end);

% Obtain XML file list
% SEG_folder = uigetdir('F:\PAS_Dataset', 'Select SEG Directory');
SEG_folder = 'F:\PAS_Dataset\SEG';
SEG_fileList = dir(SEG_folder);
SEG_fileList = SEG_fileList(3:end);

fprintf('IMG Folder\n\t%s\n\t%i files found\n', IMG_folder, length(IMG_fileList))
fprintf('SEG Folder\n\t%s\n\t%i files found\n', SEG_folder, length(SEG_fileList))

% Check that files match
IMG_names = string(regexp({IMG_fileList.name}', '^(.+?)\.', 'tokens'));
SEG_names = string(regexp({SEG_fileList.name}', '^(.+?)\.', 'tokens'));

assert(length(IMG_names) == length(SEG_names), 'Error: Number of files do not match.')
fprintf('Passed check: File count.\n')
numFiles = length(IMG_names);

assert(sum(IMG_names == SEG_names) == length(IMG_names), 'Error: File names do not match.')
fprintf('Passed check: File names.\n')
fprintf('===============================================\n')

% Allow to review the files here
flag_filecheck = input('Check file names? Y/N [N]: ', 's');
if flag_filecheck == 'Y' || flag_filecheck == 'y'
    
    disp(table(string({IMG_fileList.name}'), string({SEG_fileList.name}'), IMG_names == SEG_names, ...
        'VariableNames', {'IMG', 'SEG', 'Match'}))
    
    input('Press any key to continue...')
end
fprintf('===============================================\n')

%% Iterate through the files

Glom_files = cell(1, 41);

for file_ticker = 1:numFiles
    
    % Get curent file objects
    curr_IMG_file = IMG_fileList(file_ticker);
    curr_SEG_file = SEG_fileList(file_ticker);
    
    fprintf('File %i/%i\n', file_ticker, numFiles)
    fprintf('\tName:\t%s\n', IMG_names(file_ticker))
    
    % Convert filepaths to usable objects
    curr_IMG = openslide_open([curr_IMG_file.folder, '\', curr_IMG_file.name]);
    curr_XML = xml2struct([curr_SEG_file.folder, '\', curr_SEG_file.name]);
    [mppX, mppY, slide_width, slide_height, ~, ~, ~] ...
        = openslide_get_slide_properties(curr_IMG);
    
    fprintf('File opened successfully.\n')
    fprintf('\tmppX:\t%f\n\tmppY:\t%f\n\twidth:\t%i\n\theight:\t%i\n', ...
        mppX, mppY, slide_width, slide_height);
    
    % Extract Gloms from current XML
    Glom_files{file_ticker}.Gloms = extract_all_gloms_SINGLE(curr_XML, curr_IMG, true);
    
    % Save gloms and data for current file
    Glom_files{file_ticker}.mppX = mppX;
    Glom_files{file_ticker}.mppY = mppY;
    Glom_files{file_ticker}.FileName = IMG_names(file_ticker);
    
    % Close the file
    openslide_close(curr_IMG)
    
    fprintf('File complete.\n')
    fprintf('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n')
end
fprintf('===============================================\n')

%% Save output
save('F:\PAS_Dataset\AllGloms.mat', 'Glom_files', '-v7.3')
fprintf('Output saved to:\n\t%s\n', 'F:\PAS_Dataset\AllGloms.mat')
clear