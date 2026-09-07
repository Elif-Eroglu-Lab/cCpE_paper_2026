function step1() 

% Compile one master excel file
%
% The code prompts selection of a folder that contains the subfolders.
% NOTE: Each subfolder corresponds to a section. 
% 
% Multiple ROIs make up one continuous section.
% We had a maximum of 6 ROI in one section. The ImageJ multi-measure output csv was saved in the format 
% R1_Results.csv ... R6_Results.csv corresponding to ROI_1 ... ROI_6
%  
% NOTE: Column 10 corresponds to (Header: Slice, Row value 1, 2, 3) 
% NOTE: Column 5 corresponds to (Header: Y, ImageJ nuclei centroid Y) 
%
% Inside each group (1, 2, 3) Column 5 is sorted from smallest to largest.
% A "Source_File" column is added to check which CSV file the row comes from. 

mainFolder = uigetdir('', 'Select the folder');
if mainFolder == 0
    disp('Operation canceled.');
    return
end
fprintf('Selected folder: %s\n', mainFolder);

% --- find the subfolders ---
items = dir(mainFolder);
items = items([items.isdir]);                                  % keep folders
items = items(~ismember({items.name}, {'.', '..'}));            % drop . and ..
fprintf('Found %d subfolders.\n\n', numel(items));

fileNames = ["R1_Results.csv" "R2_Results.csv" "R3_Results.csv" ...
             "R4_Results.csv" "R5_Results.csv" "R6_Results.csv"];

timestamp  = datestr(now, 'yyyy-mm-dd_HHMMSS');
outputFile = fullfile(mainFolder, ['step1_Results_' timestamp '.xlsx']);

foldersWritten = 0;

% --- go through each subfolder ---
for k = 1:numel(items)
    folderName = items(k).name;
    fprintf('Processing folder %d/%d: %s\n', k, numel(items), folderName);

    allData    = table();
    filesFound = 0;

    for f = 1:numel(fileNames)
        filePath = fullfile(mainFolder, folderName, fileNames(f));
        if ~isfile(filePath)
            continue
        end

        try
            data = readtable(filePath, 'PreserveVariableNames', true);
        catch ME
            warning('Could not read %s: %s', fileNames(f), ME.message);
            continue
        end

        % remember which file the rows came from
        [~, label] = fileparts(fileNames(f));
        data.Source_File = repmat({char(label)}, height(data), 1);

        allData    = [allData; sortRows(data)]; %#ok<AGROW>
        filesFound = filesFound + 1;
    end

    if isempty(allData)
        continue
    end

    % sheet names: max 31 characters, no \ / * ? : [ ]
    sheetName = regexprep(folderName, '[\\/*?:\[\]]', '_');
    sheetName = sheetName(1:min(31, end));

    writetable(allData, outputFile, 'Sheet', sheetName, 'WriteMode', 'overwritesheet');
    fprintf('  Created sheet: %s (%d rows from %d files)\n', ...
        sheetName, height(allData), filesFound);
    foldersWritten = foldersWritten + 1;
end

fprintf('\nProcessing complete. %d/%d folders processed.\n', foldersWritten, numel(items));
fprintf('Output file saved to:\n%s\n', outputFile);

if ispc
    winopen(outputFile);
end
end


function out = sortRows(data)
% Keep only rows where column 10 is 1, 2 or 3, then sort by column 10,
% and by column 5 within each group. If nothing is left, return the
% data unchanged.

keep = ismember(data{:, 10}, [1 2 3]);
if any(keep)
    out = sortrows(data(keep, :), [10 5]);
else
    out = data;
end
end