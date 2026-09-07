function step3()
% Clone size analysis.
%
% When the value is below THRESHOLD the two cells are "linked", so a clone
% is a run of consecutive rows whose links are all below threshold.
%
% NOTE: Column 13 of a "step3_Results.xlsx" sheet is "coefficient" 
% The file being read is already arranged based on Y. Thus the rows are ordered and adjacent rows are the cells to compare. 
% 
% Output: one column per sheet/section, counting total cells, singlets, and the count of clones of each size.

    THRESHOLD = 0.05;   % below this, two adjacent cells are linked
    AE_COL    = 13;     % column index of the coefficient measurement

    fileNames  = {'step3_Results.xlsx'};
    outputName = 'Clone_ccpe7.xlsx';

    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir)
        scriptDir = pwd;
    end

    % Remove any existing output so a re-run cannot append to stale results.
    outputFile = fullfile(scriptDir, outputName);
    if isfile(outputFile)
        delete(outputFile);
    end

    % Go through each input file 
    for f = 1:numel(fileNames)
        inputFile = fullfile(scriptDir, fileNames{f});
        fprintf('\nProcessing file %d/%d: %s\n', f, numel(fileNames), fileNames{f});

        if ~isfile(inputFile)
            warning('File not found: %s - skipping.', fileNames{f});
            continue
        end

        try
            sheets = cellstr(sheetnames(inputFile));
        catch ME
            warning('Cannot open %s: %s - skipping.', fileNames{f}, ME.message);
            continue
        end

        nSheets = numel(sheets);
        fprintf('Found %d sheets in the file\n', nSheets);

        % Analyse each sheet 
        S       = struct('total', num2cell(zeros(nSheets, 1)), 'sizes', []);
        maxSize = 1;

        for i = 1:nSheets
            fprintf('Processing sheet %d/%d: %s\n', i, nSheets, sheets{i});

            try
                data = readtable(inputFile, 'Sheet', sheets{i}, 'VariableNamingRule', 'preserve');
            catch ME
                warning('Could not read sheet %s: %s', sheets{i}, ME.message);
                continue
            end

            if width(data) < AE_COL
                warning('Sheet %s has fewer than %d columns - skipped.', sheets{i}, AE_COL);
                continue
            end

            % Keep the column up to the last non-blank row
            ae   = asNum(data{:, AE_COL}, sheets{i});
            last = find(~isnan(ae), 1, 'last');
            if isempty(last)
                warning('No valid data in sheet %s.', sheets{i});
                continue
            end
            ae = ae(1:last);
            n  = numel(ae);

            S(i).total = n;
            S(i).sizes = cloneSizes(ae, THRESHOLD);
            maxSize    = max(maxSize, max(S(i).sizes));

            % Console summary for this sheet
            fprintf('  Total cells: %d\n', n);
            for sz = 1:max(S(i).sizes)
                k = sum(S(i).sizes == sz);
                if k > 0
                    fprintf('  %-15s: %d\n', sizeLabel(sz), k);
                end
            end

            % Every cell must be accounted for exactly once
            if sum(S(i).sizes) ~= n
                warning('Cell count mismatch in sheet %s: accounted=%d, total=%d', ...
                    sheets{i}, sum(S(i).sizes), n);
            end
        end

        % Building the output table for this file 
        % rows: total, singlets, then one row per clone size seen in any sheet
        labels = [{'Total Cells'; 'Singlets'}; ...
                  arrayfun(@sizeLabel, (2:maxSize)', 'UniformOutput', false)];

        out = cell(numel(labels), nSheets);
        for i = 1:nSheets
            if S(i).total == 0
                out(:, i) = {'No data'};
                continue
            end
            out{1, i} = S(i).total;
            for sz = 1:maxSize
                out{sz + 1, i} = sum(S(i).sizes == sz);
            end
        end

        table = [[{'Metric'}, sheets(:)']; [labels, out]];

        [~, base] = fileparts(fileNames{f});
        name = sheetName(base);
        writecell(table, outputFile, 'Sheet', name);
        fprintf('Written to output sheet: %s\n', name);
    end

    fprintf('\nProcessing complete!\n');
    fprintf('Output file: %s\n', outputFile);

    if ispc
        try, winopen(outputFile); catch, end
    end
end


% helpers

function sizes = cloneSizes(ae, threshold)
% Sizes of every clone in a column, singlets included.
% linked(i) is true when cell i is joined to cell i+1. The last value is
% never tested because there is no row below it.
    n = numel(ae);
    if n < 2
        sizes = 1;
        return
    end

    linked = ae(1:end-1) < threshold;

    % Pad with a 0 at each end and difference: +1 marks the first link of a
    % run, -1 the position just after the last one. Padding catches runs
    % that touch either end.
    d      = diff([0; linked(:); 0]);
    starts = find(d ==  1);      % first cell of each clone
    ends   = find(d == -1);      % last cell of each clone

    sizes = (ends - starts + 1)';   % L links join L+1 cells

    % anything not inside a run is a clone of one
    covered = false(n, 1);
    for r = 1:numel(starts)
        covered(starts(r):ends(r)) = true;
    end
    sizes = [sizes, ones(1, sum(~covered))];
end


function v = asNum(col, sName)
% Turn a table column into a numeric column vector.
% A well-formed step3 sheet reads back as numeric, so the cell branch
% should never run. It warns if it does, because that means something
% unexpected is sitting in the column.
    if isnumeric(col)
        v = double(col(:));
    elseif iscell(col)
        warning('Coefficient column in sheet %s is not numeric - converting from text.', sName);
        v = cellfun(@str2double, cellfun(@num2str, col, 'UniformOutput', false));
    else
        warning('Coefficient column in sheet %s is not numeric - converting from text.', sName);
        v = double(string(col));
    end
    v = v(:);
end


function label = sizeLabel(sz)
% Row label for a clone of a given size.
    if sz == 1
        label = 'Singlets';
    else
        label = sprintf('%d-cell clones', sz);
    end
end


function name = sheetName(raw)
% Valid Excel sheet name: max 31 characters, no : \ / ? * [ ]
    raw(ismember(raw, ':\/?*[]')) = '_';
    name = raw(1:min(31, end));
end