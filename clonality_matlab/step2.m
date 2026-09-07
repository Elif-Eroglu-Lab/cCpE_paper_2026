function step2()

% Background-subtract, calculate channel weight, and color spread coeffiennt
% 
% The code prompts selection of an excel file (select step1 output). Then, per sheet:
% Reads the rows belonging to R1_Results ... R6_Results (column 18), and splits them by slice number 1, 2, 3 (column 10). 
% Read column 3, and these become channel_1, 2, 3.
% 
% Sheets are grouped in blocks of 3 consecutive sheets (= one animal).
% The background of a channel is the smallest value of that channel across the 3 sheets of its group.
% 
% Subtract background for each channel: bg_sub_Mean_n = channel_n - background_n
% SUM = bg_sub_Mean_1 + bg_sub_Mean_2 + bg_sub_Mean_3
% Relative channel weight I_cer, I_GFP, I_mcherry is mean divided by SUM
% Threshold = "Yes" when all three differences are smaller than 0.05. 
% distance = sqrt of the squared differences,
% coefficient = distance / sqrt(2).

    [filename, pathname] = uigetfile('*.xlsx', 'Select Step1_output Excel file');
    if isequal(filename, 0)
        disp('Operation cancelled.');
        return
    end

    inputFile = fullfile(pathname, filename);
    timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
    outputFile = fullfile(pathname, ['step2_QC_' timestamp '.xlsx']);
    step3File = fullfile(pathname, ['step3_Results' timestamp '.xlsx']);

    sheets   = cellstr(sheetnames(inputFile));
    nSheets  = numel(sheets);
    fprintf('Found %d sheets in the file\n', nSheets);

    % Read the three channels out of every sheet 
    S = struct('name', sheets, 'ch', []);   
    maxRows = 0;

    for i = 1:nSheets
        fprintf('Processing sheet %d/%d: %s\n', i, nSheets, sheets{i});
        S(i).ch = {[] [] []};

        try
            data = readtable(inputFile, 'Sheet', sheets{i}, 'PreserveVariableNames', true);
            if width(data) < 18
                warning('Sheet %s has fewer than 18 columns - skipped.', sheets{i});
            else
                source = asText(data{:, 18});   % R1_Results ... R5_Results
                slice  = asNum(data{:, 10});    % 1, 2 or 3
                value  = asNum(data{:,  3});    % the measurement

                for r = 1:5                                  % R1 ... R5, in order
                    isSource = contains(source, sprintf('R%d_Results', r), 'IgnoreCase', true);
                    for c = 1:3                              % slice 1, 2, 3
                        v = value(isSource & slice == c);
                        S(i).ch{c} = [S(i).ch{c}; v(~isnan(v))];
                    end
                end
            end
        catch ME
            warning('Error processing sheet %s: %s', sheets{i}, ME.message);
        end

        maxRows = max(maxRows, max(cellfun(@numel, S(i).ch)));
    end

    % Calculate background per animal by looking through sections (3 sheets) 
    for first = 1:3:nSheets
        members = first:min(first + 2, nSheets);

        bg = nan(1, 3);
        for c = 1:3
            pooled = vertcat(S(members).ch);          % all sheets of the group
            pooled = vertcat(pooled{:, c});
            if ~isempty(pooled), bg(c) = min(pooled); end
        end

        for i = members
            S(i).bg    = bg;
            S(i).isTop = (i == first);                % backgrounds printed once
        end
    end

    % Calculations per sheet 
    for i = 1:nSheets
        n = max(cellfun(@numel, S(i).ch));

        ch   = nan(n, 3);
        mean = nan(n, 3);
        for c = 1:3
            ch(:, c)   = padNaN(S(i).ch{c}, n);
            mean(:, c) = ch(:, c) - S(i).bg(c);
        end

        total = sum(mean, 2);
        ratio = mean ./ total;
        ratio(~isfinite(total) | total == 0, :) = NaN;    % guard divide-by-zero

        % compare each row with the next one; after the last row use 0
        d = ratio - [ratio(2:end, :); zeros(1, 3)];

        S(i).n     = n;
        S(i).ch    = ch;
        S(i).mean  = mean;
        S(i).total = total;
        S(i).ratio = ratio;
        S(i).dist  = sqrt(sum(d.^2, 2));
        S(i).coef  = S(i).dist / sqrt(2);
        S(i).yes   = all(abs(d) < 0.05, 2);
    end

    % building the step2 sheet 
    W   = 14;                                   % columns per section
    out = cell(maxRows + 2, nSheets * W);       % row 1 = name/bg, row 2 = headers

    headers = {'channel_1', 'channel_2', 'channel_3', ...
               'bg_sub_Mean_1', 'bg_sub_Mean_2', 'bg_sub_Mean_3', 'SUM', ...
               'I_cer', 'I_GFP', 'I_mcherry', 'Threshold', 'distance', 'coefficient'};

    for i = 1:nSheets
        c0 = (i - 1) * W + 1;                   % first column of this section

        out{1, c0} = S(i).name;
        if S(i).isTop
            for c = 1:3
                out = putNum(out, 1, c0 + 3 + c, S(i).bg(c));
            end
        end
        out(2, c0 + (1:13)) = headers;

        rows = [S(i).ch, S(i).mean, S(i).total, S(i).ratio, S(i).dist, S(i).coef];
        for r = 1:S(i).n
            for c = [1:10 12 13]                % 11 is the text column
                src = c - (c > 11);             % skip Threshold in "rows"
                out = putNum(out, r + 2, c0 + c, rows(r, src));
            end
            if S(i).yes(r)
                out{r + 2, c0 + 11} = 'Yes';
            end
        end
    end

    writeSheet(out, outputFile, 'Extracted_Data');

    % ---------- 5. split the same layout into one sheet per section ----------
    used = {};
    for i = 1:nSheets
        c0      = (i - 1) * W + 1;
        block   = out(2:S(i).n + 2, c0 + (1:13));       % headers + data
        [nm, used] = sheetName(out{1, c0}, used, i);
        writeSheet(block, step3File, nm);
    end

    fprintf('\nProcessing complete!\n');
    fprintf('Input file:        %s\n', inputFile);
    fprintf('Step2 output file: %s\n', outputFile);
    fprintf('Step3 output file: %s\n', step3File);
    fprintf('Total sheets processed: %d\n', nSheets);

    if ispc
        try, winopen(outputFile); winopen(step3File); catch, end
    end
end


function t = asText(col)
% Turn a table column into text, whatever type it came in as.
    if isnumeric(col)
        t = cellstr(num2str(col(:)));
    elseif iscell(col)
        t = cellfun(@num2str, col, 'UniformOutput', false);
    else
        t = string(col);
    end
end


function v = asNum(col)
% Turn a table column into a numeric column vector.
    if isnumeric(col)
        v = double(col(:));
    elseif iscell(col)
        v = cellfun(@str2double, cellfun(@num2str, col, 'UniformOutput', false));
    else
        v = double(string(col));
    end
    v = v(:);
end


function v = padNaN(v, n)
% Column vector of length n, padded with NaN.
    v = v(:);
    v(end + 1:n, 1) = NaN;
end


function out = putNum(out, r, c, val)
% Write a number; leave the cell blank for NaN or Inf.
    if isfinite(val)
        out{r, c} = val;
    end
end


function writeSheet(c, file, name)
    writecell(c, file, 'Sheet', name);
end


function [name, used] = sheetName(raw, used, idx)
% Valid, unique Excel sheet name: max 31 characters, no : \ / ? * [ ]
    if ~ischar(raw), raw = char(string(raw)); end
    raw(ismember(raw, ':\/?*[]')) = '_';
    raw = strtrim(raw);
    if isempty(raw)
        raw = sprintf('Section_%d', idx);
    end
    raw = raw(1:min(31, end));

    name = raw;
    k = 2;
    while any(strcmpi(name, used))
        suffix = sprintf('_%d', k);
        name   = [raw(1:min(numel(raw), 31 - numel(suffix))) suffix];
        k = k + 1;
    end
    used{end + 1} = name;
end