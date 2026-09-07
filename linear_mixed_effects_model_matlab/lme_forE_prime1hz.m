close all; clc; clear all; 
filename = 'input.xlsx';
 
Eprime = [];
Condition = {};
Animal = {};
 
% naive
for s = [1 2 3 4]
    col2 = readmatrix(filename, 'Sheet', s, 'Range', 'B:B');
    col2 = col2(~isnan(col2));
    Eprime = [Eprime; col2];
    Condition = [Condition; repmat({'naive'}, numel(col2), 1)];
    Animal = [Animal; repmat({sprintf('Animal%d', s)}, numel(col2), 1)];
end
 
% 7dpci
for s = [5 6 7 8]
    col2 = readmatrix(filename, 'Sheet', s, 'Range', 'B:B');
    col2 = col2(~isnan(col2));
    Eprime = [Eprime; col2];
    Condition = [Condition; repmat({'7dpci'}, numel(col2), 1)];
    Animal = [Animal; repmat({sprintf('Animal%d', s)}, numel(col2), 1)];
end
 
% 10dpci
for s = [9 10 11 12]
    col2 = readmatrix(filename, 'Sheet', s, 'Range', 'B:B');
    col2 = col2(~isnan(col2));
    Eprime = [Eprime; col2];
    Condition = [Condition; repmat({'10dpci'}, numel(col2), 1)];
    Animal = [Animal; repmat({sprintf('Animal%d', s)}, numel(col2), 1)];
end
 
% 14dpci
for s = [13 14 15 16]
    col2 = readmatrix(filename, 'Sheet', s, 'Range', 'B:B');
    col2 = col2(~isnan(col2));
    Eprime = [Eprime; col2];
    Condition = [Condition; repmat({'14dpci'}, numel(col2), 1)];
    Animal = [Animal; repmat({sprintf('Animal%d', s)}, numel(col2), 1)];
end
 
% 18dpci
for s = [17 18 19]
    col2 = readmatrix(filename, 'Sheet', s, 'Range', 'B:B');
    col2 = col2(~isnan(col2));
    Eprime = [Eprime; col2];
    Condition = [Condition; repmat({'18dpci'}, numel(col2), 1)];
    Animal = [Animal; repmat({sprintf('Animal%d', s)}, numel(col2), 1)];
end
 
% Building the table fitlme needs
tbl = table(Eprime, Condition, Animal);
tbl.Condition = categorical(tbl.Condition, {'naive','7dpci','10dpci','14dpci','18dpci'});
tbl.Animal = categorical(tbl.Animal);
 
% linear mixed-effects model: Condition the fixed-effect, Animal as the random-effect
lme = fitlme(tbl, 'Eprime ~ Condition + (1|Animal)');
lme
 
anova(lme, 'DFMethod', 'satterthwaite') % to test whether Condition matters
 
% Pairwise comparisons between every pair of conditions using coefTest
condNames = categories(tbl.Condition);
coefNames = lme.CoefficientNames;
nComparisons = 4;  

fprintf('\n%-12s %-12s %10s %10s\n', 'Group A', 'Group B', 'pValue', 'pBonf');
for i = 2:numel(condNames)  % Start from 2 to skip naive
    H = zeros(1, numel(coefNames));
    idxA = find(strcmp(coefNames, ['Condition_' condNames{1}]));  % naive is the reference
    idxB = find(strcmp(coefNames, ['Condition_' condNames{i}]));
    if ~isempty(idxA), H(idxA) = 1; end
    if ~isempty(idxB), H(idxB) = -1; end

    pVal = coefTest(lme, H, 0, 'DFMethod', 'satterthwaite');

    pBonf = min(pVal * nComparisons, 1);
    
    fprintf('%-12s %-12s %10.5f %10.5f\n', condNames{1}, condNames{i}, pVal, pBonf);
end 