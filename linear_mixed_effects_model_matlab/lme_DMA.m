%  https://www.mathworks.com/help/stats/referencelist.html?type=function&s_tid=CRUX_topnav&category=mixed-effects-1
%
%  3 conditions: naive, veh, cCpE
%  3 outcomes:   E', E'', tan delta
%  5 frequencies: 1, 2, 4, 10, 20 Hz
%
%  sheetIdx - animal; m - E',E'',tanD, r - row/measurement within animal
%  f - frequency index
%  locID is S<sheet>_R<row> so location labels are unique per animal - Location is nested in Animal
%
%  Model:  fitlme(T, 'Y ~ Condition * FreqLog + (1+FreqLog|Animal) + (1|Location)')
%  FreqLog = log10(frequency) since frequencies are not linearly spaced.
%  log10(1 Hz) = 0, so the intercept and Condition terms are evaluated at 1 Hz
%  Fixed effects: Condition, FreqLog, Condition:FreqLog
%  Random effects: animals (intercept + slope), location (intercept)
%  naive is the reference level - fitlme default dummy coding is 'reference'
%
%  Taking varying intercept and varying slope
%  introduce Condition:FreqLog interaction (without this term the model forces every condition's frequency-sweep trend to have the same slope differing only in intercept
%  because the frequency-sweep trends (slopes) could differ across timepoints (Condition)
%
%  Using ANOVA DFMethod='Satterthwaite' instead of the default. default value is 'residual', df = n - p
%  MATLAB calls these marginal F-tests, one per fixed-effects term (printed header says "ANOVA marginal tests")
%  The Condition row is not an averaged main effect - it tests condition differences at FreqLog = 0, i.e. 1 Hz
%
%  Variance components / ICC — pulls the random-effect covariance matrices out of lme.covarianceParameters
%  (a cell array: one entry per random-effect grouping term, in the order entered in the formula) plus the residual variance
%  Only cp{1}(1,1) is used = the animal INTERCEPT variance. Animal also carries a slope,
%  so the animal share of variance changes with frequency - these are ICCs at 1 Hz, not one ICC for the sweep
%  Check the 3rd output of covarianceParameters for boundary estimates (variance ~0, |corr| ~1) - 12 animals is few for a 2x2
%
%  If p value significant, linear hypothesis tests on the fixed-effects coefficients of the fitted linear mixed-effects model
%  coefTest https://se.mathworks.com/help/stats/linearmixedmodel.coeftest.html
%  - computes an F-test for a linear combination of coefficients using Satterthwaite-approximated denominator degrees of freedom 
%  - coefTest DFMethod ALSO defaults to 'residual', so it has to be passed in explicitly to match the anova
%
%  internal test statistic is F, with numerator df = rank of contrast matrix = 1, mathematically identical to a two-tailed t-test
%  The Bonferroni alpha is the adjusted significance threshold used to prevent false positives in multiple hypothesis testing.
%  divide standard alpha (e.g., 0.05) by the total number of comparisons being performed simultaneously
%
%  Diagnostic plot: fitted(lme) is conditional (includes random effects), residuals are Pearson
%  QC Plot: mean +- SEM across animals.

clc; clear; close all;
filePath = fullfile(fileparts(mfilename('fullpath')), 'input_veh_ccpe.xlsx');    

condOrder     = {'naive','veh','cCpE'};   
nCond         = numel(condOrder);
sheetsPerCond = [4 4 4];                     
nSheets       = sum(sheetsPerCond);          
sheetCond     = repelem(1:nCond, sheetsPerCond);   

frequencies = [1 2 4 10 20];
freqLog     = log10(frequencies);
nFreq       = numel(frequencies);

cols_Ep   = 1:5;     
cols_Epp  = 6:10;    
cols_tanD = 11:15;   

LOG_Y = [false false false];

measNames = {'Eprime','Edprime','tanDelta'};
colSets   = {cols_Ep, cols_Epp, cols_tanD};

%──────────────────────────────────────────────────────────────────────────
Y_all=cell(1,3); A_all=cell(1,3); C_all=cell(1,3); F_all=cell(1,3); L_all=cell(1,3);
for m = 1:3
    Y_all{m}=[]; A_all{m}=[]; C_all{m}={}; F_all{m}=[]; L_all{m}={};
end

for sheetIdx = 1:nSheets
    cIdx = sheetCond(sheetIdx);

    raw   = readmatrix(filePath, 'Sheet', sheetIdx, 'Range', 'B:P');
    nRows = size(raw, 1);
    fprintf('Sheet %2d | %-8s (animal #%d) : %d rows\n', ...
            sheetIdx, condOrder{cIdx}, sheetIdx, nRows);

    for m = 1:3
        block = raw(:, colSets{m});       
        for r = 1:nRows
            locID = sprintf('S%02d_R%03d', sheetIdx, r);
            for f = 1:nFreq
                v = block(r, f);
                if isnan(v); continue; end
                if LOG_Y(m)
                    if v <= 0; continue; end
                    v = log10(v);
                end
                Y_all{m}(end+1,1) = v;
                A_all{m}(end+1,1) = sheetIdx;          
                C_all{m}(end+1,1) = condOrder(cIdx);
                F_all{m}(end+1,1) = freqLog(f);
                L_all{m}(end+1,1) = {locID};
            end
        end
    end
end

fprintf('\nImport complete.  Animals per condition:\n')
for c = 1:nCond
    fprintf('  %-8s : %d\n', condOrder{c}, sheetsPerCond(c))
end
fprintf('\n')

%──────────────────────────────────────────────────────────────────────────
%  FIT LME MODELS  +  KEY TEST (all Condition:FreqLog interaction)
%
%  Y ~ Condition * FreqLog + (1+FreqLog|Animal) + (1|Location)
%──────────────────────────────────────────────────────────────────────────
lmeModels = struct();
keyP      = nan(1,3);

for m = 1:3
    mn = measNames{m};

    T = table( ...
        Y_all{m}, ...
        categorical(A_all{m}), ...
        categorical(C_all{m}, condOrder), ...
        F_all{m}, ...
        categorical(L_all{m}), ...
        'VariableNames', {'Y','Animal','Condition','FreqLog','Location'} );

    lme = fitlme(T, 'Y ~ Condition * FreqLog + (1+FreqLog|Animal) + (1|Location)');
    lmeModels.(mn) = lme;

    fprintf('==============================================\n')
    fprintf('  MODEL: %s   (n = %d obs)\n', mn, height(T))
    fprintf('==============================================\n\n')

    [~, ~, fe] = fixedEffects(lme);
    fprintf('Fixed Effects:\n'); disp(fe)

    av = anova(lme, 'DFMethod', 'Satterthwaite');   % Type III, omnibus per-term F-test
    fprintf('Type-III ANOVA (omnibus test per term):\n'); disp(av)

    try
        tn    = string(av.Term);
        ixInt = find(contains(tn,'Condition') & contains(tn,'FreqLog'), 1);
        keyP(m) = av.pValue(ixInt);
        verdict = 'no detected difference in trends';
        if keyP(m) < 0.05; verdict = '>>> TRENDS DIFFER across groups'; end
        fprintf('  KEY TEST  Condition:FreqLog   p = %.4g   %s\n\n', keyP(m), verdict)
    catch
        fprintf('  (could not auto-extract interaction row)\n\n')
    end

    cp     = lme.covarianceParameters;
    s2_an  = cp{1}(1,1);
    s2_loc = cp{2}(1,1);
    s2_res = lme.MSE;
    tot    = s2_an + s2_loc + s2_res;
    fprintf('  Variance:  animal=%.4g  location=%.4g  residual=%.4g\n', s2_an, s2_loc, s2_res)
    fprintf('  ICC(animal)=%.3f   ICC(animal+location)=%.3f\n\n', s2_an/tot, (s2_an+s2_loc)/tot)

    %──────────────────────────────────────────────────────────────
    %  DIAGNOSTIC: residuals vs fitted (all observations, this outcome)
    %──────────────────────────────────────────────────────────────
    fitVals = fitted(lme);                              % conditional fitted values (incl. random effects)
    resVals = residuals(lme, 'ResidualType', 'Pearson'); % standardized residuals

    figure('Color','w','Position',[50+300*(m-1) 620 560 420]); hold on
    scatter(fitVals, resVals, 25, 'filled', ...
        'MarkerFaceColor',[0.2 0.2 0.6], 'MarkerFaceAlpha', 0.5)
    yline(0, 'r--', 'LineWidth', 1.5)
    hold off
    set(gca,'FontSize',12,'Box','on','TickDir','in')
    xlabel('Fitted values')
    ylabel('Pearson residuals')
    title(sprintf('Residuals vs Fitted — %s', mn))

end

%──────────────────────────────────────────────────────────────────────────
%  POST-HOC — pairwise TREND (slope) differences via coefTest contrasts
%  (single model fit per outcome; no releveling/refitting needed)
%──────────────────────────────────────────────────────────────────────────
nPairs  = nCond*(nCond-1)/2;
alpha_b = 0.05 / nPairs;

fprintf('\n=======================================================================================\n')
fprintf('  POST-HOC pairwise TREND differences (coefTest) — Bonferroni alpha = %.4f\n', alpha_b)
fprintf('=============+++++++++++++++++++++++++++++++++++++++++=================================\n\n')

postHocTables = struct();

for m = 1:3
    mn  = measNames{m};
    lme = lmeModels.(mn);              % model already fit earlier, naive = reference

    [beta, ~, fe] = fixedEffects(lme);
    names = fe.Name;
    nCoef = numel(names);
    covB  = lme.CoefficientCovariance; % covariance matrix of fixed-effect estimates

    % Build a contrast row per condition: which coefficients sum to that
    % condition's slope (FreqLog base slope + its own interaction offset).
    slopeVec = zeros(nCond, nCoef);
    freqIdx  = find(strcmp(names, 'FreqLog'));
    slopeVec(:, freqIdx) = 1;                       % every condition includes base slope
    for c = 2:nCond
        idx = find(strcmp(names, ['Condition_' condOrder{c} ':FreqLog']));
        slopeVec(c, idx) = 1;                        % non-reference conditions add offset
    end

    comp = strings(nPairs,1);
    dSlope = nan(nPairs,1); SEv = nan(nPairs,1); tv = nan(nPairs,1);
    DFv = nan(nPairs,1); pv = nan(nPairs,1);

    row = 0;
    for i = 1:nCond-1
        for j = i+1:nCond
            row = row + 1;
            L   = slopeVec(j,:) - slopeVec(i,:);     % contrast: slope_j - slope_i

            est = L * beta;
            se  = sqrt(L * covB * L');
            [p, ~, ~, DF2] = coefTest(lme, L, 0, 'DFMethod', 'satterthwaite');

            comp(row)   = sprintf('%s vs %s', condOrder{i}, condOrder{j});
            dSlope(row) = est;
            SEv(row)    = se;
            tv(row)     = est / se;
            DFv(row)    = DF2;
            pv(row)     = p;
        end
    end

    sigLabel = strings(nPairs,1);
    for k = 1:nPairs
        if     pv(k) < 0.001;    sigLabel(k) = "***";
        elseif pv(k) < 0.01;     sigLabel(k) = "**";
        elseif pv(k) < alpha_b;  sigLabel(k) = "*  (Bonf)";
        elseif pv(k) < 0.05;     sigLabel(k) = "~  (uncorr)";
        else;                    sigLabel(k) = "";
        end
    end

    postHocT = table(comp, dSlope, SEv, tv, DFv, pv, sigLabel, ...
        'VariableNames', {'Comparison','dSlope','SE','tStat','DF','pValue','Sig'});

    fprintf('--- %s ---\n', mn)
    disp(postHocT)
    fprintf('\n')

    postHocTables.(mn) = postHocT;    % for later export/reporting
end


%──────────────────────────────────────────────────────────────────────────
%  SUMMARY
%──────────────────────────────────────────────────────────────────────────
fprintf('==============================================\n')
fprintf('  SUMMARY — "do frequency trends differ across groups?"\n')
fprintf('==============================================\n')
for m = 1:3
    flag = ''; if keyP(m) < 0.05; flag = '  <-- differ'; end
    fprintf('  %-9s  Condition:FreqLog  p = %.4g%s\n', measNames{m}, keyP(m), flag)
end
fprintf('\n')

%──────────────────────────────────────────────────────────────────────────
%  PLOT — frequency sweep: group mean +/- SEM of per-animal means
%──────────────────────────────────────────────────────────────────────────
cmap = [0.5 0.5 0.5;
        0, 0, 1;
        1, 0.8, 0];    

baseLab = {'E'' (kPa)', 'E'''' (kPa)', 'tan \delta'};

for m = 1:3
    yl = baseLab{m};
    if LOG_Y(m); yl = ['log_{10} ' baseLab{m}]; end

    figure('Color','w','Position',[50+300*(m-1) 150 560 440]); hold on;

    for c = 1:nCond
        cMask   = strcmp(C_all{m}, condOrder{c});
        animals = unique(A_all{m}(cMask));
        nA      = numel(animals);

        anMat = NaN(nA, nFreq);
        for ai = 1:nA
            for f = 1:nFreq
                mask = cMask & A_all{m}==animals(ai) & abs(F_all{m}-freqLog(f))<1e-9;
                v = Y_all{m}(mask);
                if ~isempty(v); anMat(ai,f) = mean(v,'omitnan'); end
            end
        end

        gMean = mean(anMat, 1, 'omitnan');
        gSEM  = std(anMat, 0, 1, 'omitnan') / sqrt(nA);

        errorbar(frequencies, gMean, gSEM, '-o', ...
            'Color', cmap(c,:), 'LineWidth', 2, ...
            'MarkerFaceColor', cmap(c,:), 'MarkerSize', 7, ...
            'CapSize', 6, 'DisplayName', condOrder{c});
    end

    set(gca,'XScale','log','FontSize',12,'Box','off','TickDir','out')
    xticks(frequencies)
    xticklabels(arrayfun(@(x) sprintf('%g Hz',x), frequencies, 'UniformOutput', false))
    xlabel('Frequency (Hz)', 'FontSize', 13)
    ylabel(yl, 'FontSize', 13)
    title(yl, 'FontSize', 13)
    legend('Location','best','FontSize',10,'Box','off')
    hold off
end