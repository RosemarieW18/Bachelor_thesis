function analyze2_spheroid_attenuation()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This beneath code was refined and proof-checked with the use of AI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear; clc; close all;

    %% 1) Load attenuation sheet
    fname = 'Attenuation_Opdelt.xlsx';
    sheetName = 'Sheet1';
    T = readtable(fname, 'Sheet', sheetName);

    %% 2) Print all column names for verification
    fprintf('\nColumn names in the table:\n');
    disp(T.Properties.VariableNames);

    % If your column is not exactly "MeanMuClamped", change the string below
    colName = 'MeanMuClamped';

    %% 3) Preserve numeric Cisplatin and convert grouping variables
    T.Cisplatin_uM_num = T.Cisplatin_uM;   % numeric (0, 25, 100)
    T.Cisplatin_uM     = categorical(T.Cisplatin_uM);
    T.Histogel         = categorical(T.Histogel, {'no','yes'});
    T.mCherry          = categorical(T.mCherry,  {'no','yes'});

    %% 4) Display first few rows
    fprintf('\nFirst 8 rows of the attenuation table:\n');
    disp(head(T, 8));

    %% 5) One-way ANOVA: “No Histogel” vs. “With Histogel” (all data)
    fprintf('\n=== One-way ANOVA: Comparison of Histogel and No Histogel attenuation values ===\n');
    dataVals = T.(colName);    % use dynamic field access for the attenuation column
    [p_hist, ~, ~] = anova1(dataVals, T.Histogel, 'off');
    fprintf('  p-value (Histogel vs. No Histogel) = %.4f\n', p_hist);

    % Create a temporary categorical with reversed order: {'yes','no'}
    orderHG = categorical(T.Histogel, {'yes','no'});

    fh_hist = figure('Name','Comparison of Histogel and No Histogel');
    boxplot(dataVals, orderHG, 'Whisker', 1.5);
    ylabel('Clamped μ (1/mm)');
    xlabel('Embedding');
    title(sprintf('Comparison of Histogel and No Histogel attenuation values (p = %.3f)', p_hist));
    set(gca, 'XTickLabel', {'Histogel','No Histogel'});
    ylim([0 5]);  % Force Y-axis from 0 to 5
    grid on;
    saveas(fh_hist, 'Boxplot_Histogel_vs_NoHistogel.png');

    %% 6) One-way ANOVA: Cisplatin effect within each Histogel group
    % No Histogel
    idx_no = (T.Histogel == 'no');
    if any(idx_no)
        fprintf('\n=== One-way ANOVA: Cisplatin effect in NO Histogel ===\n');
        uniqueDoses_no = unique(T.Cisplatin_uM_num(idx_no));
        doseCat_no = categorical(T.Cisplatin_uM_num(idx_no), uniqueDoses_no, ...
                                 strcat(string(uniqueDoses_no), "µM"));
        [p_no, ~, stats_no] = anova1(dataVals(idx_no), doseCat_no, 'off');
        fprintf('  p-value (NO Histogel) = %.4f\n', p_no);

        fh_no = figure('Name','Boxplot: No Histogel, μ vs. Cisplatin');
        boxplot(dataVals(idx_no), doseCat_no, 'Whisker', 1.5);
        xlabel('Cisplatin (µM)');
        ylabel('Clamped μ (1/mm)');
        title(sprintf('NO Histogel: μ vs. Cisplatin (p = %.3f)', p_no));
        ylim([0 5]);  % Force Y-axis 0–5
        grid on;
        saveas(fh_no, 'Boxplot_NoHistogel_Cisplatin.png');

        if p_no < 0.05
            fh_post_no = figure('Name','Post-hoc NO Histogel (Tukey-Kramer)');
            multcompare(stats_no, 'CType', 'tukey-kramer');
            title('Post-hoc (NO Histogel): Cisplatin Doses');
            saveas(fh_post_no, 'Posthoc_NoHistogel_Cisplatin.png');
        end
    end

    % With Histogel
    idx_yes = (T.Histogel == 'yes');
    if any(idx_yes)
        fprintf('\n=== One-way ANOVA: Cisplatin effect in WITH Histogel ===\n');
        uniqueDoses_yes = unique(T.Cisplatin_uM_num(idx_yes));
        doseCat_yes = categorical(T.Cisplatin_uM_num(idx_yes), uniqueDoses_yes, ...
                                  strcat(string(uniqueDoses_yes), "µM"));
        [p_yes, ~, stats_yes] = anova1(dataVals(idx_yes), doseCat_yes, 'off');
        fprintf('  p-value (WITH Histogel) = %.4f\n', p_yes);

        fh_yes = figure('Name','Boxplot: With Histogel, μ vs. Cisplatin');
        boxplot(dataVals(idx_yes), doseCat_yes, 'Whisker', 1.5);
        xlabel('Cisplatin (µM)');
        ylabel('Clamped μ (1/mm)');
        title(sprintf('WITH Histogel: μ vs. Cisplatin (p = %.3f)', p_yes));
        ylim([0 5]);  % Force Y-axis 0–5
        grid on;
        saveas(fh_yes, 'Boxplot_WithHistogel_Cisplatin.png');

        if p_yes < 0.05
            fh_post_yes = figure('Name','Post-hoc WITH Histogel (Tukey-Kramer)');
            multcompare(stats_yes, 'CType', 'tukey-kramer');
            title('Post-hoc (WITH Histogel): Cisplatin Doses');
            saveas(fh_post_yes, 'Posthoc_WithHistogel_Cisplatin.png');
        end
    end

    %% 7) One-way ANOVA: mCherry effect (compare Histogel vs. No Histogel)
    fprintf('\n=== One-way ANOVA: mCherry effect (separate for mCherry=no and mCherry=yes) ===\n');

    % Without mCherry
    idx_no_mch = (T.mCherry == 'no');
    if any(idx_no_mch)
        fprintf('  mCherry = no: Comparing Histogel vs. No Histogel\n');
        subT = T(idx_no_mch,:);
        data_no_mch = subT.(colName);
        [p_nm, ~, ~] = anova1(data_no_mch, subT.Histogel, 'off');
        fprintf('    p-value (no mCherry) = %.4f\n', p_nm);

        % Create categorical with order {'yes','no'} so Histogel is left
        embedCat_no_mch = categorical(subT.Histogel, {'yes','no'});

        fh_no_mch = figure('Name','Clamped Attenuation without mCherry');
        boxplot(data_no_mch, embedCat_no_mch, 'Whisker', 1.5);
        ylabel('Clamped μ (1/mm)');
        xlabel('Embedding');
        title('Clamped Attenuation without mCherry');
        set(gca, 'XTickLabel', {'Histogel','No Histogel'});
        ylim([0 5]);
        grid on;
        saveas(fh_no_mch, 'Boxplot_ClampedAtten_without_mCherry.png');
    end

    % With mCherry
    idx_yes_mch = (T.mCherry == 'yes');
    if any(idx_yes_mch)
        fprintf('  mCherry = yes: Comparing Histogel vs. No Histogel\n');
        subT = T(idx_yes_mch,:);
        data_yes_mch = subT.(colName);
        [p_ym, ~, ~] = anova1(data_yes_mch, subT.Histogel, 'off');
        fprintf('    p-value (with mCherry) = %.4f\n', p_ym);

        % Create categorical with order {'yes','no'} so Histogel is left
        embedCat_yes_mch = categorical(subT.Histogel, {'yes','no'});

        fh_yes_mch = figure('Name','Clamped Attenuation with mCherry');
        boxplot(data_yes_mch, embedCat_yes_mch, 'Whisker', 1.5);
        ylabel('Clamped μ (1/mm)');
        xlabel('Embedding');
        title('Clamped Attenuation with mCherry');
        set(gca, 'XTickLabel', {'Histogel','No Histogel'});
        ylim([0 5]);
        grid on;
        saveas(fh_yes_mch, 'Boxplot_ClampedAtten_with_mCherry.png');
    end

    %% 8) Two-way ANOVA on clamped attenuation (Histogel × mCherry)
    fprintf('\n=== Two-way ANOVA on Clamped Attenuation (μ ≥ 0) ===\n');
    [p2, tbl2, stats2] = anovan( ...
        dataVals, ...
        { T.Histogel, T.mCherry }, ...
        'model',    'interaction', ...
        'varnames', {'Histogel','mCherry'}, ...
        'display',  'on' );

    %Compute eta-squared effect sizes
    % Sum of Squares for each row in the ANOVA table:
    ss_hist2  = tbl2{2,2};    % Histogel SS
    ss_mch2   = tbl2{3,2};    % mCherry SS
    ss_int2   = tbl2{4,2};    % Interaction SS
    ss_tot2   = tbl2{end,2};  % Total SS
    eta2_hist2 = ss_hist2 / ss_tot2;
    eta2_mch2  = ss_mch2  / ss_tot2;
    eta2_int2  = ss_int2  / ss_tot2;
    fprintf('\n=== Effect sizes (η²) ===\n');
    fprintf('  Histogel η²        = %.3f\n', eta2_hist2);
    fprintf('  mCherry η²         = %.3f\n', eta2_mch2);
    fprintf('  Interaction η²     = %.3f\n', eta2_int2);

    fh_post2 = figure('Name','Post-hoc: Histogel × mCherry');
    multcompare(stats2, 'CType', 'tukey-kramer');
    title('Pairwise Comparisons: Histogel × mCherry (Tukey–Kramer)');
    saveas(fh_post2, 'Posthoc_Histogel_mCherry.png');

    fh_box2 = figure('Name','Attenuation by Histogel & mCherry');
    groupLabels2 = strcat(string(T.Histogel), '_', string(T.mCherry));
    boxplot(dataVals, groupLabels2, 'LabelOrientation','inline', 'Whisker', 1.5);
    xtickangle(45);
    ylabel('Clamped μ (1/mm)');
    title('Clamped Attenuation by Histogel and mCherry');
    set(gca, 'XTickLabel', {'Histogel_no','Histogel_yes','NoHistogel_no','NoHistogel_yes'});
    ylim([0 5]);  % Force Y-axis 0–5
    grid on;
    saveas(fh_box2, 'Boxplot_Histogel_mCherry.png');

    %% 9) Correlation: attenuation vs. volume (Pearson)
    %if ismember('VolumeCorrected_mm3', T.Properties.VariableNames)
        %fprintf('\n=== Correlation: Attenuation (μ) vs. Volume (mm^3) ===\n');
        %[r_p, p_p] = corr(dataVals, T.VolumeCorrected_mm3, ...
        %                  'Type', 'Pearson', 'Rows', 'complete');
        %fprintf('  Pearson r = %.3f, p = %.3f\n', r_p, p_p);

        %fh_scatt = figure('Name','Volume vs Attenuation');
        %scatter(dataVals, T.VolumeCorrected_mm3, 'filled');
        %xlabel('Clamped μ (1/mm)');
        %ylabel('VolumeCorrected_mm^3');
        %title(sprintf('Volume vs Attenuation (r = %.2f, p = %.3f)', r_p, p_p));
        %ylim([0 5]);  % keep Y-axis consistent
        %grid on;
        %saveas(fh_scatt, 'Scatter_Volume_Attenuation.png');
    %else
        %warning('Column ''VolumeCorrected_mm3'' not found. Correlation skipped.');
    %end

    fprintf('\n==== All analyses complete ====\n');
end
