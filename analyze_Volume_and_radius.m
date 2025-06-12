%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This beneath code was refined and proof-checked with the use of AI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
%% ________ Volume ANALYSIS_________

function analyze_spheroid_volume()
    % Analyze spheroid volumes under different cisplatin doses and mCherry status
    clear; clc; close all;

    %% 1) Load data
    fullPath = 'D:\DTU og KU\OneDrive - Danmarks Tekniske Universitet\Lab\volumeDivided.xlsx';
    T = readtable(fullPath, ...
                  'Sheet', 1, ...
                  'ReadVariableNames', true, ...
                  'VariableNamingRule', 'modify');

    %% 2) Verify variable names
    disp('Variable names in T:');
    disp(T.Properties.VariableNames);

    %% 3) Convert grouping variables to categorical
    T.Histogel     = categorical(T.Histogel,    {'no','yes'});
    T.mCherry      = categorical(T.mCherry,     {'no','yes'});
    T.Cisplatin_uM = categorical(T.Cisplatin_uM);

    %% 4) Preview first rows
    disp('First 8 rows of T:');
    disp(head(T,8));

    %% 5) Three-way ANOVA with interactions
    fprintf('\n=== Three-way ANOVA ===\n');
    [p, tbl, stats] = anovan( ...
        T.VolumeCorrected_mm3, ...
        {T.Histogel, T.mCherry, T.Cisplatin_uM}, ...
        'model',    'interaction', ...
        'varnames', {'Histogel','mCherry','Cisplatin_uM'}, ...
        'display',  'on' );

    %% 6) Post-hoc comparisons for cisplatin doses
    figPH = figure('Name','Post-hoc: Cisplatin Levels');
    multcompare(stats, 'Dimension', 3);
    title('Pairwise Comparisons of Cisplatin µM (Dunn–Sidak)');

    %% 7) Overall boxplot for all factor combinations
    figAll = figure('Name','Overall Boxplot');
    groupVars = strcat(string(T.Histogel), '_', string(T.mCherry), '_', string(T.Cisplatin_uM));
    boxplot(T.VolumeCorrected_mm3, groupVars, ...
            'LabelOrientation','inline', 'Whisker',1.5);
    xtickangle(45);
    ylabel('Volume (mm^3)');
    title('Spheroid Volume by Histogel, mCherry & Cisplatin µM');
    grid on;

    %% 8) Dose-specific boxplots (0, 25, 100 µM)
    doses = [0 25 100];
    for i = 1:numel(doses)
        d = doses(i);
        idx = T.Cisplatin_uM == categorical(string(d));
        labels = strcat(string(T.Histogel(idx)), '_', string(T.mCherry(idx)));

        figDose = figure('Name', sprintf('Volume at %d µM cisplatin', d));
        boxplot(T.VolumeCorrected_mm3(idx), labels, ...
                'LabelOrientation','inline', 'Whisker',1.5);
        xtickangle(45);
        ylabel('Volume (mm^3)');
        title(sprintf('Spheroid Volume at %d µM cisplatin', d));
        grid on;
    end

    %% 9) Create & save 6 boxplots: 0, 25, 100 µM × without/with mCherry
    mCherryFlags  = {'no','yes'};
    mCherryLabels = {'without mCherry','with mCherry'};

    for j = 1:2
        flag  = mCherryFlags{j};
        label = mCherryLabels{j};

        for i = 1:numel(doses)
            d = doses(i);

            % Filter for this dose and mCherry status
            idx = ( string(T.Cisplatin_uM) == string(d) ) & ...
              ( string(T.mCherry)      == flag );

            % --- REMAP histogel yes/no to your custom x-labels ---
            raw = string(T.Histogel(idx));                % ["yes","no",…]
            histoLabels = repmat("", numel(raw), 1);
            histoLabels(raw=="yes") = "Histogel";
            histoLabels(raw=="no")  = "No histogel";

            % Plot
            fh = figure('Name', sprintf('%d µM cisplatin, %s', d, label) );
            boxplot(T.VolumeCorrected_mm3(idx), histoLabels, ...
               'Whisker',1.5, 'LabelOrientation','inline');
            xtickangle(45);
            ylabel('Volume (mm^3)');
            title( sprintf('Volume at %d µM cisplatin, %s', d, label) );

             % **Force same y-axis range on every plot**
            ylim([0 0.14]);

            grid on;

            % Save as PNG
            safeLabel = strrep(label, ' ', '_');  % e.g. without_mCherry
            filename = sprintf('Boxplot_%dum_cisplatin_%s.png', d, safeLabel);
            saveas(fh, filename);
        end
    end

end

%% ________ RADIUS ANALYSIS_________

%% 1) Read data from Excel
filename = 'RadiusOpdelt.xlsx';
T = readtable(filename);

% After this, the table T should look like:
%    Condition                     Dose    RadiusAtt    RadiusPix
%    ____________________________  ____   __________   __________
%    'With histogel, without ...'   0      0.210        0.129   
%    'With histogel, without ...'  25      0.225        0.134   
%    'With histogel, without ...' 100      0.240        0.123   
%    'With histogel, with ...'      0      0.221        0.149   
%    'With histogel, with ...'     25      0.228        0.157   
%    ... etc.

% Check the class of RadiusAtt og RadiusPix
if iscell(T.RadiusAtt)
    % Convert cell‐array of strings (evt. med komma) til numerisk vektor
    % Hvis dine tal i Excel bruger komma som decimaltegn, kan det være nødvendigt
    % at fjerne punktummer til korttusindtalsseparator og erstatte komma med punktum:
    %
    %   T.RadiusAtt = strrep(T.RadiusAtt, '.', '');  % fjern evt. punktummer
    %   T.RadiusAtt = strrep(T.RadiusAtt, ',', '.'); % erstat komma med punktum
    %
    T.RadiusAtt = str2double(T.RadiusAtt);
end

if iscell(T.RadiusPix)
    % Samme procedure for RadiusPix
    %
    %   T.RadiusPix = strrep(T.RadiusPix, '.', '');
    %   T.RadiusPix = strrep(T.RadiusPix, ',', '.');
    %
    T.RadiusPix = str2double(T.RadiusPix);
end

%% 2) Calculate absolute and percentage differences
T.Diff_abs_mm = T.RadiusAtt - T.RadiusPix;
T.Diff_pct    = 100 * (T.Diff_abs_mm ./ T.RadiusPix);

% Display a preview
disp(T(1:min(10,height(T)),:));

%% 3) Scatter plot: RadiusPix vs. RadiusAtt
figure('Name','Scatter: 2D-radius vs. Attenuation-radius','NumberTitle','off');
scatter(T.RadiusPix, T.RadiusAtt, 'filled');
hold on;
maxVal = max([T.RadiusPix; T.RadiusAtt]) * 1.1;
plot([0, maxVal], [0, maxVal], 'r--', 'LineWidth', 1.2);
hold off;
xlabel('Radius from 2D pixels [mm]');
ylabel('Radius from Attenuation [mm]');
title('Scatter Plot: 2D-radius vs. Attenuation-radius');
legend({'Data points','y = x'}, 'Location', 'best');
grid on;

%% 5) Save the updated table back to Excel
writetable(T, 'RadiusOpdelt_WithDiff.xlsx');
disp('Saved updated table with Diff_abs_mm and Diff_pct as RadiusOpdelt_WithDiff.xlsx');


   

  