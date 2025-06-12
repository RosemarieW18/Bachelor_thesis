%% LIGHT SHEET mCherry ANALYSIS IN MATLAB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This beneath code was refined and proof-checked with the use of AI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1) File path setup
%-----------------------------------------------------------------------------------%
baseDir    = 'D:\DTU og KU\OneDrive - Danmarks Tekniske Universitet\Lab';
subFolder  = 'Light sheet mCherry';
dataFolder = fullfile(baseDir, subFolder);
excelFile  = 'mCherry_summary.xlsx';
excelPath  = fullfile(dataFolder, excelFile);

%% 2) Import data, preserving original column names
%-----------------------------------------------------------------------------------%
opts = detectImportOptions(excelPath, 'VariableNamingRule','preserve');
T    = readtable(excelPath, opts);

% Verify the variable names
disp('Table variable names:');
disp(T.Properties.VariableNames');

%% 3) Extract relevant variables
%-----------------------------------------------------------------------------------%
fileNames  = string(T.Filename);            % image filenames
totalFluor = T.TotalFluorescence;           % total fluorescence per image
meanTop5   = T.MeanTop5;                    % mean of top 5% brightest pixels

%% 4) Assign concentration groups based on filename prefix
%-----------------------------------------------------------------------------------%
n    = height(T);
conc = nan(n,1);
for i = 1:n
    fn = fileNames(i);
    if startsWith(fn, "2B") || startsWith(fn, "4B")
        conc(i) = 0;
    elseif startsWith(fn, "7B")
        conc(i) = 25;
    elseif startsWith(fn, "9G") || startsWith(fn, "11G")
        conc(i) = 100;
    else
        warning('Unexpected filename: %s', fn);
    end
end

%% 5) Boxplot: MeanTop5 by concentration
%-----------------------------------------------------------------------------------%
fig1 = figure;
boxplot(meanTop5, conc, 'Labels', {'0 µM','25 µM','100 µM'});
xlabel('Cisplatin concentration (µM)');
ylabel('Mean of top 5% pixel intensities');
title('mCherry Top-5% Intensity by Cisplatin Dose');

% Save figure 1
saveas(fig1, fullfile(dataFolder, 'Boxplot_MeanTop5.png'));

%% 6) Boxplot: TotalFluorescence by concentration
%-----------------------------------------------------------------------------------%
fig2 = figure;
boxplot(totalFluor, conc, 'Labels', {'0 µM','25 µM','100 µM'});
xlabel('Cisplatin concentration (µM)');
ylabel('Total fluorescence (sum of all pixels)');
title('Total mCherry Fluorescence by Cisplatin Dose');

% Save figure 2
saveas(fig2, fullfile(dataFolder, 'Boxplot_TotalFluorescence.png'));

%% 7) One-way ANOVA on MeanTop5
%-----------------------------------------------------------------------------------%
[pA1, ~, statsA1] = anova1(meanTop5, conc, 'off');
fprintf('ANOVA on MeanTop5: p = %.4f\n', pA1);
if pA1 < 0.05
    multcompare(statsA1, 'CType', 'tukey-kramer');
end

%% 8) Kruskal–Wallis on MeanTop5 (do no use)
%-----------------------------------------------------------------------------------%
%[pK1, ~, statsK1] = kruskalwallis(meanTop5, conc, 'off');
%fprintf('Kruskal–Wallis on MeanTop5: p = %.4f\n', pK1);
%if pK1 < 0.05
%    multcompare(statsK1, 'CType', 'dunn-sidak');
%end

%% 9) One-way ANOVA on TotalFluorescence
%-----------------------------------------------------------------------------------%
[pA2, ~, statsA2] = anova1(totalFluor, conc, 'off');
fprintf('ANOVA on TotalFluorescence: p = %.4f\n', pA2);
if pA2 < 0.05
    multcompare(statsA2, 'CType', 'tukey-kramer');
end

%% 10) Kruskal–Wallis on TotalFluorescence
%-----------------------------------------------------------------------------------%
[pK2, ~, statsK2] = kruskalwallis(totalFluor, conc, 'off');
fprintf('Kruskal–Wallis on TotalFluorescence: p = %.4f\n', pK2);
if pK2 < 0.05
    multcompare(statsK2, 'CType', 'dunn-sidak');
end

%% 11) Simple linear regression (TotalFluorescence vs. concentration)
%-----------------------------------------------------------------------------------%
% Fit model: totalFluor = intercept + slope*conc
X        = [ones(n,1), conc];
[b,bint,~,~,stats] = regress(totalFluor, X);
slope    = b(2);
intercept= b(1);
R2       = stats(1);
p_slope  = stats(3);

fprintf('\nLinear regression on TotalFluorescence:\n');
fprintf('  slope = %.4f units per µM\n', slope);
fprintf('  intercept = %.1f\n', intercept);
fprintf('  R² = %.3f, p (slope) = %.4f\n', R2, p_slope);

% Overlay the fit line on the same boxplot
figure; 
boxplot(totalFluor,conc,'Labels',{'0','25','100'}); hold on;
xline = linspace(0,100,100)';
yline = intercept + slope*xline;
plot(xline,yline,'r-','LineWidth',2);
xlabel('Cisplatin (µM)'); ylabel('Total Fluorescence');
title('Total Fluorescence with Linear Fit');
legend('Data','Linear fit','Location','NorthEast');
saveas(gcf,fullfile(dataFolder,'TotalFluor_vs_Concentration_Regression.png'));

disp('Analysis complete. Figures saved to:');
disp(dataFolder);

