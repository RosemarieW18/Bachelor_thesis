%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This beneath code was refined and proof-checked with the use of AI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0) Opsætning

filename = 'Statistik_data_frekvens.xlsx';
sheets   = {'Ark1','Ark2','Ark3','Ark4'};

% Navne på de fire grupper i samme rækkefølge som arkene
groupNames = { ...
   'NoHistogel', ...            
   'Histogel', ...              
   'NoHistogel_mCherry', ...    
   'Histogel_mCherry' ...       
};

% 1) Læs alle fire ark ind i én samlet tabel T
T = table();
for i = 1:4
    Ti = readtable(filename, 'Sheet', sheets{i});
    Ti.Group = repmat(string(groupNames{i}), height(Ti), 1);
    T = [T; Ti];
end

% 2) Konverter til korrekte typer
T.Well      = string(T.Well);
T.Condition = categorical(T.Condition, [0 25 100], {'0','25','100'});
T.Band      = categorical(T.Band,      {'Low','Mid','High'});
T.Group     = categorical(T.Group);

% 3) Filtrer kun mCherry-grupperne ud
T12 = T(ismember(T.Group, {'NoHistogel_mCherry','Histogel_mCherry'}), :);

% 2) Filtrér kun de to relevante mCherry-grupper
T12 = T(ismember(T.Group,{'NoHistogel_mCherry','Histogel_mCherry'}),:);
T12.Group     = removecats(T12.Group);
T12.Condition = removecats(T12.Condition);
T12.Band      = removecats(T12.Band);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%% HYPOTESE 2 %%%%%%%
   % Histogel_mCherry vs. NoHistogel_mCherry (absolut energy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Three-way ANOVA med 3-vejs-interaktion
[p, tbl, stats] = anovan( ...
    T12.Absolut_energy, ...
    {T12.Group, T12.Condition, T12.Band}, ...
    'model',    'full', ...              % inkl. Group×Condition×Band
    'varnames', {'Group','Condition','Band'}, ...
    'alpha',    0.05, ...
    'display',  'off' );

%2) Formatér p-værdier til 4 decimaler
tbl2     = tbl;                         
nEffects = numel(p);                    
for k = 1:nEffects
    tbl2{k+1,6} = sprintf('%.4f', p(k)); 
end

%3) Opret figur med uitable
f = figure( ...
    'Name',        'Three-way ANOVA Results', ...
    'NumberTitle', 'off', ...
    'Color',       [1 1 1], ...
    'Position',    [200 200 800 300] );

% Hent header og data fra tbl2
colNames  = tbl2(1,2:end);    % {'DF','Sum Sq','Mean Sq','F','pValue'}
rowNames  = tbl2(2:end,1);    % {'Group','Condition',…,'Group:Cond:Band','Error','Total'}
tableData = tbl2(2:end,2:end);

% Vis tabellen i figure
uitable( ...
    'Parent',     f, ...
    'Data',       tableData, ...
    'ColumnName', colNames, ...
    'RowName',    rowNames, ...
    'Units',      'normalized', ...
    'Position',   [0 0 1 1], ...
    'FontSize',   11 );


%% SUBPLOT frequency energy 
% 10) Opsætning og indlæsning
filename   = 'Statistik_data_frekvens.xlsx';
sheets     = {'Ark1','Ark2','Ark3','Ark4'};
groupNames = { ...
   'NoHistogel', ...            
   'Histogel', ...              
   'NoHistogel_mCherry', ...    
   'Histogel_mCherry' ...       
};

% Læs alle fire ark ind i én tabel T
T = table();
for i = 1:4
    Ti = readtable(filename,'Sheet',sheets{i});
    Ti.Group = repmat(string(groupNames{i}), height(Ti), 1);
    T = [T; Ti];
end

% Konverter Well til string
T.Well = string(T.Well);

% Gør ’Condition’, ’Band’ og ’Group’ til kategoriske variable
T.Condition = categorical(T.Condition, [0 25 100], {'0','25','100'});
T.Band      = categorical(T.Band,      {'Low','Mid','High'});
T.Group     = categorical(T.Group);

% 1) Definér subsets (“NoHistogel” og “Histogel”)
T_no   = T(T.Group == 'NoHistogel_mCherry', :);
T_hist = T(T.Group == 'Histogel_mCherry',   :);

% 2) Beregn fælles y-akse-grænser for absolut energy
allValues = [T_no.Absolut_energy; T_hist.Absolut_energy];
ymin = floor(min(allValues)) - 1;
ymax = ceil( max(allValues)) + 1;

% Plot bokse side-by-side med samme y-akse
figure('Position',[100 100 1200 400]);

% ——————————————————————————————————————————
% MED Histogel
subplot(1,2,1);
boxplot( ...
    T_hist.Absolut_energy, ...
    {T_hist.Condition, T_hist.Band}, ...
    'FactorSeparator',1, ...
    'LabelVerbosity','minor', ...
    'Whisker',1.5 ...
);
xlabel('Cisplatin condition (µM) and band');
ylabel('Absolut\_energy');
title('Histogel\_mCherry');
set(gca,'XTickLabelRotation',45);
ylim([ymin ymax]);

% ——————————————————————————————————————————
% UDEN NoHistogel
subplot(1,2,2);
boxplot( ...
    T_no.Absolut_energy, ...
    {T_no.Condition, T_no.Band}, ...
    'FactorSeparator',1, ...
    'LabelVerbosity','minor', ...
    'Whisker',1.5 ...
);
xlabel('Cisplatin condition (µM) and band');
ylabel('Absolut\_energy');
title('NoHistogel\_mCherry');
set(gca,'XTickLabelRotation',45);
ylim([ymin ymax]);


%% MEAN beregning 
% 1) Filtrér kun de to relevante
T12 = T( ismember(T.Group, {'NoHistogel_mCherry','Histogel_mCherry'}) , : );
T12.Group     = removecats(T12.Group);
T12.Condition = removecats(T12.Condition);
T12.Band      = removecats(T12.Band);

% 2) Tre‐vejs ANOVA (hvis ikke allerede gjort)
[p, tbl, stats] = anovan( ...
    T12.Absolut_energy, ...
    {T12.Group, T12.Condition, T12.Band}, ...
    'model','interaction', ...
    'varnames',{'Group','Condition','Band'}, ...
    'alpha',0.05, ...
    'display','off' ...
);

% 3) Beregn ’celle‐midler’ for hver kombination af Condition × Band (mCherry)
conds = categories(T12.Condition);
bands = categories(T12.Band);
nCond = numel(conds);
nBand = numel(bands);

cellMean = nan(nCond, nBand);
for i = 1:nCond
    for j = 1:nBand
        selCell = (T12.Condition == conds{i}) & (T12.Band == bands{j});
        cellMean(i,j) = mean(T12.Absolut_energy(selCell));
    end
end

% 4) Beregn residualer
N = height(T12);
residuals = nan(N,1);
for k = 1:N
    i = find(strcmp(conds, char(T12.Condition(k))));
    j = find(strcmp(bands,  char(T12.Band(k))));
    residuals(k) = T12.Absolut_energy(k) - cellMean(i,j);
end

% Tilføj residualer som ny kolonne
T12.Residual = residuals;

% 5) Beregn grand mean
grandMean = mean(T12.Absolut_energy);

% 6) Beregn gruppespecifikke residual means
mean_no_res   = mean(T12.Residual(T12.Group == 'NoHistogel_mCherry'));
mean_hist_res = mean(T12.Residual(T12.Group == 'Histogel_mCherry'));

% 7) Justerede means
adjMean_no   = grandMean + mean_no_res;
adjMean_hist = grandMean + mean_hist_res;

% 8) Rå middelværdier for hver mCherry-gruppe
rawMean_no   = mean(T12.Absolut_energy(T12.Group == 'NoHistogel_mCherry'));
rawMean_hist = mean(T12.Absolut_energy(T12.Group == 'Histogel_mCherry'));

% 9) Opret en tabel med både rå means og justerede means
GroupNames = {'Histogel\_mCherry'; 'NoHistogel\_mCherry'};
RawMean    = [rawMean_hist;   rawMean_no];
AdjMean    = [adjMean_hist;   adjMean_no];

SummaryTbl = table(GroupNames, RawMean, AdjMean, ...
    'VariableNames',{'Group','RawMean_AbsolutEnergy','AdjMean_AbsolutEnergy'});

disp('=== Sammenligning af Rå Mean vs. Justeret Mean for mCherry Groups ===');
disp(SummaryTbl);

% 10) Konklusion på justerede means
if adjMean_hist > adjMean_no
    fprintf('\nKonklusion: HISTOGEL\_mCherry har størst justeret absolut energy (%.3e > %.3e).\n', ...
            adjMean_hist, adjMean_no);
else
    fprintf('\nKonklusion: NOHISTOGEL\_mCherry har størst justeret absolut energy (%.3e > %.3e).\n', ...
            adjMean_no, adjMean_hist);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%% HYPOTESE 3 %%%%%%%
       % Histogel_mCherry vs. NoHistogel_mCherry (percentage energy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SUBPLOT frequency energy (histogel vs no histogel)

% 2) Filtrér kun de to relevante grupper
T12 = T(ismember(T.Group, {'NoHistogel_mCherry','Histogel_mCherry'}), :);
T12.Group     = removecats(T12.Group);
T12.Condition = removecats(T12.Condition);
T12.Band      = removecats(T12.Band);

% 10) Opsætning og indlæsning
filename   = 'Statistik_data_frekvens.xlsx';
sheets     = {'Ark1','Ark2','Ark3','Ark4'};
groupNames = { ...
   'NoHistogel', ...
   'Histogel', ...
   'NoHistogel_mCherry', ...
   'Histogel_mCherry' ...
};

% Læs alle fire ark ind i én tabel T
T = table();
for i = 1:4
    Ti = readtable(filename, 'Sheet', sheets{i});
    Ti.Group = repmat(string(groupNames{i}), height(Ti), 1);
    T = [T; Ti];
end

% Konverter Well til string (giver nemmere sammenligning senere)
T.Well = string(T.Well);

% Gør ’Condition’, ’Band’ og ’Group’ til kategoriske variable
T.Condition = categorical(T.Condition, [0 25 100], {'0','25','100'});
T.Band      = categorical(T.Band,      {'Low','Mid','High'});
T.Group     = categorical(T.Group);

% 1) Definér subsets (“NoHistogel_mCherry” og “Histogel_mCherry”)
T_no   = T(T.Group == 'NoHistogel_mCherry', :);
T_hist = T(T.Group == 'Histogel_mCherry',   :);

% 2) Beregn fælles y-akse-grænser
allValues = [T_no.Percentage_energy; T_hist.Percentage_energy];
ymin = floor(min(allValues)) - 1;   % små ’pad’ for fri luft
ymax = ceil( max(allValues)) + 1;

% Opret figur og størrelsen på vinduet
figure('Position', [100 100 1200 400]);

% Vi gemmer akse-handles i et array, så vi kan flytte dem bagefter
axs = gobjects(1,2);

% ——————————————————————————————————————————
% MED histogel (venstre)
axs(1) = subplot(1,2,1);
boxplot( ...
    T_hist.Percentage_energy, ...
    {T_hist.Condition, T_hist.Band}, ...
    'FactorSeparator', 1, ...
    'LabelVerbosity', 'minor', ...
    'Whisker', 1.5 ...
);
xlabel('Cisplatin concentration (µM) and band');
ylabel('Percentage\_energy (%)');
title('Histogel_mCherry');
set(gca, 'XTickLabelRotation', 45);
ylim([0 70]);
%ylim([ymin ymax]);  % sætter samme y-akse

% ——————————————————————————————————————————
% UDEN histogel (højre)
axs(2) = subplot(1,2,2);
boxplot( ...
    T_no.Percentage_energy, ...
    {T_no.Condition, T_no.Band}, ...
    'FactorSeparator', 1, ...
    'LabelVerbosity', 'minor', ...
    'Whisker', 1.5 ...
);
xlabel('Cisplatin concentration (µM) and band');
ylabel('Percentage\_energy (%)');
title('NoHistogel_mCherry');
set(gca, 'XTickLabelRotation', 45);
ylim([0 70]);
%ylim([ymin ymax]);  % sætter samme y-akse

% ——————————————————————————————————————————
% Flyt subplot-akserne nedad for at give plads til sgtitle
for k = 1:2
    posA = axs(k).Position;   % posA = [left, bottom, width, height]
    posA(2) = posA(2) - 0.05; % træk bunden 5% nedad
    posA(4) = posA(4) - 0.03; % reducér aksens højde en smule
    axs(k).Position = posA;
end

% Tilføj fælles overordnet titel (sgtitle) og løft den
hSG = sgtitle( ...
    'Spectral Energy by Band vs. Cisplatin (Histogel\_mCherry vs. NoHistogel\_mCherry)', ...
    'FontSize', 14, 'FontWeight', 'bold' ...
);
hSG.Units = 'normalized';
posSG = hSG.Position;       % typisk [0.5 1.0 0]
posSG(2) = posSG(2) + 0.04;  % løft sgtitle ca. 4% opad
hSG.Position = posSG;
hSG.HorizontalAlignment = 'center';

% (Hvis du har annotation som følger, kan du lade den stå eller fjerne den)
annotation('textbox', [0.35,0.93,0.3,0.05], ...
    'String', '', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none' ...
);


%% 1) Three-way ANOVA med 3-vejs-interaktion
[p, tbl, stats] = anovan( ...
    T12.Percentage_energy, ...
    {T12.Group, T12.Condition, T12.Band}, ...
    'model',    'full', ...              % inkl. Group×Condition×Band
    'varnames', {'Group','Condition','Band'}, ...
    'alpha',    0.05, ...
    'display',  'off' );

%2) Formatér p-værdier til 4 decimaler
tbl2     = tbl;                         
nEffects = numel(p);                    
for k = 1:nEffects
    tbl2{k+1,6} = sprintf('%.4f', p(k)); 
end

%3) Opret figur med uitable
f = figure( ...
    'Name',        'Three-way ANOVA Results', ...
    'NumberTitle', 'off', ...
    'Color',       [1 1 1], ...
    'Position',    [200 200 800 300] );

% Hent header og data fra tbl2
colNames  = tbl2(1,2:end);    % {'DF','Sum Sq','Mean Sq','F','pValue'}
rowNames  = tbl2(2:end,1);    % {'Group','Condition',…,'Group:Cond:Band','Error','Total'}
tableData = tbl2(2:end,2:end);

% Vis tabellen i figure
uitable( ...
    'Parent',     f, ...
    'Data',       tableData, ...
    'ColumnName', colNames, ...
    'RowName',    rowNames, ...
    'Units',      'normalized', ...
    'Position',   [0 0 1 1], ...
    'FontSize',   11 );

%% POST-Hoc Group×Band 
[~,~,stats] = anovan( ...
    T12.Absolut_energy, ...
    {T12.Group, T12.Condition, T12.Band}, ...
    'model',   'full', ...
    'varnames',{'Group','Condition','Band'}, ...
    'alpha',   0.05, ...
    'display', 'off' );

%2) Multcompare for Group×Band
[cgb, ~, ~, gNames] = multcompare( ...
    stats, ...
    'Dimension', [1 3], ...
    'CType',     'tukey-kramer', ...
    'Display',   'off' );

disp('gNames fra multcompare ser sådan ud:');
disp(gNames);

%3) Post-hoc: NoHistogel vs Histogel inden for hvert Band
bandLevels = categories(T12.Band);
fprintf('\nPost-hoc: NoHistogel (mcheery) vs Histogel (mcheery) inden for hvert Band:\n');

for bi = 1:numel(bandLevels)
    B = bandLevels{bi};

    % Find alle entries for dette band
    idx_band = find(contains(gNames, ['Band=' B]));
    if numel(idx_band) ~= 2
        warning('Kun fandt %d kombinationer for band "%s".', numel(idx_band), B);
        continue;
    end

    % Identificer hvilket index er NoHistogel og Histogel
    isHistogel    = contains(gNames(idx_band), 'Group=Histogel');
    idx_Histogel  = idx_band(isHistogel);
    idx_NoHisto   = idx_band(~isHistogel);

    % Find den række i cgb med netop disse to indeks
    row = find( (cgb(:,1)==idx_NoHisto & cgb(:,2)==idx_Histogel) | ...
                (cgb(:,1)==idx_Histogel & cgb(:,2)==idx_NoHisto) );
    if isempty(row)
        fprintf('  Band = %s: Ingen signifikant sammenligning fundet.\n', B);
        continue;
    end

    % Træk resultater ud
    CIlo = cgb(row,3);
    trekant    = cgb(row,4);
    CIhi = cgb(row,5);
    pval = cgb(row,6);

    % Udskriv
    fprintf('  Band = %-4s: Δ = %.3f (CI [%.3f, %.3f]), p = %.4f\n', ...
        B, trekant, CIlo, CIhi, pval);
end

%% MEDIAN 
% 1) Filtrér kun de to mCherry-grupper
T12 = T(ismember(T.Group, {'NoHistogel_mCherry','Histogel_mCherry'}), :);
T12.Group     = removecats(T12.Group);
T12.Condition = removecats(T12.Condition);
T12.Band      = removecats(T12.Band);

% 2) Hent alle bånd‐niveauer
bands = categories(T12.Band);   % {'Low','Mid','High'}
nBands = numel(bands);

% 3) Forbered vektorer til at gemme medianer
median_no   = nan(nBands,1);
median_hist = nan(nBands,1);

% 4) Loop over hvert bånd og beregn median for hver mCherry-gruppe
for k = 1:nBands
    b = bands{k};
    
    % --- NoHistogel_mCherry for dette bånd ---
    selNo = (T12.Group == 'NoHistogel_mCherry') & (T12.Band == b);
    dataNo = T12.Percentage_energy(selNo);
    if ~isempty(dataNo)
        median_no(k) = median(dataNo);
    else
        median_no(k) = NaN;
    end
    
    % --- Histogel_mCherry for dette bånd ---
    selHist = (T12.Group == 'Histogel_mCherry') & (T12.Band == b);
    dataHist = T12.Percentage_energy(selHist);
    if ~isempty(dataHist)
        median_hist(k) = median(dataHist);
    else
        median_hist(k) = NaN;
    end
end

% 5) Byg en tabel med median‐værdierne
BandNames = bands;
MedianNo   = median_no;
MedianHist = median_hist;

MedianTbl = table(BandNames, MedianNo, MedianHist, ...
    'VariableNames',{'Band','NoHistogel_mCherry_Median','Histogel_mCherry_Median'});

disp('=== Median Percentage_energy for hver Band og mCherry Group ===');
disp(MedianTbl);

% 6) Tegn en bar‐graf for at visualisere forskellene
figure('Position',[300 200 450 300]);
bar([median_no, median_hist], 0.6);
set(gca, 'XTick', 1:nBands, 'XTickLabel', bands, 'FontSize',12);
legend({'NoHistogel\_mCherry','Histogel\_mCherry'}, 'Location','best');
ylabel('Median Percentage\_energy (%)');
title('Median Percentage\_energy per Band (mCherry Groups)');

% 7) Tilføj tekstværdier ovenpå søjlerne
yMax = max([median_no; median_hist]);
offset = 0.02 * yMax;
for k = 1:nBands
    text(k-0.15, median_no(k) + offset, sprintf('%.2f', median_no(k)), ...
         'FontSize',10, 'HorizontalAlignment','center');
    text(k+0.15, median_hist(k) + offset, sprintf('%.2f', median_hist(k)), ...
         'FontSize',10, 'HorizontalAlignment','center');
end


%% Parvise t‐tests i hvert bånd for 
T_hist_mCherry = T12(T12.Group == 'Histogel_mCherry', :);  % Undergruppe: kun Histogel_mCherry
bands  = categories(T_hist_mCherry.Band);
conds  = categories(T_hist_mCherry.Condition);    % {'0','25','100'}
pairIdx = [1 2; 1 3; 2 3];                        % 0 vs 25, 0 vs 100, 25 vs 100

fprintf('=== Parvise t‐tests for HISTOGEL\_mCherry kun ===\n');
for k = 1:numel(bands)
    b = bands{k};
    fprintf('\n-- Bånd: %s --\n', b);
    
    for m = 1:size(pairIdx,1)
        i1 = pairIdx(m,1);
        i2 = pairIdx(m,2);
        dose1 = conds{i1};  
        dose2 = conds{i2};  
        
        % Udtræk data for Histogel_mCherry i dette bånd og de to doser
        data1 = T_hist_mCherry.Percentage_energy((T_hist_mCherry.Band == b) & (T_hist_mCherry.Condition == dose1));
        data2 = T_hist_mCherry.Percentage_energy((T_hist_mCherry.Band == b) & (T_hist_mCherry.Condition == dose2));
        
        if isempty(data1) || isempty(data2)
            fprintf('  %s µM vs. %s µM: Ingen data til rådighed.\n', dose1, dose2);
        else
            [h, pval, ~, stats] = ttest2(data1, data2, 'Vartype','unequal');
            sign = '(n.s.)';
            if pval < 0.05, sign = '(*): p<0.05'; end
            fprintf('  %s µM vs. %s µM: t(%d)=%.2f, p=%.3e %s\n', ...
                    dose1, dose2, stats.df, stats.tstat, pval, sign);
        end
    end
end

%% Parvise t‐tests i hvert bånd 
T_no_mCherry = T12(T12.Group == 'NoHistogel_mCherry', :);  % Undergruppe: kun NoHistogel_mCherry
bands  = categories(T_no_mCherry.Band);
conds  = categories(T_no_mCherry.Condition);       % {'0','25','100'}
pairIdx = [1 2; 1 3; 2 3];                         % 0 vs 25, 0 vs 100, 25 vs 100

fprintf('\n=== Parvise t‐tests for NOHISTOGEL\_mCherry kun ===\n');
for k = 1:numel(bands)
    b = bands{k};
    fprintf('\n-- Bånd: %s --\n', b);
    
    for m = 1:size(pairIdx,1)
        i1 = pairIdx(m,1);
        i2 = pairIdx(m,2);
        dose1 = conds{i1};
        dose2 = conds{i2};
        
        % Udtræk NoHistogel_mCherry-data for dette bånd og de to doser
        data1 = T_no_mCherry.Percentage_energy((T_no_mCherry.Band == b) & (T_no_mCherry.Condition == dose1));
        data2 = T_no_mCherry.Percentage_energy((T_no_mCherry.Band == b) & (T_no_mCherry.Condition == dose2));
        
        if isempty(data1) || isempty(data2)
            fprintf('  %s µM vs. %s µM: Ingen data til rådighed.\n', dose1, dose2);
        else
            [h, pval, ~, stats] = ttest2(data1, data2, 'Vartype','unequal');
            sign = '(n.s.)';
            if pval < 0.05, sign = '(*): p<0.05'; end
            fprintf('  %s µM vs. %s µM: t(%d)=%.2f, p=%.3e %s\n', ...
                    dose1, dose2, stats.df, stats.tstat, pval, sign);
        end
    end
end


%% Hypotese 3: Antal subplots for alle mCherry-data
T_plot = T12;  % Brug alle mCherry-data
bands = categories(T_plot.Band);       % {'Low','Mid','High'}
conds = categories(T_plot.Condition);  % {'0','25','100'}

% 1) Find fælles y-akse-grænser
allValues = T_plot.Percentage_energy;  
ymin = floor(min(allValues)) - 1;
ymax = ceil( max(allValues)) + 1;

% 2) Opret figur med tre underplots (én for hvert bånd)
figure('Position',[100 100 1000 300]);

for b = 1:numel(bands)
    bandName = bands{b};
    subplot(1,3,b);
    hold on;
    
    % For hvert dosis‐niveau (0,25,100), hent data
    groupData = cell(1,numel(conds));
    for i = 1:numel(conds)
        c = conds{i};
        sel = (T_plot.Band == bandName) & (T_plot.Condition == c);
        groupData{i} = T_plot.Percentage_energy(sel);
    end
    
    % Kombiner i én vektor + gruppetag for boxplot
    allData = vertcat(groupData{:});
    groupLabels = [ ...
        repmat(string(conds(1)), numel(groupData{1}), 1); ...
        repmat(string(conds(2)), numel(groupData{2}), 1); ...
        repmat(string(conds(3)), numel(groupData{3}), 1)  ...
    ];
    
    boxplot(allData, groupLabels, 'Whisker',1.5, 'LabelVerbosity','minor');
    
    title(['Bånd: ' bandName], 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Cisplatin (µM)', 'FontSize', 11);
    if b == 1
        ylabel('Percentage\_energy (%)', 'FontSize', 11);
    end
     % Sæt y-akse til 0–60 % i stedet for [ymin ymax]
    ylim([0 80]);
    set(gca, 'FontSize', 10);
    
    hold off;
end

sgtitle('Boxplots: Percentage\_energy pr. bånd og dosis (mCherry)', 'FontSize', 14, 'FontWeight', 'bold');

%% BOKSPLOT kun med gruppe x band
figure('Position',[150 150 1000 300]);
bands = categories(T12.Band);

for b = 1:3
    subplot(1,3,b);
    sel   = (T12.Band == bands{b});
    data  = T12.Percentage_energy(sel);
    group = T12.Group(sel);  % fx kategorier ['NoHistogel' 'Histogel']

    % Her angiver vi custom labels:
    customLabels = {'Histogel_mCherry','NoHistogel_mCherry'};
    boxplot(data, group, 'Whisker',1.5, 'Labels', customLabels);

    title(['Band: ' bands{b}], 'FontSize', 12);
    ylabel('Percentage\_energy (%)','FontSize', 11);
    set(gca, 'FontSize', 10);
    ylim([0 70]);
end


%% BESKRIVENDE STATISTIK FOR “Percentage_energy” (Group × Condition × Band)

% 1) Filtrér kun de to mCherry-grupper
T12 = T(ismember(T.Group, {'NoHistogel_mCherry','Histogel_mCherry'}), :);
T12.Group     = removecats(T12.Group);
T12.Condition = removecats(T12.Condition);
T12.Band      = removecats(T12.Band);

% 2) Definér unikke niveauer
groups = categories(T12.Group);       % {'Histogel_mCherry','NoHistogel_mCherry'}
conds  = categories(T12.Condition);   % {'0','25','100'}
bands  = categories(T12.Band);        % {'Low','Mid','High'}

% 3) Forbered en tom tabel til resultatet
varTypes = [ repmat({'string'}, 1, 3), repmat({'double'}, 1, 7) ];
SummaryTbl = table('Size',[0, 10], ...
                   'VariableTypes', varTypes, ...
                   'VariableNames', {'Group','Condition','Band', ...
                                     'Count','Mean','Std','Median','Q1','Q3','IQR'});

% 4) Loop igennem hver kombination af Group, Condition og Band
for g = 1:numel(groups)
    grp = groups{g};
    for i = 1:numel(conds)
        cond_i = conds{i};
        for j = 1:numel(bands)
            b = bands{j};
            
            % Udvælg data for Group-Condition-Band-kombination (mCherry)
            sel = (T12.Group == grp) & ...
                  (T12.Condition == cond_i) & ...
                  (T12.Band == b);
            data = T12.Percentage_energy(sel);
            
            if isempty(data)
                cnt    = 0;
                meanV  = NaN;
                stdV   = NaN;
                medV   = NaN;
                q1V    = NaN;
                q3V    = NaN;
                iqrV   = NaN;
            else
                cnt    = numel(data);
                meanV  = mean(data);
                stdV   = std(data);
                medV   = median(data);
                q1V    = prctile(data, 25);
                q3V    = prctile(data, 75);
                iqrV   = q3V - q1V;
            end
            
            newRow = {grp, cond_i, b, cnt, meanV, stdV, medV, q1V, q3V, iqrV};
            SummaryTbl = [SummaryTbl; newRow];
        end
    end
end

% 5) Vis den samlede tabel
disp('=== Beskrivende statistik (Percentage_energy) for hver mCherry Group×Condition×Band ===');
disp(SummaryTbl);


%% Én-vejs ANOVA (Condition-effekt pr. bånd, mCherry)
for k = 1:numel(bands)
    b = bands{k};
    
    sel = (T12.Band == b);
    Y   = T12.Percentage_energy(sel);
    G   = T12.Condition(sel);  % kategorisk
    
    [p, tbl, stats] = anova1(Y, G, 'off');
    fprintf('Band = %-4s: ANOVA p = %.8e\n', b, p);
end


%% POST-HOC: Parvise sammenligninger pr. bånd (Tukey-Kramer, mCherry)
for k = 1:numel(bands)
    b = bands{k};
    sel = (T12.Band == b);
    Y   = T12.Percentage_energy(sel);
    G   = T12.Condition(sel);
    
    [p, tbl, stats] = anova1(Y, G, 'off');
    
    fprintf('\n=== Band = %s (ANOVA p = %.2e) ===\n', b, p);
    fprintf('Parvise sammenligninger (Tukey-Kramer):\n');
    
    figure('Name',['Post-hoc: ' b],'NumberTitle','off');
    multcompare(stats, 'CType','tukey-kramer');
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%% HYPOTESE 5 %%%%%%%
       % Mcherry vs band og no m cherry vs band

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
filename   = 'Statistik_data_frekvens.xlsx';
sheets     = {'Ark1','Ark2','Ark3','Ark4'};
groupNames = { ...
   'NoHistogel', ...
   'Histogel', ...
   'NoHistogel_mCherry', ...
   'Histogel_mCherry' ...
};

T = table();
for i = 1:4
    Ti = readtable(filename, 'Sheet', sheets{i});
    Ti.Group = repmat(string(groupNames{i}), height(Ti), 1);
    T = [T; Ti];
end

% Konverter ’Condition’, ’Band’ og ’Group’ til kategoriske variabler
T.Condition = categorical(T.Condition, [0 25 100], {'0','25','100'});
T.Band      = categorical(T.Band,      {'Low','Mid','High'});
T.Group     = categorical(T.Group);

% 1) Lav en ny kolonne ’mCherry’ som ’Ja’ hvis Group‐navnet indeholder '_mCherry', ellers ’Nej’
isMCherry     = contains(string(T.Group), "_mCherry");
T.mCherry     = categorical( repmat("Nej", height(T),1) );
T.mCherry(isMCherry) = categorical("Ja");


%% Two‐way ANOVA på Absolut_energy med mCherry × Band
fprintf('\n=== Two-way ANOVA på Absolut_energy (mCherry × Band) ===\n');
[p, tbl, stats] = anovan( ...
    T.Absolut_energy, ...            % respons‐variabel
    {T.mCherry, T.Band}, ...         % faktorer: mCherry, Band
    'model',    'interaction', ...   % hoved- + interaktionseffekter
    'varnames', {'mCherry','Band'}, ...
    'display',  'on' );              


%% deskriptiv Mean
% --- 0) Sørg for, at mCherry er kategorisk, og fjern NaN i Absolut_energy ---
if ~iscategorical(T.mCherry)
    T.mCherry = categorical(T.mCherry, {'Nej','Ja'});
end
idxNaN = isnan(T.Absolut_energy);
T(idxNaN,:) = [];

% --- 1) Hent de to levels (Nej, Ja) i mCherry-kolonnen ---
levels = categories(T.mCherry);   % -> {'Nej','Ja'}
nLevels = numel(levels);

% For at gemme resultater:
groupCount = nan(nLevels,1);
groupMean  = nan(nLevels,1);
groupStd   = nan(nLevels,1);
groupSEM   = nan(nLevels,1);
groupCIlo  = nan(nLevels,1);
groupCIhi  = nan(nLevels,1);

alpha = 0.05;  % for 95% CI

for i = 1:nLevels
    lvl = levels{i};
    sel = (T.mCherry == lvl);
    vals = T.Absolut_energy(sel);
    
    groupCount(i) = numel(vals);
    mu = mean(vals);
    sigma = std(vals);
    sem = sigma / sqrt(groupCount(i));
    df = groupCount(i) - 1;
    tcrit = tinv(1 - alpha/2, df);  % kritisk t-værdi
    
    ci_half = tcrit * sem;
    groupMean(i) = mu;
    groupStd(i)  = sigma;
    groupSEM(i)  = sem;
    groupCIlo(i) = mu - ci_half;
    groupCIhi(i) = mu + ci_half;
end

T_summary = table( ...
    levels, ...
    groupCount, ...
    groupMean, ...
    groupStd, ...
    groupSEM, ...
    groupCIlo, ...
    groupCIhi, ...
    'VariableNames', ...
    {'mCherry','Count','Mean_AbsEnergy','StdDev','SEM','CI_Lower','CI_Upper'} ...
);
disp('=== Deskriptiv statistik for Absolut_energy pr. mCherry-gruppe ===');
disp(T_summary);

