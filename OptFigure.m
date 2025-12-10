addpath(genpath('C:\Users\Evan\Documents\GitHub\eelib'))

%% Load Data
full_tbl = load_opt_data("Data/Opt");
full_tbl = full_tbl(~strcmp(full_tbl.Rat, 'MORT13'),:);
full_tbl.sessionType = strings([height(full_tbl),1]);
full_tbl.sessionType(:) = "off";
full_tbl.sessionType(contains(full_tbl.Protocol, 'all')) = "on";
full_tbl.sessionType(~strcmp(full_tbl.Opt,"")) = "opt";
[~,I] = sort(full_tbl.sessionType,'ascend');
full_tbl = full_tbl(I,:);
full_tbl.omission = full_tbl.RT>3;
full_tbl.DT = full_tbl.initTime - full_tbl.cueTime;
full_tbl.Sex = extractBefore(full_tbl.Rat,2);

%% Supplemental Target Figure
coords = readtable('Data/FORT_MORT_Targets.csv');
coords = coords(~isnan(coords.ML),:);
coords = coords(strcmp(coords.Experiments, 'Opt'),:);

figure
addpath(genpath('F:\AtlasPlotter'))
slice = Slice(mean(coords.AP), 'c');
slice.cut([0, -20; 0, 5], 'L')
slice.plot()
col=[0.949,0.631,0.008];
hold on
scatter(coords.ML(~strcmp(coords.Rat,'MORT13')),coords.DV(~strcmp(coords.Rat,'MORT13')),60,"filled",'MarkerFaceColor',col,'MarkerEdgeColor',[0.3,0.3,0.3])
scatter(coords.ML(strcmp(coords.Rat,'MORT13')),coords.DV(strcmp(coords.Rat,'MORT13')),60,"filled",'MarkerFaceColor','r','MarkerEdgeColor',[0.3,0.3,0.3])
axis equal
axis off

%% Create Session table
sessions = unique(full_tbl.Session);
ses_tbl = table('Size', [length(sessions), 11], 'VariableTypes', ["string","string","double","double","double","double","double","double","double","double","double"], 'VariableNames',["Rat","sessionType", "performance", "omission", "frontChoice","performance0", "omission0", "frontChoice0","rulesCompleted", "tcstd", "tcmu"]);
for i=1:length(sessions)
    sub_tbl = full_tbl(full_tbl.Session==sessions(i),:);
    ses_tbl.sessionType(i) = sub_tbl.sessionType(1);
    ses_tbl.Rat(i) = sub_tbl.Rat(1);
    ses_tbl.protocolCode(i) = extractBefore(sub_tbl.Protocol(1),3);
    rule = strcmp(sub_tbl.rule, 'L');
    ses_tbl.rulesCompleted(i) = sum(abs(diff(rule)));
    ses_tbl.performance(i) = mean(sub_tbl.performance(sub_tbl.RT < 3));
    ses_tbl.omission(i) = mean(sub_tbl.RT > 3);
    ses_tbl.frontChoice(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3));
    ses_tbl.rewards(i) = sum(sub_tbl.performance);
end

%% Create Rat Table
rats = unique(full_tbl.Rat);
rat_tbl = table('Size', [length(rats), 19], 'VariableTypes', ["string","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double"], ...
    'VariableNames',["Rat", "RT_off", "DT_off", "performance_off", "omission_off", "frontChoice_off", "rulesCompleted_off", ...
                     "RT_on", "DT_on", "performance_on", "omission_on", "frontChoice_on", "rulesCompleted_on", ...
                     "RT_opt", "DT_opt", "performance_opt", "omission_opt", "frontChoice_opt", "rulesCompleted_opt"]);
for i=1:length(rats)
    sub_tbl = full_tbl(strcmp(full_tbl.Rat,rats(i))&strcmp(full_tbl.sessionType, 'off'),:);
    rat_tbl.Rat(i) = rats(i);
    rat_tbl.RT_off(i) = median(sub_tbl.RT(sub_tbl.RT < 3));
    rat_tbl.DT_off(i) = median(sub_tbl.DT);
    rat_tbl.performance_off(i) = mean(sub_tbl.performance(sub_tbl.RT < 3));
    rat_tbl.omission_off(i) = mean(sub_tbl.RT > 3);
    rat_tbl.frontChoice_off(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3));
    rat_tbl.rulesCompleted_off(i) = mean(ses_tbl.rulesCompleted(strcmp(ses_tbl.Rat, rats(i))&strcmp(ses_tbl.sessionType, 'off')));
    sub_tbl = full_tbl(strcmp(full_tbl.Rat,rats(i))&strcmp(full_tbl.sessionType, 'on'),:);
    rat_tbl.RT_on(i) = median(sub_tbl.RT(sub_tbl.RT < 3));
    rat_tbl.DT_on(i) = median(sub_tbl.DT);
    rat_tbl.performance_on(i) = mean(sub_tbl.performance(sub_tbl.RT < 3));
    rat_tbl.omission_on(i) = mean(sub_tbl.RT > 3);
    rat_tbl.frontChoice_on(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3));
    rat_tbl.rulesCompleted_on(i) = mean(ses_tbl.rulesCompleted(strcmp(ses_tbl.Rat, rats(i))&strcmp(ses_tbl.sessionType, 'on')));
    sub_tbl = full_tbl(strcmp(full_tbl.Rat,rats(i))&strcmp(full_tbl.sessionType, 'opt'),:);
    rat_tbl.RT_opt(i) = median(sub_tbl.RT(sub_tbl.RT < 3));
    rat_tbl.DT_opt(i) = median(sub_tbl.DT);
    rat_tbl.performance_opt(i) = mean(sub_tbl.performance(sub_tbl.RT < 3));
    rat_tbl.omission_opt(i) = mean(sub_tbl.RT > 3);
    rat_tbl.frontChoice_opt(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3));
    rat_tbl.rulesCompleted_opt(i) = mean(ses_tbl.rulesCompleted(strcmp(ses_tbl.Rat, rats(i))&strcmp(ses_tbl.sessionType, 'opt')));
end

%% Run GLMs
full_tbl.rule_type = strcmp(full_tbl.rule, 'L');
full_tbl.protocolCode = extractBefore(full_tbl.Protocol, 3);
rtglm = fitglme(full_tbl(full_tbl.RT<3,:), 'RT~1+sessionType+rule_type+(1|protocolCode)+(1|Session)+(1|Rat)','distribution','gamma','link','identity');
perfglm = fitglme(full_tbl(full_tbl.RT<3,:), 'performance~1+sessionType+rule_type+(1|protocolCode)+(1|Session)+(1|Rat)','distribution','binomial');

%% Create figure
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% E
subplot_tight(4,2,[5,7])
groups = {full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.sessionType,'off')), ...
          full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.sessionType,'on')), ...
          full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.sessionType,'opt'))};
al_goodplot2(groups', 'pos', [1,2,3], 'type', {'bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0.18,0.55,0.34],'useMedian', true)
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.RT_off, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.RT_on, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)),rat_tbl.RT_opt, 'MarkerFaceColor', [0.18,0.55,0.34], 'MarkerEdgeColor', [0.18,0.55,0.34])

if coefTest(rtglm, [0,1,0,0]) < 0.05
    plot([1,2], [0.75,0.75], 'k', 'LineWidth', 2)
    text(1.5, 0.76, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(rtglm, [0,0,1,0]) < 0.05
    plot([1,3], [0.9,0.9], 'k', 'LineWidth', 2)
    text(2, 0.91, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:3)
ylim([0.2,1])
xticklabels(["OFF", "ON", "OPT"])
ylabel('RT (s)')
xlim([0.5,3.5])
set(gca,'linewidth',2)
set(gca, 'fontsize', 18)

%% F
subplot_tight(4,2,[6,8])
groups = {ses_tbl.performance(strcmp(ses_tbl.sessionType,'off')), ...
          ses_tbl.performance(strcmp(ses_tbl.sessionType,'on')), ...
          ses_tbl.performance(strcmp(ses_tbl.sessionType,'opt'))};
al_goodplot2(groups', 'pos', [1,2,3], 'type', {'bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0.18,0.55,0.34],'useMedian', true)
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.performance_off, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.performance_on, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)),rat_tbl.performance_opt, 'MarkerFaceColor', [0.18,0.55,0.34], 'MarkerEdgeColor', [0.18,0.55,0.34])

xticks(1:3)
ylim([0.55,0.77])
xticklabels(["OFF", "ON", "OPT"])
ylabel('Accuracy')
xlim([0.5,3.5])
set(gca,'linewidth',2)
set(gca, 'fontsize', 18)

%% Supplement RLDDM
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% A
subplot_tight(2,2,1)
groups = {full_tbl.DT(strcmp(full_tbl.sessionType,'off')), ...
          full_tbl.DT(strcmp(full_tbl.sessionType,'on')), ...
          full_tbl.DT(strcmp(full_tbl.sessionType,'opt'))};
al_goodplot2(groups', 'pos', [1,2,3], 'type', {'bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0.18,0.55,0.34],'useMedian', true)
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.DT_off, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.DT_on, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)),rat_tbl.DT_opt, 'MarkerFaceColor', [0.18,0.55,0.34], 'MarkerEdgeColor', [0.18,0.55,0.34])

if coefTest(dtglm, [0,0,1]) < 0.05 / 6
    plot([1,3], [12,12], 'k', 'LineWidth', 2)
    text(2, 13, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(dtglm, [0,1,0]) < 0.05 / 6
    plot([1,2], [8,8], 'k', 'LineWidth', 2)
    text(1.5, 9, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:3)
set(gca, 'YScale', 'log')
ylim([0.4,16])
yticks([0.5,1,2,4,8,16])
xticklabels(["OFF", "ON", "OPT"])
ylabel('DT (s)')
xlim([0.5,3.5])
set(gca,'linewidth',2)
set(gca, 'fontsize', 18)
set(gca, 'YMinorTick', 'off')

%% B
subplot_tight(2,2,2)
groups = {ses_tbl.omission(strcmp(ses_tbl.sessionType,'off'))*100, ...
          ses_tbl.omission(strcmp(ses_tbl.sessionType,'on'))*100, ...
          ses_tbl.omission(strcmp(ses_tbl.sessionType,'opt'))*100};
al_goodplot2(groups', 'pos', [1,2,3], 'type', {'bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0.18,0.55,0.34],'useMedian', false)
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.omission_off*100, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.omission_on*100, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)),rat_tbl.omission_opt*100, 'MarkerFaceColor', [0.18,0.55,0.34], 'MarkerEdgeColor', [0.18,0.55,0.34])

xticks(1:3)
xticklabels(["OFF", "ON", "OPT"])
ylabel('% Omissions')
xlim([0.5,3.5])
ylim([0,10])
set(gca,'linewidth',2)
set(gca, 'fontsize', 18)

%% C
subplot_tight(2,2,3)
groups = {ses_tbl.frontChoice(strcmp(ses_tbl.sessionType,'off'))*100, ...
          ses_tbl.frontChoice(strcmp(ses_tbl.sessionType,'on'))*100, ...
          ses_tbl.frontChoice(strcmp(ses_tbl.sessionType,'opt'))*100};
al_goodplot2(groups', 'pos', [1,2,3], 'type', {'bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0.18,0.55,0.34],'useMedian', false)
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.frontChoice_off*100, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.frontChoice_on*100, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)),rat_tbl.frontChoice_opt*100, 'MarkerFaceColor', [0.18,0.55,0.34], 'MarkerEdgeColor', [0.18,0.55,0.34])

if coefTest(fcglm, [1,0,0]) < 0.05 / 5
    text(1, 64, '#', 'HorizontalAlignment','center','fontsize', 20)
end

ylim([35,65])
xticks(1:3)
xticklabels(["OFF", "ON", "OPT"])
ylabel('% Left Choice')
xlim([0.5,3.5])
set(gca,'linewidth',2)
set(gca, 'fontsize', 18)