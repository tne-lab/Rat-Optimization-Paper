addpath(genpath('C:\Users\Evan\Documents\GitHub\eelib'))

%% Load Data
full_tbl = load_opt_data("Data/BB");
full_tbl.sessionType = strings([height(full_tbl),1]);
full_tbl.sessionType(:) = "pr";
full_tbl.sessionType(contains(full_tbl.Protocol, 'off')) = "off";
full_tbl.sessionType(contains(full_tbl.Protocol, 'all')) = "on";
[~,I] = sort(full_tbl.sessionType,'ascend');
full_tbl = full_tbl(I,:);
full_tbl.omission = full_tbl.RT>3;
full_tbl.DT = full_tbl.initTime - full_tbl.cueTime;
full_tbl.Sex = extractBefore(full_tbl.Rat,2);
full_tbl = full_tbl(~strcmp(full_tbl.Rat, 'MORT09'),:);

%% Add trial in rule information
trial = 1;
prevRule = "";
prevSes = "";
progress = 0;
prevStim = 0;
block_trial = 0;
for i=1:height(full_tbl)
    if ~strcmp(prevSes,full_tbl.Session(i))
        prevStim = 0;
        block_trial = 0;
    end
    if ~strcmp(prevRule, full_tbl.rule(i)) || ~strcmp(prevSes,full_tbl.Session(i))
        progress = 0;
        trial = 1;
        prevRule = full_tbl.rule(i);
        prevSes = full_tbl.Session(i);
    end
    full_tbl.tinb(i) = block_trial;
    if full_tbl.stim(i) ~= prevStim
        block_trial = 0;
        prevStim = full_tbl.stim(i);
    elseif strcmp(full_tbl.sessionType(i), 'pr')
        block_trial = block_trial + 1;
    end
    full_tbl.tinrule(i) = trial;
    trial = trial + 1;
    if full_tbl.performance(i) == 1
        progress = progress + 1;
        full_tbl.progress(i) = 1;
    else
        full_tbl.progress(i) = -progress;
        progress = 0;
    end
end

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
    if strcmp(ses_tbl.sessionType(i), 'pr')
        ses_tbl.performance(i) = mean(sub_tbl.performance(sub_tbl.RT < 3&sub_tbl.stim==0));
        ses_tbl.omission(i) = mean(sub_tbl.RT(sub_tbl.stim==0) > 3);
        ses_tbl.frontChoice(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3 &sub_tbl.stim==0));
        ses_tbl.rewards(i) = sum(sub_tbl.performance)/nnz(sub_tbl.stim==0)*height(sub_tbl);
        ses_tbl.progress(i) = mean(sub_tbl.progress(sub_tbl.stim==0));
        ses_tbl.performance0(i) = mean(sub_tbl.performance(sub_tbl.RT < 3&sub_tbl.stim==1));
        ses_tbl.omission0(i) = mean(sub_tbl.RT(sub_tbl.stim==1) > 3);
        ses_tbl.frontChoice0(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3 &sub_tbl.stim==1));
        ses_tbl.rewards0(i) = sum(sub_tbl.performance)/nnz(sub_tbl.stim==1)*height(sub_tbl);
        ses_tbl.progress0(i) = mean(sub_tbl.progress(sub_tbl.stim==1));
    else
        ses_tbl.performance(i) = mean(sub_tbl.performance(sub_tbl.RT < 3));
        ses_tbl.omission(i) = mean(sub_tbl.RT > 3);
        ses_tbl.frontChoice(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3));
        ses_tbl.rewards(i) = sum(sub_tbl.performance);
        ses_tbl.progress(i) = mean(sub_tbl.progress);
    end
end

progress_tbl = table('Size', [length(sessions)+nnz(ses_tbl.progress0), 1], 'VariableTypes', ["string"], 'VariableNames', ["Rat"]);
j = 1;
for i=1:height(ses_tbl)
    progress_tbl.Rat(j) = ses_tbl.Rat(i);
    progress_tbl.sessionType(j) = ses_tbl.sessionType(i);
    progress_tbl.protocolCode(j) = ses_tbl.protocolCode(i);
    progress_tbl.stim(j) = 0;
    progress_tbl.progress(j) = ses_tbl.progress(i);
    if strcmp(ses_tbl.sessionType(i), 'pr')
        progress_tbl.Rat(j+1) = ses_tbl.Rat(i);
        progress_tbl.protocolCode(j+1) = ses_tbl.protocolCode(i);
        progress_tbl.sessionType(j+1) = ses_tbl.sessionType(i);
        progress_tbl.progress(j+1) = ses_tbl.progress0(i);
        progress_tbl.stim(j+1) = 1;
        j = j + 1;
    end
    j = j + 1;
end

%% Create Rat Table
rats = unique(full_tbl.Rat);
rat_tbl = table('Size', [length(rats), 24], 'VariableTypes', ["string","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double"], ...
    'VariableNames',["Rat", "RT_off", "DT_off", "performance_off", "omission_off", "frontChoice_off", "rulesCompleted_off", ...
                     "RT_on", "DT_on", "performance_on", "omission_on", "frontChoice_on", "rulesCompleted_on", ...
                     "RT_pr", "DT_pr", "performance_pr", "omission_pr", "frontChoice_pr", "rulesCompleted_pr", ...
                     "RT_pr0", "DT_pr0", "performance_pr0", "omission_pr0", "frontChoice_pr0"]);
for i=1:length(rats)
    sub_tbl = full_tbl(strcmp(full_tbl.Rat,rats(i))&strcmp(full_tbl.sessionType, 'off'),:);
    rat_tbl.Rat(i) = rats(i);
    rat_tbl.RT_off(i) = median(sub_tbl.RT(sub_tbl.RT < 3));
    rat_tbl.DT_off(i) = median(sub_tbl.DT);
    rat_tbl.performance_off(i) = mean(sub_tbl.performance(sub_tbl.RT < 3));
    rat_tbl.omission_off(i) = mean(sub_tbl.RT > 3);
    rat_tbl.frontChoice_off(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3));
    rat_tbl.rulesCompleted_off(i) = mean(ses_tbl.rulesCompleted(strcmp(ses_tbl.Rat, rats(i))&strcmp(ses_tbl.sessionType, 'off')));
    rat_tbl.progress_off(i) = mean(progress_tbl.progress(strcmp(progress_tbl.Rat, rats(i))&strcmp(progress_tbl.sessionType, 'off')));
    sub_tbl = full_tbl(strcmp(full_tbl.Rat,rats(i))&strcmp(full_tbl.sessionType, 'on'),:);
    rat_tbl.RT_on(i) = median(sub_tbl.RT(sub_tbl.RT < 3));
    rat_tbl.DT_on(i) = median(sub_tbl.DT);
    rat_tbl.performance_on(i) = mean(sub_tbl.performance(sub_tbl.RT < 3));
    rat_tbl.omission_on(i) = mean(sub_tbl.RT > 3);
    rat_tbl.frontChoice_on(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3));
    rat_tbl.rulesCompleted_on(i) = mean(ses_tbl.rulesCompleted(strcmp(ses_tbl.Rat, rats(i))&strcmp(ses_tbl.sessionType, 'on')));
    rat_tbl.progress_on(i) = mean(progress_tbl.progress(strcmp(progress_tbl.Rat, rats(i))&strcmp(progress_tbl.sessionType, 'on')));
    sub_tbl = full_tbl(strcmp(full_tbl.Rat,rats(i))&strcmp(full_tbl.sessionType, 'pr')&full_tbl.stim==0,:);
    rat_tbl.RT_pr(i) = median(sub_tbl.RT(sub_tbl.RT < 3));
    rat_tbl.DT_pr(i) = median(sub_tbl.DT);
    rat_tbl.performance_pr(i) = mean(sub_tbl.performance(sub_tbl.RT < 3));
    rat_tbl.omission_pr(i) = mean(sub_tbl.RT > 3);
    rat_tbl.frontChoice_pr(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3));
    rat_tbl.progress_pr(i) = mean(progress_tbl.progress(strcmp(progress_tbl.Rat, rats(i))&strcmp(progress_tbl.sessionType, 'pr')&progress_tbl.stim==0));
    sub_tbl = full_tbl(strcmp(full_tbl.Rat,rats(i))&strcmp(full_tbl.sessionType, 'pr')&full_tbl.stim==1,:);
    rat_tbl.RT_pr0(i) = median(sub_tbl.RT(sub_tbl.RT < 3));
    rat_tbl.DT_pr0(i) = median(sub_tbl.DT);
    rat_tbl.performance_pr0(i) = mean(sub_tbl.performance(sub_tbl.RT < 3));
    rat_tbl.omission_pr0(i) = mean(sub_tbl.RT > 3);
    rat_tbl.frontChoice_pr0(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3));
    rat_tbl.progress_pr0(i) = mean(progress_tbl.progress(strcmp(progress_tbl.Rat, rats(i))&strcmp(progress_tbl.sessionType, 'pr')&progress_tbl.stim==1));
end

%% Run GLMs
full_tbl.rule_type = strcmp(full_tbl.rule, 'L');
full_tbl.protocolCode = extractBefore(full_tbl.Protocol, 3);
rtglm = fitglme(full_tbl(full_tbl.RT<3,:), 'RT~1+stim+sessionType+rule_type+(1|Session)+(1|Rat)+(1|protocolCode)','distribution','gamma','link','identity');
perfglm = fitglme(full_tbl(full_tbl.RT<3,:), 'performance~1+stim+sessionType+rule_type+(1|Session)+(1|Rat)+(1|protocolCode)','distribution','binomial');
omglm = fitglme(full_tbl(:,:), 'omission~1+sessionType+stim+rule_type+(1|Session)+(1|Rat)+(1|protocolCode)','distribution','binomial');
dtglm = fitglme(full_tbl, 'DT~1+sessionType+stim+rule_type+(1|Session)+(1|Rat)+(1|protocolCode)','distribution','gamma','link','log');
fcglm = fitglme(full_tbl(full_tbl.RT<3,:), 'frontChoice~1+sessionType+stim+rule_type+(1|Session)+(1|Rat)+(1|protocolCode)','distribution','binomial');
progglm = fitglme(sortrows(progress_tbl, 2), 'progress~1+sessionType+stim+(1|Rat)+(1|protocolCode)');


%%
tt = full_tbl(full_tbl.RT<3,:);
test = tt.RT - predict(rtglm, tt) + tt.stim * double(rtglm.Coefficients(2,2));
g1 = test(strcmp(tt.sessionType,'pr')&tt.stim==0);
n1 = length(g1);
g2 = test(strcmp(tt.sessionType,'pr')&tt.stim==1);
n2 = length(g2);
(mean(g2) - mean(g1)) / sqrt(((n1-1)*var(g1)+(n2-1)*var(g2))/(n1+n2-2))
%%
tt = full_tbl(full_tbl.RT<3,:);
test = tt.RT;
g1 = test(strcmp(tt.sessionType,'off'));
n1 = length(g1);
g2 = test(strcmp(tt.sessionType,'on'));
n2 = length(g2);
(mean(g2) - mean(g1)) / sqrt(((n1-1)*var(g1)+(n2-1)*var(g2))/(n1+n2-2))
%% Setup figure
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% A
subplot_tight(8,9,[1,2,3,10,11,12,19,20,21,28,29,30])
coords = readtable('Data/FORT_MORT_Targets.csv');
coords = coords(~strcmp(coords.Rat, 'MORT09'), :);
coords = coords(~isnan(coords.ML),:);
coords = coords(~strcmp(coords.Experiments, 'Opt'));
addpath(genpath('F:\AtlasPlotter'))
slice = Slice(mean(coords.AP), 'c');
slice.cut([0, -20; 0, 5], 'L')
slice.plot()
col=[0.949,0.631,0.008];
hold on
scatter(coords.ML,coords.DV,60,"filled",'MarkerFaceColor',col,'MarkerEdgeColor',[0.3,0.3,0.3])
axis equal
axis off

%% B
subplot_tight(8,9,4:6)
stim_runs = full_tbl.stim(strcmp(full_tbl.Session, "123"));
run_length = 0;
condition = 0;
for i=1:length(stim_runs)
    if stim_runs(i) ~= condition
        if condition == 1
            color = 'r';
        else
            color = 'b';
        end
        patch('XData', [i - run_length+0.3, i - run_length + 0.3, i - 0.3, i - 0.3], 'YData', [0, 1, 1, 0], 'FaceColor', color, 'FaceAlpha', 0.5,'EdgeColor','none')
        condition = stim_runs(i);
        run_length = 0;
    end
    run_length = run_length + 1;
end
patch('XData', [i - run_length+0.3, i - run_length + 0.3, i - 0.3, i - 0.3], 'YData', [0, 1, 1, 0], 'FaceColor', 'b', 'FaceAlpha', 0.5,'EdgeColor','none')
axis tight
axis off
ylim([-1,1])

%% D & E
valid = contains(full_tbl.Protocol, 'on')&~contains(full_tbl.Protocol, 'all');
sessions = unique(full_tbl.Session(valid));
valid_tbl = full_tbl(valid,:);
blocks = nan(16*length(sessions),9,2);
inds = [1, 0];
stim = 0;
i = 1;
for session=sessions'
    sub_tbl = valid_tbl(strcmp(valid_tbl.Session, session), :);
    tinb = 1;
    for j=1:height(sub_tbl)
        if stim ~= sub_tbl.stim(j)
            stim = ~stim;
            inds(stim + 1) = inds(stim + 1) + 1;
            tinb = 1;
        end
        valid_tbl.trial_num(i) = j;
        valid_tbl.tinb(i) = tinb;
        i = i + 1;
        blocks(inds(stim+1), tinb, stim+1) = sub_tbl.RT(j);
        tinb = tinb + 1;
    end
    inds = inds + 1;
end

blocks2 = nan(size(blocks));
for i=1:length(blocks)
    last = find(~isnan(blocks(i,:,1)),1,'last');
    blocks2(i,(10 - last):9,1) = blocks(i,1:last,1);
    last = find(~isnan(blocks(i,:,2)),1,'last');
    blocks2(i,(10 - last):9,2) = blocks(i,1:last,2);
end

subplot_tight(8,9,[13,14,15,31,32,33])
on_off_to_on = blocks(:,:,2);
off_off_to_on = blocks2(mod(1:length(blocks2),16) ~= 0,:,1);
hold on
rc1 = ribbon_ci(1:9, on_off_to_on, @nanmean);
rc1.FaceColor = [1,0,0];
rc1.FaceAlpha = 0.3;
rc1.EdgeColor = 'none';
plot(nanmean(on_off_to_on),'linewidth',2, 'Color', [1,0,0])
rc2 = ribbon_ci(-9:-1, off_off_to_on, @nanmean);
rc2.FaceColor = [0,0,1];
rc2.FaceAlpha = 0.3;
rc2.EdgeColor = 'none';
plot(-9:-1, nanmean(off_off_to_on),'linewidth',2, 'Color', [0,0,1])
ylabel('RT (s)')
title('Off to On')
xlabel('Trials Since Change')
yticks(0.5:0.1:0.9)
xlim([-9, 9])
ylim([0.5,0.9])
xticks([-9,0,9])
set(gca,'fontsize',18)
set(gca,'linewidth',2)
set(gca,'fontname','Helvetica')

subplot_tight(8,9,[16,17,18,34,35,36])
off_on_to_off = blocks(mod(1:length(blocks), 16)~=1,:,1);
on_on_to_off = blocks2(:,:,2);
hold on
rc3 = ribbon_ci(-9:-1, on_on_to_off, @nanmean);
rc3.FaceColor = [1,0,0];
rc3.FaceAlpha = 0.3;
rc3.EdgeColor = 'none';
plot(-9:-1, nanmean(on_on_to_off),'linewidth',2,'Color',[1,0,0])
rc4 = ribbon_ci(1:9, off_on_to_off, @nanmean);
rc4.FaceColor = [0,0,1];
rc4.FaceAlpha = 0.3;
rc4.EdgeColor = 'none';
plot(nanmean(off_on_to_off),'linewidth',2, 'Color', [0,0,1])
title('On to Off')
xlabel('Trials Since Change')
xlim([-9, 9])
ylim([0.5,0.9])
yticklabels([])
xticks([-9,0,9])
yticks(0.5:0.1:0.9)
set(gca,'fontsize',18)
set(gca,'linewidth',2)
set(gca,'fontname','Helvetica')

%% F
subplot_tight(8,9,[37,38,39,46,47,48,55,56,57,64,64,66])
groups = {full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.sessionType,'off')), ...
          full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.sessionType,'on')), ...
          full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.sessionType,'pr')&full_tbl.stim==0), ...
          full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.sessionType,'pr')&full_tbl.stim==1)};
al_goodplot2(groups', 'pos', [1,2,3,3], 'type', {'bilateral','bilateral','left','right'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0,0,1;1,0,0],'useMedian', true)
plot([ones(size(rat_tbl.RT_off)), 2*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.RT_off, rat_tbl.RT_on, nan(size(rat_tbl.RT_off))]', 'Color', [0.6,0.6,0.6], 'LineStyle','--')
plot([2.9*ones(size(rat_tbl.RT_off)), 3.1*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.RT_pr, rat_tbl.RT_pr0, nan(size(rat_tbl.RT_off))]', 'Color', [0.6,0.6,0.6], 'LineStyle','--')
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.RT_off, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.RT_on, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)) - 0.1,rat_tbl.RT_pr, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(3*ones(size(rat_tbl.RT_off)) + 0.1,rat_tbl.RT_pr0, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])

if coefTest(rtglm, [0,0,1,0,0]) < 0.05 / 5
    plot([1,2], [0.9,0.9], 'k', 'LineWidth', 2)
    text(1.5, 0.91, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(rtglm, [0,1,0,0,0]) < 0.05 / 5
    plot([2.9,3.1], [0.9,0.9], 'k', 'LineWidth', 2)
    text(3, 0.91, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:3)
ylim([0.2,1])
xticklabels(["OFF", "ON", "BB"])
ylabel('RT (s)')
xlim([0.5,3.5])
set(gca,'linewidth',2)
set(gca, 'fontsize', 18)

%% G
subplot_tight(8,9,[40,41,42,49,50,51,58,59,60,67,68,69])
groups = {ses_tbl.performance(strcmp(ses_tbl.sessionType,'off')), ...
          ses_tbl.performance(strcmp(ses_tbl.sessionType,'on')), ...
          ses_tbl.performance(strcmp(ses_tbl.sessionType,'pr')), ...
          ses_tbl.performance0(strcmp(ses_tbl.sessionType,'pr'))};
al_goodplot2(groups', 'pos', [1,2,3,3], 'type', {'bilateral','bilateral','left','right'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0,0,1;1,0,0])
plot([ones(size(rat_tbl.RT_off)), 2*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.performance_off, rat_tbl.performance_on, nan(size(rat_tbl.RT_off))]', 'Color', [0.6,0.6,0.6], 'LineStyle','--')
plot([2.9*ones(size(rat_tbl.RT_off)), 3.1*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.performance_pr, rat_tbl.performance_pr0, nan(size(rat_tbl.RT_off))]', 'Color', [0.6,0.6,0.6], 'LineStyle','--')
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.performance_off, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.performance_on, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)) - 0.1,rat_tbl.performance_pr, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(3*ones(size(rat_tbl.RT_off)) + 0.1,rat_tbl.performance_pr0, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])

if coefTest(perfglm, [0,0,1,0,0]) < 0.05 / 5
    plot([1,2], [0.75,0.75], 'k', 'LineWidth', 2)
    text(1.5, 0.755, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(perfglm, [0,0,0,1,0]) < 0.05 / 5
    plot([1,2.9], [0.73,0.73], 'k', 'LineWidth', 2)
    text(2.05, 0.735, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(perfglm, [0,1,0,0,0]) < 0.05 / 5
    plot([2.9,3.1], [0.73,0.73], 'k', 'LineWidth', 2)
    text(3, 0.735, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(perfglm, [0,1,-1,1,0]) < 0.05 / 5
    plot([2,3.1], [0.71,0.71], 'k', 'LineWidth', 2)
    text(2.55, 0.715, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:3)
xticklabels(["OFF", "ON", "BB"])
ylabel('Accuracy')
ylim([0.5,0.77])
xlim([0.5,3.5])
set(gca,'linewidth',2)
set(gca, 'fontsize', 18)

%% RLDDM
chains = dir(strcat("Data/BBRLDDM"));
chains = chains(3:end);
traces = readtable(fullfile(chains(1).folder, chains(1).name));
for i=2:length(chains)
    traces = [traces; readtable(fullfile(chains(i).folder, chains(i).name))];
end
names = traces.Properties.VariableNames;

%% Gelman Rubin
L = height(traces) / length(chains);
for i=2:length(names)
    chain_means = zeros(size(chains));
    within_var = zeros(size(chains));
    for j=1:length(chain_means)
        chain_means(j) = mean(traces{1+L*(j-1):L*j,i});
        within_var(j) = std(traces{1+L*(j-1):L*j,i});
    end
    grand_mean = mean(chain_means);
    B = L/(length(chains)-1)*sum((chain_means-grand_mean).^2);
    W = mean(within_var);
    disp(strcat(names(i), ": ", num2str(((L-1)/L*W+B/L)/W)))
end

%% Make traces
a_trace = table2array(traces(:, startsWith(names, 'a_subj')));
v_trace = table2array(traces(:, startsWith(names, 'v_subj')));
n_trace = table2array(traces(:, startsWith(names, 't_subj')));
l_trace = table2array(traces(:, startsWith(names, 'alpha_subj')));
b_trace = table2array(traces(:, startsWith(names, 'zt_subj')));
f_trace = table2array(traces(:, startsWith(names, 'forg_subj')));
s_trace = table2array(traces(:, startsWith(names, 'surp_subj')));
st_trace = table2array(traces(:, startsWith(names, 'st')));
a_stim = table2array(traces(:, startsWith(names, 'a_stim_subj')));
v_stim = table2array(traces(:, startsWith(names, 'v_stim_subj')));
t_stim = table2array(traces(:, startsWith(names, 't_stim_subj')));
b_stim = table2array(traces(:, startsWith(names, 'zt_stim_subj')));
a_ON = table2array(traces(:, startsWith(names, 'a_stimON_subj')));
v_ON = table2array(traces(:, startsWith(names, 'v_stimON_subj')));
t_ON = table2array(traces(:, startsWith(names, 't_stimON_subj')));
b_ON = table2array(traces(:, startsWith(names, 'zt_stimON_subj')));
a_BB = table2array(traces(:, startsWith(names, 'a_stimBB_subj')));
v_BB = table2array(traces(:, startsWith(names, 'v_stimBB_subj')));
t_BB = table2array(traces(:, startsWith(names, 't_stimBB_subj')));
b_BB = table2array(traces(:, startsWith(names, 'zt_stimBB_subj')));
z_trace = table2array(traces(:, strcmp(names, 'z')));
qS_trace = table2array(traces(:, startsWith(names, 'qS')));
qL_trace = table2array(traces(:, startsWith(names, 'qL')));

%% Make table
hssm_tbl = table('Size', [height(full_tbl),1], 'VariableTypes', ["double"], 'VariableNames',["rt"]);
hssm_tbl.rt = full_tbl.RT;
hssm_tbl.response = full_tbl.frontChoice + 1;
hssm_tbl.light = full_tbl.light + 1;
hssm_tbl.feedback = full_tbl.performance;
hssm_tbl.omission = full_tbl.omission;
hssm_tbl.response(hssm_tbl.omission) = -1;
hssm_tbl.split_by = str2double(full_tbl.Session);
hssm_tbl.sex = strcmp(full_tbl.Sex, "F");
hssm_tbl.condition = ones(height(hssm_tbl), 1);
hssm_tbl.condition(strcmp(full_tbl.sessionType, 'on')) = 2;
hssm_tbl.condition(strcmp(full_tbl.sessionType, 'pr')) = 3;
hssm_tbl.stim = full_tbl.stim;
hssm_tbl.subj_idx = grp2idx(categorical(full_tbl.Rat));
hssm_tbl.protocol = str2double(extractBefore(full_tbl.Protocol, 3)) - 43;
hssm_tbl.rule = full_tbl.rule;
hssm_tbl.rt(hssm_tbl.response==0) = -hssm_tbl.rt(hssm_tbl.response==0);

%%
lapse_rate=0.01;
ndr = mean(hssm_tbl.omission);
light_seq = [0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1;
             1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1;
             1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1;
             0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1;
             0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1] + 1;
rule_seq = [2, 0, 1, 0, 2, 0, 1, 0;
            0, 2, 0, 1, 0, 2, 0, 1;
            1, 0, 2, 0, 1, 0, 2, 0;
            0, 1, 0, 2, 0, 1, 0, 2;
            0, 1, 0, 2, 0, 2, 0, 1] + 1;
T = height(traces);
N = height(hssm_tbl);
S = length(unique(hssm_tbl.split_by));

slices = 32;
slice = reshape(1:T, [], slices)';
RT = cell(slices,1);
v_full = cell(slices,1);
a_full = cell(slices,1);
choices = cell(slices,1);
correct_full = cell(slices,1);
rule = cell(slices,1);
rule_num = cell(slices,1);
Vs = cell(slices,1);
Vs_hat = cell(slices,1);
b_full = cell(slices,1);

includes = [1,1,1];
RTMats = cell(size(includes,1),1);
correctMats = cell(size(includes,1),1);
rng(626)
parpool(16)
for I=1:size(includes,1)
    include = includes(I,:);
    parfor s=1:slices
        inds = slice(s,:);
        a_temp = a_trace(inds,:);
        l_temp = l_trace(inds,:);
        n_temp = n_trace(inds,:);
        v_temp = v_trace(inds,:);
        b_temp = b_trace(inds,:);
        f_temp = f_trace(inds,:);
        s_temp = s_trace(inds,:);
        st_temp = st_trace(inds);
        z_temp = z_trace(inds);
        qs_temp = qS_trace(inds);
        ql_temp = qL_trace(inds);
        a_stim_temp = a_stim(inds,:);
        v_stim_temp = v_stim(inds,:);
        b_stim_temp = b_stim(inds,:);
        t_stim_temp = t_stim(inds,:);
        a_stimON_temp = a_ON(inds,:);
        v_stimON_temp = v_ON(inds,:);
        b_stimON_temp = b_ON(inds,:);
        t_stimON_temp = t_ON(inds,:);
        a_stimBB_temp = a_BB(inds,:);
        v_stimBB_temp = v_BB(inds,:);
        b_stimBB_temp = b_BB(inds,:);
        t_stimBB_temp = t_BB(inds,:);
        Vs_temp = zeros(N,T/slices,2);
        Vs_hat_temp = zeros(N,T/slices,2);
        rn_temp = zeros(T/slices,S);
        V = [0,0];
        V_hat = [0,0];
        v_t = zeros(N,T/slices);
        a_t = zeros(N,T/slices);
        b_t = zeros(N,T/slices);
        v_that = zeros(N,T/slices);
        b_that = zeros(N,T/slices);
        n_t = zeros(N,T/slices);
        choices_temp = zeros(T/slices,N);
        rule_temp = false(T/slices,N);
        rt_temp = zeros(T/slices, N);
        correct_temp = false(T/slices, N);
        Q = [qs_temp(1), qs_temp(1), ql_temp(1)];
        Q_hat = [qs_temp(1), qs_temp(1), ql_temp(1)];
        blockNum = 1;
        for j=1:T/slices
            ses = 0;
            for i=1:N
                if i == 1 || hssm_tbl.split_by(i) ~= hssm_tbl.split_by(i-1)
                    if i~= 1
                        rn_temp(j,ses) = rn_temp(j,ses) + (blockNum - 1) / 5;
                    end
                    Q = [qs_temp(j), qs_temp(j), ql_temp(j)];
                    Q_hat = [qs_temp(j), qs_temp(j), ql_temp(j)];
                    blockNum = 1;
                    ses = ses + 1;
                end
                if hssm_tbl.condition(i) == 2
                    alpha_t = a_temp(j, hssm_tbl.subj_idx(i)) + a_stimON_temp(j, hssm_tbl.subj_idx(i));
                    delta_t = v_temp(j, hssm_tbl.subj_idx(i)) + v_stimON_temp(j, hssm_tbl.subj_idx(i));
                    bias_t = b_temp(j, hssm_tbl.subj_idx(i)) + b_stimON_temp(j, hssm_tbl.subj_idx(i));
                    ndt_t = n_temp(j, hssm_tbl.subj_idx(i)) + t_stimON_temp(j, hssm_tbl.subj_idx(i));
                elseif hssm_tbl.condition(i) == 3
                    alpha_t = a_temp(j, hssm_tbl.subj_idx(i)) + a_stimBB_temp(j, hssm_tbl.subj_idx(i));
                    delta_t = v_temp(j, hssm_tbl.subj_idx(i)) + v_stimBB_temp(j, hssm_tbl.subj_idx(i));
                    bias_t = b_temp(j, hssm_tbl.subj_idx(i)) + b_stimBB_temp(j, hssm_tbl.subj_idx(i));
                    ndt_t = n_temp(j, hssm_tbl.subj_idx(i)) + t_stimBB_temp(j, hssm_tbl.subj_idx(i));
                    if hssm_tbl.stim(i) == 1
                        alpha_t = alpha_t + a_stim_temp(j, hssm_tbl.subj_idx(i));
                        delta_t = delta_t + v_stim_temp(j, hssm_tbl.subj_idx(i));
                        bias_t = bias_t + b_stim_temp(j, hssm_tbl.subj_idx(i));
                        ndt_t = ndt_t + t_stim_temp(j, hssm_tbl.subj_idx(i));
                    end
                else
                    alpha_t = a_temp(j, hssm_tbl.subj_idx(i));
                    delta_t = v_temp(j, hssm_tbl.subj_idx(i));
                    bias_t = b_temp(j, hssm_tbl.subj_idx(i));
                    ndt_t = n_temp(j, hssm_tbl.subj_idx(i));
                end
                lr_t = logit(l_temp(j, hssm_tbl.subj_idx(i)));
                forg_t = logit(f_temp(j, hssm_tbl.subj_idx(i)));
                surp_t = exp(s_temp(j, hssm_tbl.subj_idx(i)));
                n_t(i,j) = ndt_t;
                a_t(i,j) = alpha_t;
                if hssm_tbl.response(i) ~= -1
                    V(1) = Q(1) + (hssm_tbl.light(i) == 1) * (Q(3));
                    V(2) = Q(2) + (hssm_tbl.light(i) == 2) * (Q(3));
                    Vs_temp(i,j,:) = V;
                    v_t(i,j) = (V(1) - V(2)) * delta_t;
                    b_t(i,j) = logit((Q(1) - Q(2)) * bias_t + z_temp(j));
                    pe = hssm_tbl.feedback(i) - Q(hssm_tbl.response(i));
                    adj_lr = lr_t * abs(pe) ^ surp_t;
                    Q(hssm_tbl.response(i)) = Q(hssm_tbl.response(i)) + adj_lr * pe;
                    Q(abs(hssm_tbl.response(i)-3)) = Q(abs(hssm_tbl.response(i)-3)) * (1-forg_t);
                    if hssm_tbl.response(i) == hssm_tbl.light(i)
                        pe = hssm_tbl.feedback(i) - Q(3);
                        adj_lr = lr_t * abs(pe) ^ surp_t;
                        Q(3) = Q(3) + adj_lr * pe;
                    else
                        Q(3) = Q(3) * (1-forg_t);
                    end
                end
                if rand < ndr
                    choices_temp(j,i) = -1;
                    rt_temp(j,i) = -1;
                    correct_temp(j,i) = 0;
                else
                    V_hat(1) = Q_hat(1) + (light_seq(hssm_tbl.protocol(i),mod(rn_temp(j,ses)-1,8)*5+blockNum) == 1) * (Q_hat(3));
                    V_hat(2) = Q_hat(2) + (light_seq(hssm_tbl.protocol(i),mod(rn_temp(j,ses)-1,8)*5+blockNum) == 2) * (Q_hat(3));
                    Vs_hat_temp(i,j,:) = V_hat;
                    v_that(i,j) = (V_hat(1) - V_hat(2)) * delta_t;
                    b_that(i,j) = logit((Q_hat(1) - Q_hat(2)) * bias_t + z_temp(j));
                    if rand < lapse_rate
                        choices_temp(j,i) = -(round(rand)-2);
                        rt_temp(j,i) = unifrnd(0,3.005);
                    else
                        rt_temp(j,i) = wienerrng(a_t(i,j), n_t(i) + st_temp(j) * (rand - 0.5),b_that(i,j)*a_t(i,j),v_that(i,j));
                        choices_temp(j, i) = -((rt_temp(j,i)>0)-2);
                        rt_temp(j,i) = abs(rt_temp(j,i));
                    end
                    if choices_temp(j,i) == rule_seq(hssm_tbl.protocol(i), mod(rn_temp(j,ses)-1,8)+1) || (choices_temp(j,i) == light_seq(hssm_tbl.protocol(i), mod(rn_temp(j,ses)-1,8)*5+blockNum) && 3 == rule_seq(hssm_tbl.protocol(i), mod(rn_temp(j,ses)-1,8)+1))
                        correct_temp(j,i) = 1;
                    else
                        correct_temp(j,i) = 0;
                    end
                    pe = correct_temp(j,i) - Q_hat(choices_temp(j,i));
                    adj_lr = lr_t * abs(pe) ^ surp_t;
                    Q_hat(choices_temp(j,i)) = Q_hat(choices_temp(j,i)) + adj_lr * pe;
                    Q_hat(abs(choices_temp(j,i)-3)) = Q_hat(abs(choices_temp(j,i)-3)) * (1-forg_t);
                    if choices_temp(j,i) == light_seq(hssm_tbl.protocol(i), mod(rn_temp(j,ses)-1, 8)*5 + blockNum)
                        pe = correct_temp(j,i) - Q_hat(3);
                        adj_lr = lr_t * abs(pe) ^ surp_t;
                        Q_hat(3) = Q_hat(3) + adj_lr * pe;
                    else
                        Q_hat(3) = Q_hat(3) * (1-forg_t);
                    end
                end
                rule_temp(j,i) = rule_seq(hssm_tbl.protocol(i), mod(rn_temp(j,ses)-1,8)+1) > 2;
                if correct_temp(j,i)
                    blockNum = blockNum + 1;
                    if blockNum == 6
                        rn_temp(j,ses) = rn_temp(j,ses) + 1;
                    end
                else
                    blockNum = 1;
                end
                blockNum = mod(blockNum - 1, 5) + 1;
            end
            rn_temp(j,end) = rn_temp(j,end) + (blockNum - 1) / 5;
            disp(j)
        end
        RT{s} = rt_temp;
        rule_num{s} = rn_temp;
        choices{s} = choices_temp;
        correct_full{s} = correct_temp;
        v_full{s} = v_t;
        a_full{s} = a_t;
        b_full{s} = b_t;
        Vs{s} = Vs_temp;
        Vs_hat{s} = Vs_hat_temp;
        rule{s} = rule_temp;
    end
    RTMats{I} = cell2mat(cellfun(@(x) transpose(x),RT,'UniformOutput',false)');
    correctMats{I} = cell2mat(cellfun(@(x) transpose(x),correct_full,'UniformOutput',false)');
end
delete(gcp('nocreate'))
a_full_mat = cell2mat(a_full');
b_full_mat = cell2mat(b_full');
RT_mat=cell2mat(cellfun(@(x) transpose(x),RT,'UniformOutput',false)');
dV_full_mat = cell2mat(cellfun(@(x) x(:,:,1)',Vs,'UniformOutput',false)) - cell2mat(cellfun(@(x) x(:,:,2)',Vs,'UniformOutput',false));
dV_hat_full_mat = cell2mat(cellfun(@(x) x(:,:,1)',Vs_hat,'UniformOutput',false)) - cell2mat(cellfun(@(x) x(:,:,2)',Vs_hat,'UniformOutput',false));
choices_mat = cell2mat(cellfun(@(x) transpose(x),choices,'UniformOutput',false)');
rule_num_mat = cell2mat(rule_num);
correct_full_mat = cell2mat(cellfun(@(x) transpose(x),correct_full,'UniformOutput',false)');
rule_mat = cell2mat(cellfun(@(x) transpose(x),rule,'UniformOutput',false)');

%% H
subplot_tight(8,9,[43,44,45,52,53,54,61,62,63,70,71,72])
mid_stds=[traces.a_std,traces.a_std,traces.v_std,traces.v_std];
mid_stds=sqrt((mid_stds.^2+[traces.a_stimON_std,traces.a_stim_std,traces.v_stimON_std,traces.v_stim_std].^2)/2);
mid_params=[traces.a_stimON,traces.a_stim,traces.v_stimON,traces.v_stim];
normed = mid_params./mid_stds;
al_goodplot2({normed(:,1)',normed(:,2)',normed(:,3)',normed(:,4)'}, 'pos', [1,1,2,2], 'type', {'left','right','left','right'},'boxw',0.4, 'col',[0.949,0.631,0.008;1,0,0;0.949,0.631,0.008;1,0,0],'useMedian', true)
patch('XData',[0,0,3,3],'YData',effect_size*[-1,1,1,-1],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.3)
ylim([-2,3])
xticks([1,2])
xlim([0.5,2.5])
xticklabels(["a", "v"])
ylabel('Effect Size')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

for i=1:4
    disp(nnz(abs(normed(:,i))<effect_size)/length(normed(:,i))*100)
    disp(nnz(normed(:,i)>0)/length(normed(:,i))*100)
    disp(median(normed(:,i)))
    disp(' ')
end

%%
mid_params=[traces.a_stimON,traces.v_stimON,traces.zt_stimON];
mid_stds=[traces.a_std,traces.v_std,traces.zt_std];
mid_stds=sqrt((mid_stds.^2+[traces.a_stimON_std,traces.v_stimON_std,traces.zt_stimON_std].^2)/2);
effect_size=0.1;
figure
hold on
for i=1:size(mid_params,2)
    normed=mid_params(:,i)./mid_stds(:,i);
    plot_95kd(normed,'c',[247,119,110]/255,@(x)0.9*x/max(x)-(i-1))
    sorted=sort(normed);
end
patch('XData',effect_size*[-1,1,1,-1],'YData',[1,1,1-size(mid_params,2),1-size(mid_params,2)],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.3)
yticks(1-size(mid_params,2):0)
yticklabels(flip(["Boundary Separation","Drift Rate", "Bias"]))
ylim([0.5-size(mid_params,2), 1])
xlim([-2,2])
xticks(-2:2)
xlabel('Effect Size')
set(gca,'fontsize',18)

%%
figure
points = -0.75:0.001:0.75;
choice_curves = zeros(length(points),T);
choice_curves_hat = zeros(length(points),T);
rng(626)
for i=1:T
    data_choice=hssm_tbl.response(hssm_tbl.response>0);
    [temp,I]=sort(dV_full_mat(i,hssm_tbl.response>0)+rand(1,nnz(hssm_tbl.response>0))*10e-6);
    choice_curves(:,i)=interp1(temp,movmean(data_choice(I),0.1,'samplepoints',temp),points);
    hat_rt=choices_mat(RT_mat(:,i)>0,i);
    [temp,I]=sort(dV_hat_full_mat(i,RT_mat(:,i)>0)+rand(1,nnz(RT_mat(:,i)>0))*10e-6);
    choice_curves_hat(:,i)=interp1(temp,movmean(hat_rt(I),0.1,'samplepoints',temp),points);
    if mod(i,100)==0
        disp(i)
    end
end
mcurve = 100*nanmedian(-choice_curves+2,2);
sor_curve = 100*(sort(-choice_curves+2,2));
lcurve = sor_curve(:,0.025*4000);
hcurve = sor_curve(:, 0.975*4000);
mcurve2=100*nanmedian(-choice_curves_hat+2,2);
sor_curve2 = 100*(sort(-choice_curves_hat+2,2));
lcurve2 = sor_curve2(:,0.025*4000);
hcurve2 = sor_curve2(:, 0.975*4000);
hold on
patch('XData',[points,flip(points)], 'YData', [lcurve;flip(hcurve)],'FaceColor',[180,0,255]/255,'EdgeColor','none','FaceAlpha',0.3)
patch('XData',[points,flip(points)], 'YData', [lcurve2;flip(hcurve2)],'FaceColor',[255,154,0]/255,'EdgeColor','none','FaceAlpha',0.3)
plot(points,mcurve,'-','Color',[180,0,255]/255,'LineWidth',2)
plot(points,mcurve2,'-','Color',[255,154,0]/255,'LineWidth',2)
xlabel('\DeltaV')
ylabel('% Left Choice')
xlim([points(1),points(end)])
xticks(-1:0.5:1)
yticks(0:50:100)
ylim([0,100])
set(gca,'fontsize',18)

%%
figure
points = -0.75:0.001:0.75;
choice_curves = zeros(length(points),T);
choice_curves_hat = zeros(length(points),T);
rng(626)
for i=1:T
    data_rt=hssm_tbl.rt(hssm_tbl.response>0);
    [temp,I]=sort(dV_full_mat(i,hssm_tbl.response>0)+rand(1,nnz(hssm_tbl.response>0))*10e-6);
    choice_curves(:,i)=interp1(temp,movmean(data_rt(I),0.1,'samplepoints',temp),points);
    hat_rt=RT_mat(RT_mat(:,i)>0,i);
    [temp,I]=sort(dV_hat_full_mat(i,RT_mat(:,i)>0)+rand(1,nnz(RT_mat(:,i)>0))*10e-6);
    choice_curves_hat(:,i)=interp1(temp,movmean(hat_rt(I),0.1,'samplepoints',temp),points);
    if mod(i,100)==0
        disp(i)
    end
end
mcurve = nanmedian(choice_curves,2);
sor_curve = (sort(choice_curves,2));
lcurve = sor_curve(:,0.025*4000);
hcurve = sor_curve(:, 0.975*4000);
mcurve2=nanmedian(choice_curves_hat,2);
sor_curve2 = (sort(choice_curves_hat,2));
lcurve2 = sor_curve2(:,0.025*4000);
hcurve2 = sor_curve2(:, 0.975*4000);
hold on
patch('XData',[points,flip(points)], 'YData', [lcurve;flip(hcurve)],'FaceColor',[180,0,255]/255,'EdgeColor','none','FaceAlpha',0.3)
patch('XData',[points,flip(points)], 'YData', [lcurve2;flip(hcurve2)],'FaceColor',[255,154,0]/255,'EdgeColor','none','FaceAlpha',0.3)
plot(points,mcurve,'-','Color',[180,0,255]/255,'LineWidth',2)
plot(points,mcurve2,'-','Color',[255,154,0]/255,'LineWidth',2)
xlabel('\DeltaV')
ylabel('RT (s)')
xlim([points(1),points(end)])
xticks(-1:0.5:1)
yticks(0.55:0.1:0.75)
ylim([0.55,0.75])
set(gca,'fontsize',18)

%%
figure
hold on
plot_95kd(rule_num_mat(:), 'c',[247,119,110]/255)
yticks([])
plot_95kd(ses_tbl.rulesCompleted, 'c',[247,119,110]/255)
xlabel('Rules Completed')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%%
figure
hold on
pts = 0:0.001:3;
% stims = hssm_tbl.condition==2;
% no_stims = hssm_tbl.condition==1;
stims = hssm_tbl.condition==3&hssm_tbl.stim==1;
no_stims = hssm_tbl.condition==3&hssm_tbl.stim==0;
bw = 0.05;
for s = 1:slices
    for i=1:50:T/slices
        temp = RT{s}(i,:);
        temp1 = squeeze(temp(stims));
        temp1=rmmissing(temp1(temp1>0));
        temp2 = squeeze(temp(no_stims));
        temp2=rmmissing(temp2(temp2>0));
        if ~isempty(temp1)
            vals = ksdensity(temp2,pts,'bandwidth',bw);
            plot(pts, vals, "Color",[0.6,0.6,0.6,0.05])
            vals = ksdensity(temp1,pts,'bandwidth',bw);
            plot(pts, vals, "Color",[0.949,0.631,0.008,0.05])
        end
    end
end
hold on
vals=ksdensity(squeeze(hssm_tbl.rt(no_stims&hssm_tbl.response>-1)),pts,'bandwidth',bw);
plot(pts, vals, "Color",[0.6,0.6,0.6],'LineWidth',1.5)
vals=ksdensity(squeeze(hssm_tbl.rt(stims&hssm_tbl.response>-1)),pts,'bandwidth',bw);
plot(pts, vals, "Color",[0.949,0.631,0.008],'LineWidth',1.5)
xticks(0:3)
xlim([0,3])
xlabel('RT (s)')
yticks([])
set(gca,'fontsize',18)
set(gca,'linewidth',2)
%%

mid_params=[traces.a_stim,traces.v_stim,traces.zt_stim,traces.t_stim];
mid_stds=[traces.a_std,traces.v_std,traces.zt_std,traces.t_std];
mid_stds=sqrt((mid_stds.^2+[traces.a_stim_std,traces.v_stim_std,traces.zt_stim_std,traces.t_stim_std].^2)/2);
effect_size=0.1;
figure
hold on
for i=1:size(mid_params,2)
    normed=mid_params(:,i)./mid_stds(:,i);
    plot_95kd(normed,'c',[247,119,110]/255,@(x)0.9*x/max(x)-(i-1))
    sorted=sort(normed);
    disp(nnz(abs(sorted)<effect_size)/length(sorted)*100)
end
patch('XData',effect_size*[-1,1,1,-1],'YData',[1,1,1-size(mid_params,2),1-size(mid_params,2)],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.3)
yticks(1-size(mid_params,2):0)
yticklabels(flip(["Boundary Separation","Drift Rate", "Bias","Non-Decision Time"]))
ylim([0.5-size(mid_params,2), 1])
xlim([-2,2])
xticks(-2:2)
xlabel('Effect Size')
set(gca,'fontsize',18)

%% Supplement RLDDM
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% A
subplot_tight(2,2,1)
hold on
pts = 0:0.001:3;
stims = hssm_tbl.condition==2;
no_stims = hssm_tbl.condition==1;
% stims = hssm_tbl.condition==3&hssm_tbl.stim==1;
%no_stims = hssm_tbl.condition==3&hssm_tbl.stim==0;
bw = 0.05;
for s = 1:slices
    for i=1:50:T/slices
        temp = RT{s}(i,:);
        temp1 = squeeze(temp(stims));
        temp1=rmmissing(temp1(temp1>0));
        temp2 = squeeze(temp(no_stims));
        temp2=rmmissing(temp2(temp2>0));
        if ~isempty(temp1)
            vals = ksdensity(temp2,pts,'bandwidth',bw);
            plot(pts, vals, "Color",[0.6,0.6,0.6,0.05])
            vals = ksdensity(temp1,pts,'bandwidth',bw);
            plot(pts, vals, "Color",[0.949,0.631,0.008,0.05])
        end
    end
end
hold on
vals=ksdensity(squeeze(hssm_tbl.rt(no_stims&hssm_tbl.response>-1)),pts,'bandwidth',bw);
plot(pts, vals, "Color",[0.6,0.6,0.6],'LineWidth',1.5)
vals=ksdensity(squeeze(hssm_tbl.rt(stims&hssm_tbl.response>-1)),pts,'bandwidth',bw);
plot(pts, vals, "Color",[0.949,0.631,0.008],'LineWidth',1.5)
xticks(0:3)
xlim([0,3])
xlabel('RT (s)')
yticks([])
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% B
subplot_tight(2,2,2)
hold on
pts = 0:0.001:3;
stims = hssm_tbl.condition==3&hssm_tbl.stim==1;
no_stims = hssm_tbl.condition==3&hssm_tbl.stim==0;
bw = 0.05;
for s = 1:slices
    for i=1:50:T/slices
        temp = RT{s}(i,:);
        temp1 = squeeze(temp(stims));
        temp1=rmmissing(temp1(temp1>0));
        temp2 = squeeze(temp(no_stims));
        temp2=rmmissing(temp2(temp2>0));
        if ~isempty(temp1)
            vals = ksdensity(temp2,pts,'bandwidth',bw);
            plot(pts, vals, "Color",[0,0,1,0.05])
            vals = ksdensity(temp1,pts,'bandwidth',bw);
            plot(pts, vals, "Color",[1,0,0,0.05])
        end
    end
end
hold on
vals=ksdensity(squeeze(hssm_tbl.rt(no_stims&hssm_tbl.response>-1)),pts,'bandwidth',bw);
plot(pts, vals, "Color",[0,0,1],'LineWidth',1.5)
vals=ksdensity(squeeze(hssm_tbl.rt(stims&hssm_tbl.response>-1)),pts,'bandwidth',bw);
plot(pts, vals, "Color",[1,0,0],'LineWidth',1.5)
xticks(0:3)
xlim([0,3])
xlabel('RT (s)')
yticks([])
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%% C
subplot_tight(2,2,3)
points = -0.75:0.001:0.75;
choice_curves = zeros(length(points),T);
choice_curves_hat = zeros(length(points),T);
rng(626)
for i=1:T
    data_choice=hssm_tbl.response(hssm_tbl.response>0);
    [temp,I]=sort(dV_full_mat(i,hssm_tbl.response>0)+rand(1,nnz(hssm_tbl.response>0))*10e-6);
    choice_curves(:,i)=interp1(temp,movmean(data_choice(I),0.1,'samplepoints',temp),points);
    hat_rt=choices_mat(RT_mat(:,i)>0,i);
    [temp,I]=sort(dV_hat_full_mat(i,RT_mat(:,i)>0)+rand(1,nnz(RT_mat(:,i)>0))*10e-6);
    choice_curves_hat(:,i)=interp1(temp,movmean(hat_rt(I),0.1,'samplepoints',temp),points);
    if mod(i,100)==0
        disp(i)
    end
end
mcurve = 100*nanmedian(-choice_curves+2,2);
sor_curve = 100*(sort(-choice_curves+2,2));
lcurve = sor_curve(:,0.025*4000);
hcurve = sor_curve(:, 0.975*4000);
mcurve2=100*nanmedian(-choice_curves_hat+2,2);
sor_curve2 = 100*(sort(-choice_curves_hat+2,2));
lcurve2 = sor_curve2(:,0.025*4000);
hcurve2 = sor_curve2(:, 0.975*4000);
hold on
patch('XData',[points,flip(points)], 'YData', [lcurve;flip(hcurve)],'FaceColor',[180,0,255]/255,'EdgeColor','none','FaceAlpha',0.3)
patch('XData',[points,flip(points)], 'YData', [lcurve2;flip(hcurve2)],'FaceColor',[255,154,0]/255,'EdgeColor','none','FaceAlpha',0.3)
plot(points,mcurve,'-','Color',[180,0,255]/255,'LineWidth',2)
plot(points,mcurve2,'-','Color',[255,154,0]/255,'LineWidth',2)
xlabel('\DeltaV')
ylabel('% Left Choice')
xlim([points(1),points(end)])
xticks(-1:0.5:1)
yticks(0:50:100)
ylim([0,100])
set(gca,'fontsize',18)
set(gca,'linewidth',2)


%% D
subplot_tight(2,2,4)
mid_stds=[traces.zt_std,traces.zt_std,traces.t_std,traces.t_std];
mid_stds=sqrt((mid_stds.^2+[traces.zt_stimON_std,traces.zt_stim_std,traces.t_stimON_std,traces.t_stim_std].^2)/2);
mid_params=[traces.zt_stimON,traces.zt_stim,traces.t_stimON,traces.t_stim];
normed = mid_params./mid_stds;
effect_size=0.1;
al_goodplot2({normed(:,1)',normed(:,2)',normed(:,3)',normed(:,4)'}, 'pos', [1,1,2,2], 'type', {'left','right','left','right'},'boxw',0.4, 'col',[0.949,0.631,0.008;1,0,0;0.949,0.631,0.008;1,0,0],'useMedian', true)
patch('XData',[0,0,3,3],'YData',effect_size*[-1,1,1,-1],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.3)
ylim([-2,2])
xticks([1,2])
xlim([0.5,2.5])
xticklabels(["B", "t"])
ylabel('Effect Size')
set(gca,'fontsize',18)
set(gca,'linewidth',2)


%% Supplement Behavior
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% A
subplot_tight(2,2,1)
groups = {full_tbl.DT(strcmp(full_tbl.sessionType,'off')), ...
          full_tbl.DT(strcmp(full_tbl.sessionType,'on')), ...
          full_tbl.DT(strcmp(full_tbl.sessionType,'pr')&full_tbl.stim==0), ...
          full_tbl.DT(strcmp(full_tbl.sessionType,'pr')&full_tbl.stim==1)};
al_goodplot2(groups', 'pos', [1,2,3,3], 'type', {'bilateral','bilateral','left','right'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0,0,1;1,0,0], 'useMedian', true)
plot([ones(size(rat_tbl.RT_off)), 2*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.DT_off, rat_tbl.DT_on, nan(size(rat_tbl.RT_off))]', 'Color', [0.6,0.6,0.6], 'LineStyle','--')
plot([2.9*ones(size(rat_tbl.RT_off)), 3.1*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.DT_pr, rat_tbl.DT_pr0, nan(size(rat_tbl.RT_off))]', 'Color', [0.6,0.6,0.6], 'LineStyle','--')
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.DT_off, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.DT_on, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)) - 0.1,rat_tbl.DT_pr, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(3*ones(size(rat_tbl.RT_off)) + 0.1,rat_tbl.DT_pr0, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])

if coefTest(dtglm, [0,0,1,0,0]) < 0.05 / 5
    plot([1,2], [6,6], 'k', 'LineWidth', 2)
    text(1.5, 6.15, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(dtglm, [0,0,0,1,0]) < 0.05 / 5
    plot([1,2.9], [8,8], 'k', 'LineWidth', 2)
    text(2.05, 8.3, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(dtglm, [0,1,0,0,0]) < 0.05 / 5
    plot([2.9,3.1], [7,7], 'k', 'LineWidth', 2)
    text(3, 7.2, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(dtglm, [0,1,-1,1,0]) < 0.05 / 5
    plot([2,3.1], [5,5], 'k', 'LineWidth', 2)
    text(2.55, 5.1, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:3)
xticklabels(["OFF", "ON", "BB"])
ylabel('DT (s)')
set(gca, 'YScale', 'log')
ylim([0.4,16])
yticks([0.5,1,2,4,8,16])
xlim([0.5,3.5])
set(gca, 'fontsize', 18)
set(gca, 'YMinorTick', 'off')

%% B
subplot_tight(2,2,2)
groups = {ses_tbl.omission(strcmp(ses_tbl.sessionType,'off'))*100, ...
          ses_tbl.omission(strcmp(ses_tbl.sessionType,'on'))*100, ...
          ses_tbl.omission(strcmp(ses_tbl.sessionType,'pr'))*100, ...
          ses_tbl.omission0(strcmp(ses_tbl.sessionType,'pr'))*100};
al_goodplot2(groups', 'pos', [1,2,3,3], 'type', {'bilateral','bilateral','left','right'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0,0,1;1,0,0])
plot([ones(size(rat_tbl.RT_off)), 2*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.omission_off, rat_tbl.omission_on, nan(size(rat_tbl.RT_off))]'*100, 'Color', [0.6,0.6,0.6], 'LineStyle','--')
plot([2.9*ones(size(rat_tbl.RT_off)), 3.1*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.omission_pr, rat_tbl.omission_pr0, nan(size(rat_tbl.RT_off))]'*100, 'Color', [0.6,0.6,0.6], 'LineStyle','--')
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.omission_off*100, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.omission_on*100, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)) - 0.1,rat_tbl.omission_pr*100, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(3*ones(size(rat_tbl.RT_off)) + 0.1,rat_tbl.omission_pr0*100, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])

if coefTest(omglm, [0,0,1,0]) < 0.05 / 5
    plot([1,2], [17.5,17.5], 'k', 'LineWidth', 2)
    text(1.5, 17.6, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(omglm, [0,0,0,1]) < 0.05 / 5
    plot([1,2.9], [17,17], 'k', 'LineWidth', 2)
    text(2.05, 17.1, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(omglm, [0,1,0,0]) < 0.05 / 5
    plot([2.9,3.1], [17,17], 'k', 'LineWidth', 2)
    text(3, 17.1, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(omglm, [0,1,-1,1]) < 0.05 / 5
    plot([2,3.1], [16.5,16.5], 'k', 'LineWidth', 2)
    text(2.55, 16.6, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:3)
xticklabels(["OFF", "ON", "BB"])
ylabel('% Omissions')
ylim([0,20])
xlim([0.5,3.5])
set(gca, 'fontsize', 18)

%% D
subplot_tight(2,2,4)
groups = {ses_tbl.frontChoice(strcmp(ses_tbl.sessionType,'off'))*100, ...
          ses_tbl.frontChoice(strcmp(ses_tbl.sessionType,'on'))*100, ...
          ses_tbl.frontChoice(strcmp(ses_tbl.sessionType,'pr'))*100, ...
          ses_tbl.frontChoice0(strcmp(ses_tbl.sessionType,'pr'))*100};
al_goodplot2(groups', 'pos', [1,2,3,3], 'type', {'bilateral','bilateral','left','right'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0,0,1;1,0,0])
plot([ones(size(rat_tbl.RT_off)), 2*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.frontChoice_off, rat_tbl.frontChoice_on, nan(size(rat_tbl.RT_off))]'*100, 'Color', [0.6,0.6,0.6], 'LineStyle','--')
plot([2.9*ones(size(rat_tbl.RT_off)), 3.1*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.frontChoice_pr, rat_tbl.frontChoice_pr0, nan(size(rat_tbl.RT_off))]'*100, 'Color', [0.6,0.6,0.6], 'LineStyle','--')
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.frontChoice_off*100, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.frontChoice_on*100, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)) - 0.1,rat_tbl.frontChoice_pr*100, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(3*ones(size(rat_tbl.RT_off)) + 0.1,rat_tbl.frontChoice_pr0*100, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])

if coefTest(fcglm, [0,1,0,0,0]) < 0.05 / 5
    plot([2.9,3.1], [63,63], 'k', 'LineWidth', 2)
    text(3, 64, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(fcglm, [0,0,0,1,0]) < 0.05 / 5
    plot([1,2.9], [63,63], 'k', 'LineWidth', 2)
    text(2.05, 64, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(fcglm, [0,0,1,0,0]) < 0.05 / 5
    plot([1,2], [66,66], 'k', 'LineWidth', 2)
    text(1.5, 0.67, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(fcglm, [0,1,-1,1,0]) < 0.05 / 5
    plot([2,3.1], [6,6], 'k', 'LineWidth', 2)
    text(2.55, 61, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(fcglm, [1,0,0,0,0]) < 0.05 / 5
    text(1, 68, '#', 'HorizontalAlignment','center','fontsize', 20)
end
if coefTest(fcglm, [1,0,0,1,0]) < 0.05 / 5
    text(2.9, 68, '#', 'HorizontalAlignment','center','fontsize', 20)
end

xticks(1:3)
xticklabels(["OFF", "ON", "BB"])
ylabel('% Left Choice')
ylim([40,70])
yticks([40,50,60,70])
xlim([0.5,3.5])
set(gca, 'fontsize', 18)

%%
subplot_tight(2,2,4)
groups = {progress_tbl.progress(strcmp(progress_tbl.sessionType,'off')), ...
          progress_tbl.progress(strcmp(progress_tbl.sessionType,'on')), ...
          progress_tbl.progress(strcmp(progress_tbl.sessionType,'pr')&progress_tbl.stim==0), ...
          progress_tbl.progress(strcmp(progress_tbl.sessionType,'pr')&progress_tbl.stim==1)};
al_goodplot2(groups', 'pos', [1,2,3,3], 'type', {'bilateral','bilateral','left','right'},'boxw',0.4, 'col',[0.6,0.6,0.6;0.949,0.631,0.008;0,0,1;1,0,0])
plot([ones(size(rat_tbl.RT_off)), 2*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.progress_off, rat_tbl.progress_on, nan(size(rat_tbl.RT_off))]', 'Color', [0.6,0.6,0.6], 'LineStyle','--')
plot([2.9*ones(size(rat_tbl.RT_off)), 3.1*ones(size(rat_tbl.RT_off)), ones(size(rat_tbl.RT_off))]', [rat_tbl.progress_pr, rat_tbl.progress_pr0, nan(size(rat_tbl.RT_off))]', 'Color', [0.6,0.6,0.6], 'LineStyle','--')
scatter(ones(size(rat_tbl.RT_off)),rat_tbl.progress_off, 'MarkerFaceColor', [0.6,0.6,0.6], 'MarkerEdgeColor', [0.6,0.6,0.6])
scatter(2*ones(size(rat_tbl.RT_off)),rat_tbl.progress_on, 'MarkerFaceColor', [0.949,0.631,0.008], 'MarkerEdgeColor', [0.949,0.631,0.008])
scatter(3*ones(size(rat_tbl.RT_off)) - 0.1,rat_tbl.progress_pr, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(3*ones(size(rat_tbl.RT_off)) + 0.1,rat_tbl.progress_pr0, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])

if coefTest(progglm, [0,1,0,0]) < 0.05 / 6
    plot([1,2.9], [0.5,0.5], 'k', 'LineWidth', 2)
    text(2, 0.53, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:3)
xticklabels(["OFF", "ON", "BB"])
ylabel('Rule Progress')
xlim([0.5,3.5])
set(gca, 'fontsize', 18)