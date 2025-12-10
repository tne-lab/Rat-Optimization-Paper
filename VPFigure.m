addpath(genpath('C:\Users\Evan\Documents\GitHub\eelib'))

%% Load Data
full_tbl = load_pr_data("Data/VP", @read_ss_pyb_vp);
full_tbl.omission = full_tbl.RT>3;
full_tbl.DT = full_tbl.initTime - full_tbl.cueTime;
full_tbl = full_tbl(~strcmp(full_tbl.Rat, 'FORT07'), :);
full_tbl = full_tbl(~strcmp(full_tbl.Rat, 'MORT09'),:);

%% Add Rule info
trial = 1;
prevRule = "";
prevSes = "";
progress = 0;
for i=1:height(full_tbl)
    if ~strcmp(prevRule, full_tbl.rule(i)) || ~strcmp(prevSes,full_tbl.Session(i))
        progress = 0;
        trial = 1;
        prevRule = full_tbl.rule(i);
        prevSes = full_tbl.Session(i);
    end
    if full_tbl.performance(i) == 1
        progress = progress + 1;
        full_tbl.progress(i) = 1;
    else
        full_tbl.progress(i) = -progress;
        progress = 0;
    end
    full_tbl.tinrule(i) = trial;
    trial = trial + 1;
end

%% Create session table
sessions = unique(full_tbl.Session);
ses_tbl = table('Size', [length(sessions), 9], 'VariableTypes', ["string","string","double","double","double","double","double","double","double"], 'VariableNames',["Rat","sessionType", "performance", "omission", "frontChoice","performance0", "omission0", "frontChoice0","rulesCompleted"]);
for i=1:length(sessions)
    sub_tbl = full_tbl(full_tbl.Session==sessions(i),:);
    ses_tbl.Rat(i) = sub_tbl.Rat(1);
    ses_tbl.protocolCode(i) = extractBefore(sub_tbl.Protocol(1),3);
    rule = strcmp(sub_tbl.rule, 'L');
    ses_tbl.progress(i) = mean(sub_tbl.progress(strcmp(sub_tbl.stim,'0')));
    ses_tbl.performance(i) = mean(sub_tbl.performance(sub_tbl.RT < 3&strcmp(sub_tbl.stim,'0')));
    ses_tbl.omission(i) = mean(sub_tbl.RT(strcmp(sub_tbl.stim,'0')) > 3);
    ses_tbl.frontChoice(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3 &strcmp(sub_tbl.stim,'0')));
    ses_tbl.performance1(i) = mean(sub_tbl.performance(sub_tbl.RT < 3&strcmp(sub_tbl.stim,'1')));
    ses_tbl.progress1(i) = mean(sub_tbl.progress(strcmp(sub_tbl.stim,'1')));
    ses_tbl.omission1(i) = mean(sub_tbl.RT(strcmp(sub_tbl.stim,'1')) > 3);
    ses_tbl.frontChoice1(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3 &strcmp(sub_tbl.stim,'1')));
    ses_tbl.progress2(i) = mean(sub_tbl.progress(strcmp(sub_tbl.stim,'2')));
    ses_tbl.performance2(i) = mean(sub_tbl.performance(sub_tbl.RT < 3&strcmp(sub_tbl.stim,'2')));
    ses_tbl.omission2(i) = mean(sub_tbl.RT(strcmp(sub_tbl.stim,'2')) > 3);
    ses_tbl.frontChoice2(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3 &strcmp(sub_tbl.stim,'2')));
    ses_tbl.progress3(i) = mean(sub_tbl.progress(strcmp(sub_tbl.stim,'3')));
    ses_tbl.performance3(i) = mean(sub_tbl.performance(sub_tbl.RT < 3&strcmp(sub_tbl.stim,'3')));
    ses_tbl.omission3(i) = mean(sub_tbl.RT(strcmp(sub_tbl.stim,'3')) > 3);
    ses_tbl.frontChoice3(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3 &strcmp(sub_tbl.stim,'3')));
    ses_tbl.progressB(i) = mean(sub_tbl.progress(strcmp(sub_tbl.stim,'B')));
    ses_tbl.performanceB(i) = mean(sub_tbl.performance(sub_tbl.RT < 3&strcmp(sub_tbl.stim,'B')));
    ses_tbl.omissionB(i) = mean(sub_tbl.RT(strcmp(sub_tbl.stim,'B')) > 3);
    ses_tbl.frontChoiceB(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3 &strcmp(sub_tbl.stim,'B')));
end

progress_tbl = table('Size', [length(sessions)*5, 1], 'VariableTypes', ["string"], 'VariableNames', ["Rat"]);
j = 1;
for i=1:height(ses_tbl)
    progress_tbl.Rat(j:j+4) = ses_tbl.Rat(i);
    progress_tbl.stim(j) = '0';
    progress_tbl.stim(j+1) = '1';
    progress_tbl.stim(j+2) = '2';
    progress_tbl.stim(j+3) = '3';
    progress_tbl.stim(j+4) = 'B';
    progress_tbl.progress(j) = ses_tbl.progress(i);
    progress_tbl.progress(j+1) = ses_tbl.progress1(i);
    progress_tbl.progress(j+2) = ses_tbl.progress2(i);
    progress_tbl.progress(j+3) = ses_tbl.progress3(i);
    progress_tbl.progress(j+4) = ses_tbl.progressB(i);
    for k=0:4
        progress_tbl.protocolCode(j+k) = ses_tbl.protocolCode(i);
    end
    j = j + 5;
end

%% Create Rat Table
rats = unique(full_tbl.Rat);
rat_tbl = table('Size', [length(rats), 1], 'VariableTypes', ["string"], 'VariableNames',["Rat"]);
for i=1:length(rats)
    rat_tbl.Rat(i) = rats(i);
    conditions = ["0","1","2","3","B"];
    for j=1:length(conditions)
        sub_tbl = full_tbl(strcmp(full_tbl.Rat,rats(i))&strcmp(full_tbl.stim, conditions(j)),:);
        rat_tbl.(strcat('RT_', conditions(j)))(i) = median(sub_tbl.RT(sub_tbl.RT < 3));
        rat_tbl.(strcat('DT_', conditions(j)))(i) = median(sub_tbl.DT);
        rat_tbl.(strcat('performance_', conditions(j)))(i) = mean(sub_tbl.performance(sub_tbl.RT < 3));
        rat_tbl.(strcat('omission_', conditions(j)))(i) = mean(sub_tbl.RT > 3);
        rat_tbl.(strcat('frontChoice_', conditions(j)))(i) = mean(sub_tbl.frontChoice(sub_tbl.RT < 3));
        rat_tbl.(strcat('progress_', conditions(j)))(i) = mean(sub_tbl.progress);
    end
end

%% Run VP GLM
full_tbl.rule_type = strcmp(full_tbl.rule, 'L');
full_tbl.protocolCode = extractBefore(full_tbl.Protocol, 3);
vp_rtglm = fitglme(full_tbl(full_tbl.RT<3&full_tbl.RT>0.01,:), 'RT~1+stim+rule_type+(1|Rat)+(1|Session)+(1|protocolCode)','distribution','gamma','link','identity');
vp_accglm = fitglme(full_tbl(full_tbl.RT<3,:), 'performance~1+stim+rule_type+(1|Rat)+(1|Session)+(1|protocolCode)','distribution','binomial','fitmethod','laplace');
vp_omglm = fitglme(full_tbl(:,:), 'omission~1+stim+rule_type+(1|Rat)+(1|Session)+(1|protocolCode)','distribution','binomial','fitmethod','laplace');
vp_dtglm = fitglme(full_tbl, 'DT~1+stim+rule_type+(1|Rat)+(1|Session)+(1|protocolCode)','distribution','gamma','link','log');
vp_fcglm = fitglme(full_tbl(full_tbl.RT<3,:), 'frontChoice~1+stim+rule_type+(1|Rat)+(1|Session)+(1|protocolCode)','distribution','binomial');
progglm = fitglme(progress_tbl, 'progress~1+stim+(1|Rat)+(1|protocolCode)');

%%
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% A
subplot_tight(7,2,1)
stim_runs = full_tbl.stim(strcmp(full_tbl.Session, "1"));
run_length = 0;
condition = string(stim_runs(1));
for i=1:length(stim_runs)
    if ~strcmp(string(stim_runs(i)), condition)
        if strcmp(condition, '0')
            color = 'b';
        elseif strcmp(condition, '1')
            color = ([0,0,1] * 2 + [1,0,0])/3;
        elseif strcmp(condition, '2')
            color = ([0,0,1] + [1,0,0] * 2)/3;
        elseif strcmp(condition, '3')
            color = 'r';
        else
            color = [240 100 10]/256;
        end
        patch('XData', [i - run_length+0.3, i - run_length + 0.3, i - 0.3, i - 0.3], 'YData', [0, 1, 1, 0], 'FaceColor', color, 'FaceAlpha', 0.5,'EdgeColor','none')
        condition = stim_runs(i);
        run_length = 0;
    end
    run_length = run_length + 1;
end
if strcmp(condition, '0')
    color = 'b';
elseif strcmp(condition, '1')
    color = ([0,0,1] * 2 + [1,0,0])/3;
elseif strcmp(condition, '2')
    color = ([0,0,1] + [1,0,0] * 2)/3;
elseif strcmp(condition, '3')
    color = 'r';
else
    color = [0.94,0.39,0.04];
end
patch('XData', [i - run_length+0.3, i - run_length + 0.3, i - 0.3, i - 0.3], 'YData', [0, 1, 1, 0], 'FaceColor', color, 'FaceAlpha', 0.5,'EdgeColor','none')
axis tight
axis off
ylim([-1,1])

%% B
subplot_tight(7,2,[3,5,7])
groups = {full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.stim,'0')), ...
          full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.stim,'1')), ...
          full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.stim,'2')), ...
          full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.stim,'3')), ...
          full_tbl.RT(full_tbl.RT<3&strcmp(full_tbl.stim,'B'))};
c1 = ([0,0,1] * 2 + [1,0,0])/3;
c2 = ([0,0,1] + [1,0,0] * 2)/3;
al_goodplot2(groups', 'pos', [1,2,3,4,5], 'type', {'bilateral','bilateral','bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0,0,1;c1;c2;1,0,0;0.94,0.39,0.04],'useMedian', true)
scatter(ones(size(rat_tbl.RT_0)),rat_tbl.RT_0, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(2*ones(size(rat_tbl.RT_0)),rat_tbl.RT_1, 'MarkerFaceColor', c1, 'MarkerEdgeColor', c1)
scatter(3*ones(size(rat_tbl.RT_0)),rat_tbl.RT_2, 'MarkerFaceColor', c2, 'MarkerEdgeColor', c2)
scatter(4*ones(size(rat_tbl.RT_0)),rat_tbl.RT_3, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])
scatter(5*ones(size(rat_tbl.RT_0)),rat_tbl.RT_B, 'MarkerFaceColor', [0.94,0.39,0.04], 'MarkerEdgeColor', [0.94,0.39,0.04])

if coefTest(vp_rtglm, [0,0,1,0,0,0]) < 0.05
    plot([1,5], [1.1,1.1], 'k', 'LineWidth', 2)
    text(3, 1.11, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(vp_rtglm, [0,0,0,1,0,0]) < 0.05
    plot([1,4], [1.025,1.025], 'k', 'LineWidth', 2)
    text(2.5, 1.035, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(vp_rtglm, [0,0,0,0,1,0]) < 0.05
    plot([1,3], [0.95,0.95], 'k', 'LineWidth', 2)
    text(2, 0.96, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:5)
ylim([0.2,1.2])
xticklabels(["OFF", "100", "200", "300", "300@20"])
ylabel('RT (s)')
xlim([0.5,5.5])
set(gca, 'fontsize', 18)
set(gca,'linewidth',2)

%% C
subplot_tight(7,2,[4,6,8])
groups = {ses_tbl.performance, ...
          ses_tbl.performance1, ...
          ses_tbl.performance2, ...
          ses_tbl.performance3, ...
          ses_tbl.performanceB};
c1 = ([0,0,1] * 2 + [1,0,0])/3;
c2 = ([0,0,1] + [1,0,0] * 2)/3;
al_goodplot2(groups', 'pos', [1,2,3,4,5], 'type', {'bilateral','bilateral','bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0,0,1;c1;c2;1,0,0;0.94,0.39,0.04],'useMedian', true)
scatter(ones(size(rat_tbl.RT_0)),rat_tbl.performance_0, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(2*ones(size(rat_tbl.RT_0)),rat_tbl.performance_1, 'MarkerFaceColor', c1, 'MarkerEdgeColor', c1)
scatter(3*ones(size(rat_tbl.RT_0)),rat_tbl.performance_2, 'MarkerFaceColor', c2, 'MarkerEdgeColor', c2)
scatter(4*ones(size(rat_tbl.RT_0)),rat_tbl.performance_3, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])
scatter(5*ones(size(rat_tbl.RT_0)),rat_tbl.performance_B, 'MarkerFaceColor', [0.94,0.39,0.04], 'MarkerEdgeColor', [0.94,0.39,0.04])

if coefTest(vp_accglm, [0,0,0,0,1,0]) < 0.05
    plot([1,3], [0.85,0.85], 'k', 'LineWidth', 2)
    text(2, 0.86, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:5)
ylim([0.4,0.9])
xticklabels(["OFF", "100", "200", "300", "300@20"])
ylabel('Accuracy')
xlim([0.5,5.5])
set(gca,'linewidth',2)
set(gca, 'fontsize', 18)

%% D
ts = [-7,-6,-5,-4,-3,-2,-1];
ts2 = [1,2,3,4,5,6,7];
sessions = unique(full_tbl.Session);
blocks = nan(39*length(sessions),7,2);
conditions = strings(39*length(sessions),2);

ind = 1;
for session=sessions'
    sub_tbl = full_tbl(strcmp(full_tbl.Session, session), :);
    stim = string(sub_tbl.stim(1));
    tinb = 1;
    block = nan(1, 7);
    first_block= true;
    for j=1:height(sub_tbl)
        if ~strcmp(stim, string(sub_tbl.stim(j)))
            if first_block
                blocks(ind, end-(tinb - 2):end,1) = block(1:tinb-1);
                conditions(ind, 1) = stim;
            else
                blocks(ind-1, 1:tinb-1,2) = block(1:tinb-1);
                blocks(ind, end-(tinb - 2):end,1) = block(1:tinb-1);
                conditions(ind-1, 2) = stim;
                conditions(ind,1) = stim;
            end
            stim = string(sub_tbl.stim(j));
            tinb = 1;
            ind = ind + 1;
            first_block = false;
        end
        block(tinb) = sub_tbl.RT(j);
        tinb = tinb + 1;
    end
    blocks(ind-1, 1:tinb-1,2) = block(1:tinb-1);
end

colors = [0,0,1;1/3,0,2/3;2/3,0,1/3;1,0,0;0.94,0.39,0.04];
subplot_tight(7,3,[13,16,19])
hold on
plot(ts,nanmean(blocks(strcmp(conditions(:,1), '0')&strcmp(conditions(:,2), 'B'),:,1)),'Color',colors(1,:), 'LineWidth',2)
plot(ts,nanmean(blocks(strcmp(conditions(:,1), '1')&strcmp(conditions(:,2), 'B'),:,1)),'Color',colors(2,:), 'LineWidth',2)
plot(ts,nanmean(blocks(strcmp(conditions(:,1), '2')&strcmp(conditions(:,2), 'B'),:,1)),'Color',colors(3,:), 'LineWidth',2)
plot(ts,nanmean(blocks(strcmp(conditions(:,1), '3')&strcmp(conditions(:,2), 'B'),:,1)),'Color',colors(4,:), 'LineWidth',2)
plot(ts2,nanmean(blocks(strcmp(conditions(:,1), '0')&strcmp(conditions(:,2), 'B'),:,2)),'Color',colors(5,:), 'LineWidth',2)
plot(ts2,nanmean(blocks(strcmp(conditions(:,1), '1')&strcmp(conditions(:,2), 'B'),:,2)),'Color',colors(5,:), 'LineWidth',2)
plot(ts2,nanmean(blocks(strcmp(conditions(:,1), '2')&strcmp(conditions(:,2), 'B'),:,2)),'Color',colors(5,:), 'LineWidth',2)
plot(ts2,nanmean(blocks(strcmp(conditions(:,1), '3')&strcmp(conditions(:,2), 'B'),:,2)),'Color',colors(5,:), 'LineWidth',2)
ylabel('RT (s)')
set(gca, 'fontsize', 18)
set(gca,'linewidth',2)
xlabel('Trials Since Change')

%% RLDDM
chains = dir(strcat("Data/fRLDDM_vp"));
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
a_100 = table2array(traces(:, startsWith(names, 'a_stim100_subj')));
v_100 = table2array(traces(:, startsWith(names, 'v_stim100_subj')));
t_100 = table2array(traces(:, startsWith(names, 't_stim100_subj')));
b_100 = table2array(traces(:, startsWith(names, 'zt_stim100_subj')));
a_200 = table2array(traces(:, startsWith(names, 'a_stim200_subj')));
v_200 = table2array(traces(:, startsWith(names, 'v_stim200_subj')));
t_200 = table2array(traces(:, startsWith(names, 't_stim200_subj')));
b_200 = table2array(traces(:, startsWith(names, 'zt_stim200_subj')));
a_300 = table2array(traces(:, startsWith(names, 'a_stim300_subj')));
v_300 = table2array(traces(:, startsWith(names, 'v_stim300_subj')));
t_300 = table2array(traces(:, startsWith(names, 't_stim300_subj')));
b_300 = table2array(traces(:, startsWith(names, 'zt_stim300_subj')));
a_300_20 = table2array(traces(:, startsWith(names, 'a_stim300_20_subj')));
v_300_20 = table2array(traces(:, startsWith(names, 'v_stim300_20_subj')));
t_300_20 = table2array(traces(:, startsWith(names, 't_stim300_20_subj')));
b_300_20 = table2array(traces(:, startsWith(names, 'zt_stim300_20_subj')));
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
% hssm_tbl.sex = strcmp(full_tbl.Sex, "F");
hssm_tbl.condition = ones(size(full_tbl.stim));
hssm_tbl.condition(strcmp(full_tbl.stim, '1')) = 2;
hssm_tbl.condition(strcmp(full_tbl.stim, '2')) = 3;
hssm_tbl.condition(strcmp(full_tbl.stim, '3')) = 4;
hssm_tbl.condition(strcmp(full_tbl.stim, 'B')) = 5;
hssm_tbl.protocol = str2double(extractBefore(full_tbl.Protocol, 3)) - 43;
hssm_tbl.rule = full_tbl.rule;
hssm_tbl.subj_idx = grp2idx(categorical(full_tbl.Rat));

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
% parpool(16)
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
        a_stim100_temp = a_100(inds,:);
        v_stim100_temp = v_100(inds,:);
        b_stim100_temp = b_100(inds,:);
        t_stim100_temp = t_100(inds,:);
        a_stim200_temp = a_200(inds,:);
        v_stim200_temp = v_200(inds,:);
        b_stim200_temp = b_200(inds,:);
        t_stim200_temp = t_200(inds,:);
        a_stim300_temp = a_300(inds,:);
        v_stim300_temp = v_300(inds,:);
        b_stim300_temp = b_300(inds,:);
        t_stim300_temp = t_300(inds,:);
        a_stim300_20_temp = a_300_20(inds,:);
        v_stim300_20_temp = v_300_20(inds,:);
        b_stim300_20_temp = b_300_20(inds,:);
        t_stim300_20_temp = t_300_20(inds,:);
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
                    alpha_t = a_temp(j, hssm_tbl.subj_idx(i)) + a_stim100_temp(j, hssm_tbl.subj_idx(i));
                    delta_t = v_temp(j, hssm_tbl.subj_idx(i)) + v_stim100_temp(j, hssm_tbl.subj_idx(i));
                    bias_t = b_temp(j, hssm_tbl.subj_idx(i)) + b_stim100_temp(j, hssm_tbl.subj_idx(i));
                    ndt_t = n_temp(j, hssm_tbl.subj_idx(i)) + t_stim100_temp(j, hssm_tbl.subj_idx(i));
                elseif hssm_tbl.condition(i) == 3
                    alpha_t = a_temp(j, hssm_tbl.subj_idx(i)) + a_stim200_temp(j, hssm_tbl.subj_idx(i));
                    delta_t = v_temp(j, hssm_tbl.subj_idx(i)) + v_stim200_temp(j, hssm_tbl.subj_idx(i));
                    bias_t = b_temp(j, hssm_tbl.subj_idx(i)) + b_stim200_temp(j, hssm_tbl.subj_idx(i));
                    ndt_t = n_temp(j, hssm_tbl.subj_idx(i)) + t_stim200_temp(j, hssm_tbl.subj_idx(i));
                elseif hssm_tbl.condition(i) == 4
                    alpha_t = a_temp(j, hssm_tbl.subj_idx(i)) + a_stim300_temp(j, hssm_tbl.subj_idx(i));
                    delta_t = v_temp(j, hssm_tbl.subj_idx(i)) + v_stim300_temp(j, hssm_tbl.subj_idx(i));
                    bias_t = b_temp(j, hssm_tbl.subj_idx(i)) + b_stim300_temp(j, hssm_tbl.subj_idx(i));
                    ndt_t = n_temp(j, hssm_tbl.subj_idx(i)) + t_stim300_temp(j, hssm_tbl.subj_idx(i));
                elseif hssm_tbl.condition(i) == 5
                    alpha_t = a_temp(j, hssm_tbl.subj_idx(i)) + a_stim300_20_temp(j, hssm_tbl.subj_idx(i));
                    delta_t = v_temp(j, hssm_tbl.subj_idx(i)) + v_stim300_20_temp(j, hssm_tbl.subj_idx(i));
                    bias_t = b_temp(j, hssm_tbl.subj_idx(i)) + b_stim300_20_temp(j, hssm_tbl.subj_idx(i));
                    ndt_t = n_temp(j, hssm_tbl.subj_idx(i)) + t_stim300_20_temp(j, hssm_tbl.subj_idx(i));
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


%% E
subplot_tight(7,3,[14,17,20])
mid_params=[traces.a_stim100,traces.a_stim200,traces.a_stim300,traces.a_stim300_20];
mid_stds=[traces.a_std,traces.a_std,traces.a_std,traces.a_std];
mid_stds=sqrt((mid_stds.^2+[traces.a_stim100_std,traces.a_stim200_std,traces.a_stim300_std,traces.a_stim300_20_std].^2)/2);
normed = mid_params./mid_stds;
al_goodplot2({normed(:,1)',normed(:,2)',normed(:,3)',normed(:,4)'}, 'pos', [1,2,3,4], 'type', {'bilateral','bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[c1;c2;1,0,0;0.94,0.39,0.04], 'useMedian', true)
patch('XData',[0,0,6,6],'YData',effect_size*[-1,1,1,-1],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.6)
ylim([-2,3])
xticklabels(["100", "200", "300", "300@20"])
xlim([0.5,4.5])
xticks(1:4)
ylabel('Effect Size')
set(gca,'fontsize',18)
set(gca,'linewidth',2)
title('Boundary Separation')

for i=1:4
    disp(nnz(abs(normed(:,i))<effect_size)/length(normed(:,i))*100)
    disp(nnz(normed(:,i)>0)/length(normed(:,i))*100)
    disp(median(normed(:,i)))
    disp(' ')
end


%% F
subplot_tight(7,3,[15,18,21])
mid_params=[traces.v_stim100,traces.v_stim200,traces.v_stim300,traces.v_stim300_20];
mid_stds=[traces.v_std,traces.v_std,traces.v_std,traces.v_std];
mid_stds=sqrt((mid_stds.^2+[traces.v_stim100_std,traces.v_stim200_std,traces.v_stim300_std,traces.v_stim300_20_std].^2)/2);
normed = mid_params./mid_stds;
al_goodplot2({normed(:,1)',normed(:,2)',normed(:,3)',normed(:,4)'}, 'pos', [1,2,3,4], 'type', {'bilateral','bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[c1;c2;1,0,0;0.94,0.39,0.04], 'useMedian', true)
patch('XData',[0,0,6,6],'YData',effect_size*[-1,1,1,-1],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.6)
ylim([-2,3])
xticklabels(["100", "200", "300", "300@20"])
xlim([0.5,4.5])
xticks(1:4)
title('Drift Rate')
ylabel('Effect Size')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

for i=1:4
    disp(nnz(abs(normed(:,i))<effect_size)/length(normed(:,i))*100)
    disp(nnz(normed(:,i)>0)/length(normed(:,i))*100)
    disp(median(normed(:,i)))
    disp(' ')
end

%% Supplement Behavior
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% A
subplot_tight(2,2,1)
groups = {full_tbl.DT(strcmp(full_tbl.stim,'0')), ...
          full_tbl.DT(strcmp(full_tbl.stim,'1')), ...
          full_tbl.DT(strcmp(full_tbl.stim,'2')), ...
          full_tbl.DT(strcmp(full_tbl.stim,'3')), ...
          full_tbl.DT(strcmp(full_tbl.stim,'B'))};
c1 = ([0,0,1] * 2 + [1,0,0])/3;
c2 = ([0,0,1] + [1,0,0] * 2)/3;
al_goodplot2(groups', 'pos', [1,2,3,4,5], 'type', {'bilateral','bilateral','bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0,0,1;c1;c2;1,0,0;0.94,0.39,0.04],'useMedian', true)
scatter(ones(size(rat_tbl.RT_0)),rat_tbl.DT_0, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(2*ones(size(rat_tbl.RT_0)),rat_tbl.DT_1, 'MarkerFaceColor', c1, 'MarkerEdgeColor', c1)
scatter(3*ones(size(rat_tbl.RT_0)),rat_tbl.DT_2, 'MarkerFaceColor', c2, 'MarkerEdgeColor', c2)
scatter(4*ones(size(rat_tbl.RT_0)),rat_tbl.DT_3, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])
scatter(5*ones(size(rat_tbl.RT_0)),rat_tbl.DT_B, 'MarkerFaceColor', [0.94,0.39,0.04], 'MarkerEdgeColor', [0.94,0.39,0.04])

if coefTest(vp_dtglm, [0,0,1,0,0]) < 0.05
    plot([1,5], [24,24], 'k', 'LineWidth', 2)
    text(3, 25, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(vp_dtglm, [0,0,0,1,0]) < 0.05
    plot([1,4], [18,18], 'k', 'LineWidth', 2)
    text(2.5, 19, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(vp_dtglm, [0,0,0,0,1]) < 0.05
    plot([1,3], [14,14], 'k', 'LineWidth', 2)
    text(2, 15, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:5)
xticklabels(["OFF", "100", "200", "300", "300@20"])
ylabel('DT (s)')
set(gca, 'YScale', 'log')
ylim([0.4,32])
yticks([0.5,1,2,4,8,16,32])
xlim([0.5,5.5])
set(gca, 'fontsize', 18)
set(gca, 'linewidth', 2)
set(gca, 'YMinorTick', 'off')

%% B
subplot_tight(2,2,2)
groups = {ses_tbl.omission*100, ...
          ses_tbl.omission1*100, ...
          ses_tbl.omission2*100, ...
          ses_tbl.omission3*100, ...
          ses_tbl.omissionB*100};
c1 = ([0,0,1] * 2 + [1,0,0])/3;
c2 = ([0,0,1] + [1,0,0] * 2)/3;
al_goodplot2(groups', 'pos', [1,2,3,4,5], 'type', {'bilateral','bilateral','bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0,0,1;c1;c2;1,0,0;0.94,0.39,0.04])
scatter(ones(size(rat_tbl.RT_0)),rat_tbl.omission_0*100, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(2*ones(size(rat_tbl.RT_0)),rat_tbl.omission_1*100, 'MarkerFaceColor', c1, 'MarkerEdgeColor', c1)
scatter(3*ones(size(rat_tbl.RT_0)),rat_tbl.omission_2*100, 'MarkerFaceColor', c2, 'MarkerEdgeColor', c2)
scatter(4*ones(size(rat_tbl.RT_0)),rat_tbl.omission_3*100, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])
scatter(5*ones(size(rat_tbl.RT_0)),rat_tbl.omission_B*100, 'MarkerFaceColor', [0.94,0.39,0.04], 'MarkerEdgeColor', [0.94,0.39,0.04])

if coefTest(vp_omglm, [0,0,1,0,0]) < 0.05
    plot([1,5], [18,18], 'k', 'LineWidth', 2)
    text(3, 19, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(vp_omglm, [0,0,0,1,0]) < 0.05
    plot([1,4], [16,16], 'k', 'LineWidth', 2)
    text(2.5, 17, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(vp_omglm, [0,0,0,0,1]) < 0.05
    plot([1,3], [14,14], 'k', 'LineWidth', 2)
    text(2, 15, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:5)
xticklabels(["OFF", "100", "200", "300", "300@20"])
ylabel('% Omissions')
ylim([0,20])
xlim([0.5,5.5])
set(gca, 'fontsize', 18)
set(gca, 'linewidth', 2)

%% C
subplot_tight(2,2,3)
groups = {ses_tbl.frontChoice*100, ...
          ses_tbl.frontChoice1*100, ...
          ses_tbl.frontChoice2*100, ...
          ses_tbl.frontChoice3*100, ...
          ses_tbl.frontChoiceB*100};
c1 = ([0,0,1] * 2 + [1,0,0])/3;
c2 = ([0,0,1] + [1,0,0] * 2)/3;
al_goodplot2(groups', 'pos', [1,2,3,4,5], 'type', {'bilateral','bilateral','bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0,0,1;c1;c2;1,0,0;0.94,0.39,0.04])
scatter(ones(size(rat_tbl.RT_0)),rat_tbl.frontChoice_0*100, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(2*ones(size(rat_tbl.RT_0)),rat_tbl.frontChoice_1*100, 'MarkerFaceColor', c1, 'MarkerEdgeColor', c1)
scatter(3*ones(size(rat_tbl.RT_0)),rat_tbl.frontChoice_2*100, 'MarkerFaceColor', c2, 'MarkerEdgeColor', c2)
scatter(4*ones(size(rat_tbl.RT_0)),rat_tbl.frontChoice_3*100, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])
scatter(5*ones(size(rat_tbl.RT_0)),rat_tbl.frontChoice_B*100, 'MarkerFaceColor', [0.94,0.39,0.04], 'MarkerEdgeColor', [0.94,0.39,0.04])

if coefTest(vp_fcglm, [0,0,1,0,0]) < 0.05
    plot([1,5], [75,75], 'k', 'LineWidth', 2)
    text(3, 77, '*', 'HorizontalAlignment','center','fontsize', 30)
end
if coefTest(vp_fcglm, [0,0,0,1,0]) < 0.05
    plot([1,4], [70,70], 'k', 'LineWidth', 2)
    text(2.5, 72, '*', 'HorizontalAlignment','center','fontsize', 30)
end

xticks(1:5)
xticklabels(["OFF", "100", "200", "300", "300@20"])
ylabel('% Left Choice')
ylim([30,80])
yticks([30, 40,50,60,70,80])
xlim([0.5,5.5])
set(gca, 'fontsize', 18)
set(gca, 'linewidth', 2)

%% D
subplot_tight(2,2,4)
groups = {progress_tbl.progress(strcmp(progress_tbl.stim,"0")), ...
          progress_tbl.progress(strcmp(progress_tbl.stim,"1")), ...
          progress_tbl.progress(strcmp(progress_tbl.stim,"2")), ...
          progress_tbl.progress(strcmp(progress_tbl.stim,"3")), ...
          progress_tbl.progress(strcmp(progress_tbl.stim,"B"))};
c1 = ([0,0,1] * 2 + [1,0,0])/3;
c2 = ([0,0,1] + [1,0,0] * 2)/3;
al_goodplot2(groups', 'pos', [1,2,3,4,5], 'type', {'bilateral','bilateral','bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[0,0,1;c1;c2;1,0,0;0.94,0.39,0.04])
scatter(ones(size(rat_tbl.RT_0)),rat_tbl.progress_0, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1])
scatter(2*ones(size(rat_tbl.RT_0)),rat_tbl.progress_1, 'MarkerFaceColor', c1, 'MarkerEdgeColor', c1)
scatter(3*ones(size(rat_tbl.RT_0)),rat_tbl.progress_2, 'MarkerFaceColor', c2, 'MarkerEdgeColor', c2)
scatter(4*ones(size(rat_tbl.RT_0)),rat_tbl.progress_3, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])
scatter(5*ones(size(rat_tbl.RT_0)),rat_tbl.progress_B, 'MarkerFaceColor', [0.94,0.39,0.04], 'MarkerEdgeColor', [0.94,0.39,0.04])


xticks(1:5)
xticklabels(["OFF", "100", "200", "300", "300@20"])
ylabel('Rule Progress')
xlim([0.5,5.5])
set(gca, 'fontsize', 18)
set(gca, 'linewidth', 2)

%% Supplement block transitions
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

codes = ["0", "1", "2", "3", "B"];
linestyles = ["-.", "-", ":", "--"];
for i=1:4
    subplot_tight(2,2,i)
    hold on
    ind = 1;
    for j=1:length(codes)
        if i ~= j
            plot(ts,nanmean(blocks(strcmp(conditions(:,1), codes(j))&strcmp(conditions(:,2), codes(i)),:,1)),'Color',colors(j,:), 'LineWidth',2,'LineStyle', linestyles(ind))
            plot(ts2,nanmean(blocks(strcmp(conditions(:,1), codes(j))&strcmp(conditions(:,2), codes(i)),:,2)),'Color',colors(i,:), 'LineWidth',2,'LineStyle', linestyles(ind))
            ind = ind+1;
        end
    end
    ylim([0.4,1.2])
    ylabel('RT (s)')
    set(gca, 'fontsize', 18)
    set(gca,'linewidth',2)
    xlabel('Trials Since Change')
end

%% Supplement RLDDM
figure('Renderer', 'painters', 'Position', [100 1 1080 1080])
%ha = tight_subplot(2, 2, 0.05, 0.1, 0.05);
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

%% A-E
for j=1:5
    subplot_tight(2,5,j)
    hold on
    pts = 0:0.001:3;
    valid = hssm_tbl.condition==j;
    bw = 0.05;
    for s = 1:slices
        for i=1:50:T/slices
            temp = RT{s}(i,:);
            temp1 = squeeze(temp(valid));
            temp1=rmmissing(temp1(temp1>0));
            if ~isempty(temp1)
                vals = ksdensity(temp1,pts,'bandwidth',bw);
                plot(pts, vals, "Color",[colors(j,:), 0.05])
            end
        end
    end
    vals=ksdensity(squeeze(hssm_tbl.rt(valid&hssm_tbl.response>-1)),pts,'bandwidth',bw);
    plot(pts, vals, "Color",colors(j,:),'LineWidth',1.5)
    xticks(0:3)
    axis tight
    xlim([0,3])
    xlabel('RT (s)')
    yticks([])
    set(gca,'fontsize',18)
    set(gca,'linewidth',2)
end

%% F
subplot_tight(2,3,4)
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

%% G
subplot_tight(2,3,5)
mid_params=[traces.zt_stim100,traces.zt_stim200,traces.zt_stim300,traces.zt_stim300_20];
mid_stds=[traces.zt_std,traces.zt_std,traces.zt_std,traces.zt_std];
mid_stds=sqrt((mid_stds.^2+[traces.zt_stim100_std,traces.zt_stim200_std,traces.zt_stim300_std,traces.zt_stim300_20_std].^2)/2);
normed = mid_params./mid_stds;
al_goodplot2({normed(:,1)',normed(:,2)',normed(:,3)',normed(:,4)'}, 'pos', [1,2,3,4], 'type', {'bilateral','bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[c1;c2;1,0,0;0.94,0.39,0.04], 'useMedian', true)
patch('XData',[0,0,6,6],'YData',effect_size*[-1,1,1,-1],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.6)
ylim([-2,2])
xticklabels(["100", "200", "300", "300@20"])
xlim([0.5,4.5])
xticks(1:4)
ylabel('Effect Size')
set(gca,'fontsize',18)
set(gca,'linewidth',2)
title('Bias')

%% H
subplot_tight(2,3,6)
mid_params=[traces.t_stim100,traces.t_stim200,traces.t_stim300,traces.t_stim300_20];
mid_stds=[traces.t_std,traces.t_std,traces.t_std,traces.t_std];
mid_stds=sqrt((mid_stds.^2+[traces.t_stim100_std,traces.t_stim200_std,traces.t_stim300_std,traces.t_stim300_20_std].^2)/2);
normed = mid_params./mid_stds;
al_goodplot2({normed(:,1)',normed(:,2)',normed(:,3)',normed(:,4)'}, 'pos', [1,2,3,4], 'type', {'bilateral','bilateral','bilateral','bilateral'},'boxw',0.4, 'col',[c1;c2;1,0,0;0.94,0.39,0.04], 'useMedian', true)
patch('XData',[0,0,6,6],'YData',effect_size*[-1,1,1,-1],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.6)
ylim([-2,2])
xticklabels(["100", "200", "300", "300@20"])
xlim([0.5,4.5])
xticks(1:4)
ylabel('Effect Size')
set(gca,'fontsize',18)
set(gca,'linewidth',2)
title('Non-decision Time')