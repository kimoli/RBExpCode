close all
clear all

machine = 'COMPUPITAR';

if strcmpi(machine, 'OREK')
    basedir = 'E:\pcp2ChR2 data\rebound';
elseif strcmpi(machine, 'COMPUPITAR')
    basedir = 'D:\pcp2ChR2 data\rebound';
else
    disp('Please specify computer so we know what directory to use')
    basedir = '';
end
cd(basedir)
load('overallData.mat')

mice = unique(pretestData.mouse);

% code phase in number instead of in string
%       phase == 1 == pretest
%       phase == 2 == post training
%       phase == 3 == post extinction
rbstats.mouse = [];
rbstats.phase = [];
rbstats.amp = [];
rbstats.lat = [];
rbstats.lasamp = [];
rbstats.laspow = [];
rbstats.winstart = [];

crstats.mouse = [];
crstats.phase = [];
crstats.amp = [];

daycrstats.mouse = [];
daycrstats.phase = [];
daycrstats.amp = [];
daycrstats.prob = [];
daycrstats.meantr = [];

daystats.mouse = [];
daystats.phase = [];
daystats.rb.amp = [];
daystats.rb.lat = [];
daystats.rb.prob = [];
daystats.lasamp = [];
daystats.laspow = [];
daystats.meanRBTr = [];
daystats.winstart = [];


timeVector = 1:size(posttestData.eyelidpos,2);
timeVector = timeVector * 0.00488372;
timeVector = timeVector - 0.2;


for m = 1:length(mice)
    [rbstats, crstats, daystats, daycrstats] = getrbprops(pretestData,1, rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(lastacqData,1.5, rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(posttestData,2, rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(lastextData,2.5, rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(postextData,3, rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    %plot(daystats.meanRBTr')
end


%% plot the eyelid traces
% plot for lowest intensity
pretestidx = daystats.phase == 1 & daystats.laspow==15;
posttestidx = daystats.phase == 2  & daystats.laspow==15;
postextidx = daystats.phase == 3  & daystats.laspow==15;

figure
subplot(4,3,1)
plot(timeVector, daystats.meanRBTr(pretestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
ylabel('15 mW')
title('pre training')
xlim([0.85 max(timeVector)])

subplot(4,3,2)
plot(timeVector, daystats.meanRBTr(posttestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
title('post training')
xlim([0.85 max(timeVector)])

subplot(4,3,3)
plot(timeVector, daystats.meanRBTr(postextidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
title('post extinction')
xlim([0.85 max(timeVector)])

% plot for middle intensity
pretestidx = daystats.phase == 1 & daystats.laspow==30;
posttestidx = daystats.phase == 2  & daystats.laspow==30;
postextidx = daystats.phase == 3  & daystats.laspow==30;

subplot(4,3,4)
plot(timeVector, daystats.meanRBTr(pretestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([0.85 max(timeVector)])
ylabel('30 mW')

subplot(4,3,5)
plot(timeVector, daystats.meanRBTr(posttestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([0.85 max(timeVector)])

subplot(4,3,6)
plot(timeVector, daystats.meanRBTr(postextidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([0.85 max(timeVector)])

% plot for highest intensity
pretestidx = daystats.phase == 1 & daystats.laspow==60;
posttestidx = daystats.phase == 2  & daystats.laspow==60;
postextidx = daystats.phase == 3  & daystats.laspow==60;

subplot(4,3,7)
plot(timeVector, daystats.meanRBTr(pretestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([0.85 max(timeVector)])
ylabel('60 mW')

subplot(4,3,8)
plot(timeVector, daystats.meanRBTr(posttestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([0.85 max(timeVector)])

subplot(4,3,9)
plot(timeVector, daystats.meanRBTr(postextidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([0.85 max(timeVector)])

% plot crs
pretestidx = daycrstats.phase == 1;
posttestidx = daycrstats.phase == 1.5;
postextidx = daycrstats.phase == 2.5;

subplot(4,3,10)
plot(timeVector, daycrstats.meantr(pretestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daycrstats.meantr(pretestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([0 0.6105])
ylabel('CS')

subplot(4,3,11)
plot(timeVector, daycrstats.meantr(posttestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daycrstats.meantr(posttestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([0 0.6105])

subplot(4,3,12)
plot(timeVector, daycrstats.meantr(postextidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daycrstats.meantr(postextidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([0 0.6105])

%% plot relationship between CR and RB amplitudes across mice
rbidx = daystats.phase == 2  & daystats.laspow==60;
cridx = daycrstats.phase == 1.5;
cridxnums = find(cridx);
figure
scatter(daystats.rb.amp(rbidx), daycrstats.amp(cridxnums(5:end)))
xlim([0 1])
ylim([0 1])
% lol it is so flat

%% plot relationship between CR and RB probabilities across mice
rbidx = daystats.phase == 2  & daystats.laspow==60;
cridx = daycrstats.phase == 1.5;
cridxnums = find(cridx);
figure
scatter(daystats.rb.prob(rbidx), daycrstats.prob(cridxnums(5:end)))
xlim([0 1])
ylim([0 1])
% lol it is so flat

%% stats
% RB before vs after training
preidx = daystats.phase == 1 & daystats.laspow==60;
postidx = daystats.phase == 2  & daystats.laspow==60;
