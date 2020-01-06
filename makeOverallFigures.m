close all
clear all

machine = 'OREK';

if strcmpi(machine, 'OREK')
    basedir = 'E:\pcp2ChR2 data\rebound';
elseif strcmpi(machine, 'COMPUPITAR')
    basedir = 'D:\pcp2ChR2 data\rebound';
else
    disp('Please specify computer so we know what directory to use')
    basedir = '';
end
cd(basedir)
load('overallData_191219.mat')

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
daystats.rb.hitamp = [];
daystats.rb.lat = [];
daystats.rb.prob = [];
daystats.lasamp = [];
daystats.laspow = [];
daystats.meanRBTr = [];
daystats.meanRBTrHit = [];
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

%% compare animals performance within group across light intensities
[signrankp_paired, ttestp_paired] = makeSubplotSummary(daystats, 211, 218, 'Paired');
[signrankp_unpaired, ttestp_unpaired] = makeSubplotSummary(daystats, 234, 241, 'Unpaired');

[prevals, postvals, premad, postmad]=makePlots_RBxPower(daystats, daystats.rb.amp, [15,30,60], [211,213,214,215,216,217,218], 1, 1, 2);
xlabel('Laser Power (mW)')
ylabel('Rebound Amplitude (FEC)')
ylim([0 1])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')
checkShapiroWilk = prevals(:,1)-postvals(:,1); % R shapiro.test: W = 0.85763, p-value = 0.1442
checkShapiroWilk = prevals(:,2)-postvals(:,2); % R shapiro.test: W = 0.86774, p-value = 0.1773
checkShapiroWilk = prevals(:,3)-postvals(:,3); % R shapiro.test: W = 0.91333, p-value = 0.4194
[h, p, ci, stats] = ttest(prevals(:,1),postvals(:,1)) % p = 0.0461
[h, p, ci, stats] = ttest(prevals(:,2),postvals(:,2)) % p = 0.0135
[h, p, ci, stats] = ttest(prevals(:,3),postvals(:,3)) % p = 0.0198

[pretestprobs, trainprobs, pretestprobmad, trainprobmad]=makePlots_RBxPower(daystats, daystats.rb.prob, [15,30,60], [211,213,214,215,216,217,218], 1, 1, 2);
xlabel('Laser Power (mW)')
ylabel('Rebound Probability (FEC)')
ylim([0 1])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')
    writeTestDiffs([pretestprobs(:,1),trainprobs(:,1)]) % R Shapiro.test: W = 0.66444, p-value = 0.001497
    [p,h,stats] = signrank(pretestprobs(:,1),trainprobs(:,1)) % p = 0.25
    writeTestDiffs([pretestprobs(:,2),trainprobs(:,2)]) % R Shapiro.test: W = 0.75174, p-value = 0.01323
    [p,h,stats] = signrank(pretestprobs(:,2),trainprobs(:,2)) % p = 0.0156
    writeTestDiffs([pretestprobs(:,3),trainprobs(:,3)]) % R Shapiro.test: W = 0.91532, p-value = 0.4339
    [h,p,ci,stats] = ttest(pretestprobs(:,2),trainprobs(:,2)) % p = 0.0054

checkShapiroWilk = median(prevals,2)-median(postvals,2); % R shapiro.test: W = 0.90486, p-value = 0.3614
[h,p,ci,stats] = ttest(median(prevals,2),median(postvals,2),'Tail','left');
% plot boxplot
predata = median(prevals,2);
postdata = median(postvals,2);
prequant = quantile(predata,3);
postquant = quantile(postdata,3);
figure
plot([0.75 1.25], [prequant(2) prequant(2)], 'Color', [0 0 1])
hold on
plot([0.75, 1.25], [prequant(1) prequant(1)], 'Color', [0 1 1])
plot([0.75, 1.25], [prequant(3) prequant(3)], 'Color', [0 1 1])
plot([1.75 2.25], [postquant(2) postquant(2)], 'Color', [1 0 0])
plot([1.75 2.25], [postquant(1) postquant(1)], 'Color', [1 0 1])
plot([1.75 2.25], [postquant(3) postquant(3)], 'Color', [1 0 1])
xlim([0.5 2.5])
ylim([0 1])
for i = 1:7
    plot([1 2],[predata(i) postdata(i)], 'LineStyle', ':', 'Color', [0 0 0])
end
scatter(ones(7,1), predata, 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[0 0 1])
scatter(ones(7,1)*2, postdata, 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[1 0 0])
ylabel('RB probability')

[preamps, postTrainAmps, premad, postmad]=makePlots_RBxPower(daystats, daystats.rb.hitamp, [15,30,60], [211,213,214,215,216,217,218], 1, 1, 2);
xlabel('Laser Power (mW)')
ylabel('Rebound Amplitude (FEC)')
ylim([0 1])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')

figure
quants30 = quantile(postTrainAmps(:,2),3);
quants60 = quantile(postTrainAmps(:,3),3);
plot([0.75 1.25], [quants30(2) quants30(2)], 'Color', [0 0 1])
hold on
plot([0.75 1.25], [quants30(1) quants30(1)], 'Color', [0 1 1])
plot([0.75 1.25], [quants30(3) quants30(3)], 'Color', [0 1 1])
plot([1.75 2.25], [quants60(2) quants60(2)], 'Color', [1 0 0])
plot([1.75 2.25], [quants60(1) quants60(1)], 'Color', [1 0 1])
plot([1.75 2.25], [quants60(3) quants60(3)], 'Color', [1 0 1])
for i = 1:7
    plot([1 2],[postvals(i,2) postvals(i,3)], 'LineStyle', ':', 'Color', [0 0 0])
end
scatter(ones(7,1), postvals(:,2), 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[0 0 1])
scatter(ones(7,1)*2, postvals(:,3), 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[1 0 0])
ylim([0 1])
xlim([0.5 2.5])
xlabel('Laser Power (mW)')
ylabel('Rebound Amplitude (FEC)')

% for experimental group extinction: looking at rebound probabilities
[trainprobs, extprobs, trainprobmad, extprobmad]=makePlots_RBxPower(daystats, daystats.rb.prob, [15,30,60], [211,213,214,215,216,217,218], 1, 2, 3);
xlabel('Laser Power (mW)')
ylabel('Rebound Probability')
ylim([0 1])
xlim([0.5 3.5])
[pretestprobs, ~, pretestprobmad, ~]=makePlots_RBxPower(daystats, daystats.rb.prob, [15,30,60], [211,213,214,215,216,217,218], 0, 1, nan);
legend('Post-Training','Post-Training','Post-Extinction','Post-Extinction','Pre-Test','Pre-Test', 'Location', 'NorthWest')

    % check stats - pretest vs postext there is no statistically significant difference between
    % RB probabilities during pretest and extinction when we do bonferroni
    % correction, might be able to argue that you wouldn't bother doing the
    % comparison between values where the mean is so obviously 0
    checkShapiroWilk = pretestprobs(:,1)-extprobs(:,1); % R shapiro.test: uh there is no difference here
    [p,h,stats] = signrank(pretestprobs(:,1),extprobs(:,1)); % p = 1
    checkShapiroWilk = pretestprobs(:,2)-extprobs(:,2); % shapiro.test: W = 0.58193, p-value = 0.0001664
    [p,h,stats] = signrank(pretestprobs(:,2),extprobs(:,2)); % p = 0.5
    checkShapiroWilk = pretestprobs(:,3)-extprobs(:,3); % shapiro.test: W = 0.86235, p-value = 0.1974
    [h,p,ci,stats]=ttest(pretestprobs(:,3),extprobs(:,3)) % p = 0.0772
        
     % check stats - traintest vs postext (diference at 30 mw but not 15
     % and 60)
     writeTestDiffs([trainprobs(:,1),extprobs(:,1)]) % R shapiro.test: W = 0.81244, p-value = 0.05426
     [h,p,ci,stats] = ttest(trainprobs(:,1),extprobs(:,1)) % p = 0.078
     writeTestDiffs([trainprobs(:,2), extprobs(:,2)])  % R shapiro.test: W = 0.86828, p-value = 0.1793
     [h,p,ci,stats] = ttest(trainprobs(:,2),extprobs(:,2)) % p = 0.0054
     writeTestDiffs([trainprobs(:,3),extprobs(:,3)]) % R shapiro.test: W = 0.91649, p-value = 0.4426
     [h,p,ci,stats] = ttest(trainprobs(:,3),extprobs(:,3)) % p = 0.1082
     
% for experimental group extinction: looking at rebound amplitudes
[trainamps, extamps, trainampmad, extampmad]=makePlots_RBxPower(daystats, daystats.rb.amp, [15,30,60], [211,213,214,215,216,217,218], 1, 2, 3);
xlabel('Laser Power (mW)')
ylabel('Rebound Amplitude (FEC)')
ylim([0 1])
xlim([0.5 3.5])
[pretestamps, ~, pretestampmad, ~]=makePlots_RBxPower(daystats, daystats.rb.amp, [15,30,60], [211,213,214,215,216,217,218], 0, 1, nan);
legend('Post-Training','Post-Training','Post-Extinction','Post-Extinction','Pre-Test','Pre-Test', 'Location', 'NorthWest')

    % check stats - pretest vs postext there is no statistically significant difference between
    % RB probabilities during pretest and extinction when we do bonferroni
    % correction, might be able to argue that you wouldn't bother doing the
    % comparison between values where the mean is so obviously 0
    writeTestDiffs([pretestamps(:,1),extamps(:,1)]) % R shapiro.test: W = 0.94186, p-value = 0.6555
    [h,p,ci,stats] = ttest(pretestamps(:,1),extamps(:,1)) % p = 0.0701
    writeTestDiffs([pretestamps(:,2),extamps(:,2)]) % R Shapiro.test: W = 0.61637, p-value = 0.000423
    [p,h,stats] = signrank(pretestamps(:,2),extamps(:,2)) % p = 0.0469
    writeTestDiffs([pretestamps(:,3),extamps(:,3)]) % R SHapiro.test: W = 0.88522, p-value = 0.2506
    [h,p,ci,stats]=ttest(pretestamps(:,3),extamps(:,3), 'Tail', 'left') % p = 0.0484, too low for this set of comparisons
        % might be able to get away with just reporting this one (at least for
        % the thesis, will see what Javier thinks)
        
     % check stats - traintest vs postext (diference at 30 and 60 mW)
     writeTestDiffs([trainamps(:,1),extamps(:,1)]) % R shapiro.test: W = 0.88348, p-value = 0.2424
     [h,p,ci,stats] = ttest(trainamps(:,1),extamps(:,1)) % p = 0.0821, this may be a nonsense comparison because we know there are generally no rebounds at 15 mW so there isn't going to be much data here
     writeTestDiffs([trainamps(:,2), extamps(:,2)])  % R shapiro.test: W = 0.91532, p-value = 0.434
     [h,p,ci,stats] = ttest(trainamps(:,2),extamps(:,2)) % p = 0.0054
     writeTestDiffs([trainamps(:,3),extamps(:,3)]) % R shapiro.test: W = 0.89872, p-value = 0.3233
     [h,p,ci,stats] = ttest(trainamps(:,3),extamps(:,3)) % p = 0.0147
     
% make a CSV for doing a repeated measures ANOVA on rebound probabilities
% in R
% let col1 be experiment phase, col2 be mouse, col3 be laser power, and col
% 4 be RB probability
theseMice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218'};
powers = [15;30;60];
bigcell = cell(63,4); % 7*3*3 rows (n*phases*powers), 4 columns
iter = 1;
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'pretest';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = pretestprobs(i,c);
        iter = iter+1;
    end
end
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'post_training';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = trainprobs(i,c);
        iter = iter+1;
    end
end
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'post_extinction';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = extprobs(i,c);
        iter = iter+1;
    end
end
bigtable = cell2table(bigcell,'VariableNames',{'phase','mouse','power_mW','RB_prob'});
writetable(bigtable, 'RBProbs_20200104.csv')

% make a CSV for doing a repeated measures ANOVA on rebound probabilities
% in R
% let col1 be experiment phase, col2 be mouse, col3 be laser power, and col
% 4 be RB probability
theseMice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218'};
powers = [15;30;60];
bigcell = cell(63,4); % 7*3*3 rows (n*phases*powers), 4 columns
iter = 1;
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'pretest';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = pretestamps(i,c);
        iter = iter+1;
    end
end
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'post_training';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = trainamps(i,c);
        iter = iter+1;
    end
end
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'post_extinction';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = extamps(i,c);
        iter = iter+1;
    end
end
bigtable = cell2table(bigcell,'VariableNames',{'phase','mouse','power_mW','RB_amp'});
writetable(bigtable, 'RBAmps_20200104.csv')


% for control mice
[preunpprob, postunpprob, preunpmad, postunpmad]=makePlots_RBxPower(daystats, daystats.rb.prob, [15,30,60], [234,235,236,237,238,239,240,241], 1, 1, 2);
xlabel('Laser Power (mW)')
ylabel('Rebound Probability (FEC)')
ylim([0 1])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')
% checkShapiroWilk = prevals(:,1)-postvals(:,1); % R shapiro.test: W = 0.66444, p-value = 0.001497
% checkShapiroWilk = prevals(:,2)-postvals(:,2); % R shapiro.test: W = 0.75174, p-value = 0.01323
% checkShapiroWilk = prevals(:,3)-postvals(:,3); % R shapiro.test: W = 0.91532, p-value = 0.4339
% [p,h,stats] = signrank(prevals(:,1),postvals(:,1));
% [p,h,stats] = signrank(prevals(:,2),postvals(:,2));
% [h, p, ci, stats] = ttest(prevals(:,3),postvals(:,3));
% 

% compare pre and post regardless of laser power, for control mice
checkShapiroWilk = median(preunpprob,2)-median(postunpprob,2); % R shapiro.test: W = 0.81042, p-value = 0.03697
[p,h,stats] = signrank(median(preunpprob,2),median(postunpprob,2));
% plot boxplot
predata = median(preunpprob,2);
postdata = median(postunpprob,2);
prequant = quantile(predata,3);
postquant = quantile(postdata,3);
figure
plot([0.75 1.25], [prequant(2) prequant(2)], 'Color', [0 0 1])
hold on
plot([0.75, 1.25], [prequant(1) prequant(1)], 'Color', [0 1 1])
plot([0.75, 1.25], [prequant(3) prequant(3)], 'Color', [0 1 1])
plot([1.75 2.25], [postquant(2) postquant(2)], 'Color', [1 0 0])
plot([1.75 2.25], [postquant(1) postquant(1)], 'Color', [1 0 1])
plot([1.75 2.25], [postquant(3) postquant(3)], 'Color', [1 0 1])
xlim([0.5 2.5])
ylim([0 1])
for i = 1:8
    plot([1 2],[predata(i) postdata(i)], 'LineStyle', ':', 'Color', [0 0 0])
end
scatter(ones(8,1), predata, 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[0 0 1])
scatter(ones(8,1)*2, postdata, 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[1 0 0])
ylabel('RB probability')

% save a csv for comparing rebounds after training in the experimental and
% the control animals
theseMice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218'};
powers = [15;30;60];
bigcell = cell(42,4); % 7*2*3 rows (n*phases*powers), 4 columns
iter = 1;
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'post_training';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = trainprobs(i,c);
        iter = iter+1;
    end
end
theseMice = {'OK234';'OK235';'OK236';'OK237';'OK238';'OK239';'OK241';'OK242'};
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'post_unpaired';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = postunpprob(i,c);
        iter = iter+1;
    end
end
bigtable = cell2table(bigcell,'VariableNames',{'phase','mouse','power_mW','RB_prob'});
writetable(bigtable, 'RBProbs_expCont_20200104.csv')


%% plot the eyelid traces centered on event of interest
makeEyetraceSubplots(daystats, daycrstats, 211, 218, 'paired',timeVector)
makeEyetraceSubplots(daystats, daycrstats, 234, 241, 'unpaired',timeVector)
figure
makeEyetraceSubplots_collapsed(daystats, daycrstats, 211, 218, 'paired',timeVector)
figure
makeEyetraceSubplots_collapsed(daystats, daycrstats, 234, 241, 'unpaired',timeVector)


%% plot the eyelid traces centered on event of interest
a = find(daycrstats.phase==1.5);
b = find(daycrstats.mouse>=211);
c = find(daycrstats.mouse<=218);
temp = intersect(a,b);
posttestidx = intersect(temp,c);
posttestrbidx60 = find(daystats.phase == 2  & daystats.laspow==60& daystats.mouse>=211 & daystats.mouse<=218);
posttestrbidx30 = find(daystats.phase == 2  & daystats.laspow==30& daystats.mouse>=211 & daystats.mouse<=218);
posttestrbidx15 = find(daystats.phase == 2  & daystats.laspow==15& daystats.mouse>=211 & daystats.mouse<=218);

colordef black
figure
for i = 1:7
    subplot(2,4,i)
    plot(timeVector, daycrstats.meantr(posttestidx(i),:)-mean(daycrstats.meantr(posttestidx(i),1:40)), 'Color', [1 1 1])
    hold on
    plot(timeVector-0.85, daystats.meanRBTr(posttestrbidx60(i),:)-mean(daystats.meanRBTr(posttestrbidx60(i),1:40)), 'Color', [0 1 0])
    plot(timeVector-0.85, daystats.meanRBTr(posttestrbidx30(i),:)-mean(daystats.meanRBTr(posttestrbidx30(i),1:40)), 'Color', [1 0 1])
    plot(timeVector-0.85, daystats.meanRBTr(posttestrbidx15(i),:)-mean(daystats.meanRBTr(posttestrbidx15(i),1:40)), 'Color', [0 1 1])
    ylim([-0.05 1])
    xlim([-0.1 0.6])
end

a = find(daycrstats.phase==1.5);
b = find(daycrstats.mouse>=234);
c = find(daycrstats.mouse<=241);
temp = intersect(a,b);
posttestidx = intersect(temp,c);
posttestrbidx60 = find(daystats.phase == 2  & daystats.laspow==60& daystats.mouse>=234 & daystats.mouse<=241);
posttestrbidx30 = find(daystats.phase == 2  & daystats.laspow==30& daystats.mouse>=234 & daystats.mouse<=241);
posttestrbidx15 = find(daystats.phase == 2  & daystats.laspow==15& daystats.mouse>=234 & daystats.mouse<=241);

figure
for i = 1:8
    subplot(2,4,i)
    plot(timeVector, daycrstats.meantr(posttestidx(i),:)-mean(daycrstats.meantr(posttestidx(i),1:40)), 'Color', [1 1 1])
    hold on
    plot(timeVector-0.85, daystats.meanRBTr(posttestrbidx60(i),:)-mean(daystats.meanRBTr(posttestrbidx60(i),1:40)), 'Color', [0 1 0])
    plot(timeVector-0.85, daystats.meanRBTr(posttestrbidx30(i),:)-mean(daystats.meanRBTr(posttestrbidx30(i),1:40)), 'Color', [1 0 1])
    plot(timeVector-0.85, daystats.meanRBTr(posttestrbidx15(i),:)-mean(daystats.meanRBTr(posttestrbidx15(i),1:40)), 'Color', [0 1 1])
    ylim([-0.05 1])
    xlim([-0.1 0.6])
end

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

%% check if there is a relationship between days to reacquire and the rebounds after extinction
basedir = 'E:\pcp2ChR2 data\rebound';
mice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218'};

daily.mouse = {};
daily.crprob = [];
daily.cradjamp = [];
daily.cradjampHit = [];
daily.sesstype = []; % 0 = training, 1 = extinction, 2 = laser, 3 = laser extinction trial
for m = 1:length(mice)
    mousedir = [basedir,'\',mice{m,1}];
    cd(mousedir)
    
    days = dir('19*');
    for d = 1:length(days)
        
        daydir = [mousedir,'\',days(d,1).name];
        cd(daydir)
        
        % check if there is a trialdata file present
        loaded = false;
        if exist('newTrialdata.mat','file')==2
            load('newTrialdata.mat')
            loaded = true;
        elseif exist('trialdata.mat','file')==2
            load('trialdata.mat')
            loaded = true;
        end
        
        if loaded == true
            % check if this day had CS trials or if it was a laser test day
            numCSTrials = sum(trials.c_csdur>0);
            
            if numCSTrials > 0
                CSUSidx = find(trials.c_csdur>0 & trials.c_usdur>0);
                CSOnlyidx = find(trials.c_csdur>0 & trials.c_usdur==0);
                laserIdx = find(trials.laser.dur>0);               
                if ~isempty(CSUSidx)
                    checkidx = CSUSidx;
                    if length(laserIdx)>length(CSUSidx)                        
                        daily.sesstype(end+1,1) = 3;
                    else
                        daily.sesstype(end+1,1) = 0;
                    end
                else
                    checkidx = CSOnlyidx;
                    daily.sesstype(end+1,1) = 1;
                end
                
                baseline = nan(length(checkidx),1);
                cradjamp = nan(length(checkidx),1);
                stable = nan(length(checkidx),1);
                %eyeadj = nan(length(checkidx),size(trials.eyelidpos,2));
                for t = 1:length(checkidx)
                    baseline(t,1) = mean(trials.eyelidpos(checkidx(t,1),1:40));
                    stable(t,1) = max(trials.eyelidpos(checkidx(t,1),1:40))<0.3;
                    cradjamp(t,1) = max(trials.eyelidpos(checkidx(t,1),76:85))-baseline(t,1);
                    %eyeadj(t,:)=trials.eyelidpos(checkidx(t,1),:)-baseline(t,1);
                end
                
                meancradjamp = median(cradjamp(stable==1,1));
                meancradjampHit = median(cradjamp(stable==1 & cradjamp>=0.1));
                crprob = sum(cradjamp(stable==1,1)>=0.1)./sum(stable==1);
                
%                 figure
%                 plot(eyeadj(stable==1,:)')
%                 hold on
%                 plot(mean(eyeadj(stable==1,:)),'Color',[0 0 0],'LineWidth',3)
%                 title(num2str(meancradjamp))
%                 pause
                
                daily.mouse{end+1,1} = mice{m,1};
                daily.cradjamp(end+1,1) = meancradjamp;
                daily.cradjampHit(end+1,1) = meancradjampHit;
                daily.crprob(end+1,1) = crprob;
                
                clear trials meancradjamp crprob baseline cradjamp stable CSUSidx...
                    checkidx CSOnlyidx
            else
                daily.mouse{end+1,1} = mice{m,1};
                daily.sesstype(end+1,1) = 2;
                daily.cradjamp(end+1,1) = nan;
                daily.crprob(end+1,1) = nan;
                daily.cradjampHit(end+1,1) = nan;
            end
            
            clear numCSTrials
            
        end
        
        clear loaded daydir
    end
    
    clear mousedir
end


% for m = 1:length(mice)
%     thisMouse = mice(m,1);
%     idx = [];
%     for i = 1:length(daily.mouse)
%         if strcmpi(thisMouse, daily.mouse{i,1})
%             idx(end+1)=i;
%         end
%     end
%     viewData = [daily.sesstype(idx), daily.crprob(idx), daily.cradjamp(idx), daily.cradjampHit(idx)];
%     pause
%     clear idx
% end

% let performance threshold be >80% CRs and >50% median amp of hit trials
threshProb = 0.8;
threshAmp = 0.6; %0.5 if use col 4, this is more accurate for asymptotes though
numdays.mouse = mice;
numdays.acq = nan(7,1);
numdays.reacq = nan(7,1);
for m = 1:length(mice)
    thisMouse = mice(m,1);
    idx = [];
    for i = 1:length(daily.mouse)
        if strcmpi(thisMouse, daily.mouse{i,1})
            idx(end+1)=i;
        end
    end
    
    tempArr = [daily.sesstype(idx), daily.crprob(idx), daily.cradjamp(idx), ...
        daily.cradjampHit(idx)];
    
    % get days to acquisition
    counter = 1;
    iteridx = 2;
    while tempArr(iteridx,1)==0
        if tempArr(iteridx,2)>=threshProb && tempArr(iteridx,3)>=threshAmp
            break
        else
            counter = counter+1;
        end
        iteridx = iteridx + 1;
        if iteridx > length(tempArr) || ~tempArr(iteridx,1)==0
            counter = nan;
            break
        end
    end
    numdays.acq(m,1) = counter;
    

    % get days to reacquisition
    trainidx = find(tempArr(:,1)==2);
    iteridx = trainidx(3)+1;
    counter = 1;    
    while tempArr(iteridx,1)==0
        if tempArr(iteridx,2)>=threshProb && tempArr(iteridx,3)>=threshAmp
            break
        else
            counter = counter+1;
        end
        iteridx = iteridx + 1;
        if iteridx > length(tempArr) || ~tempArr(iteridx,1)==0
            counter = nan;
            break
        end
    end
    numdays.reacq(m,1) = counter;
    
end

savingsRatio = numdays.reacq./numdays.acq;

%% plot relationship between savings and rb after extinction
rbidx = daystats.phase == 3  & daystats.laspow==60 & daystats.mouse>=211 & daystats.mouse<=218;
figure
scatter(daystats.rb.prob(rbidx), savingsRatio)
xlabel('RB probability after extinction')
ylabel('days to reacq/days to acq')
ylim([0 0.6])
xlim([0 1])
temprb = daystats.rb.prob(rbidx);
[r, p] = corr(temprb(~isnan(savingsRatio)), savingsRatio(~isnan(savingsRatio)));
lsline
titlestring = ['Pearsons r = ', num2str(r), '; p = ', num2str(p)];
title(titlestring)

figure
scatter(daystats.rb.amp(rbidx), savingsRatio)
xlabel('RB amplitude after extinction')
ylabel('days to reacq/days to acq')
ylim([0 0.6])
xlim([0 0.5])
temprb = daystats.rb.amp(rbidx);
[r, p] = corr(temprb(~isnan(savingsRatio)), savingsRatio(~isnan(savingsRatio)));
lsline
titlestring = ['Pearsons r = ', num2str(r), '; p = ', num2str(p)];
title(titlestring)

% plot experimental group animals' performance across training
CRProbsXDays = [];
CRAmpsXDays = [];
for m = 1:length(mice)
    thisMouse = mice(m,1);
    idx = [];
    for i = 1:length(daily.mouse)
        if strcmpi(thisMouse, daily.mouse{i,1})
            idx(end+1)=i;
        end
    end
    temp = find(daily.sesstype==2);
    temp2 = ismember(idx,temp);
    testidx = idx(temp2);
    pretestidx = testidx(1);
    posttestidx = testidx(2);
    crprobvals = daily.crprob(pretestidx+1:posttestidx-1,1)';
    crampvals = daily.cradjampHit(pretestidx+1:posttestidx-1,1)';
    if isempty(CRProbsXDays)
        CRProbsXDays = crprobvals;
        CRAmpsXDays = crampvals;
    elseif size(CRProbsXDays,2) > size(crprobvals,2)
        while size(CRProbsXDays,2) > size(crprobvals,2)
            crprobvals = [crprobvals,NaN];
            crampvals = [crampvals,NaN];
        end
        CRProbsXDays = [CRProbsXDays; crprobvals];
        CRAmpsXDays = [CRAmpsXDays; crampvals];
    elseif size(CRProbsXDays,2)<size(crprobvals,2)
        while size(CRProbsXDays,2)<size(crprobvals,2)
            CRProbsXDays = [CRProbsXDays,nan(size(CRProbsXDays,1),1)];
            CRAmpsXDays = [CRAmpsXDays,nan(size(CRAmpsXDays,1),1)];
        end
        CRProbsXDays = [CRProbsXDays; crprobvals];
        CRAmpsXDays = [CRAmpsXDays; crampvals];
    else
        CRProbsXDays = [CRProbsXDays; crprobvals];
        CRAmpsXDays = [CRAmpsXDays; crampvals];
    end
end

colordef white
figure
errorbar([1:19], nanmedian(CRProbsXDays), mad(CRProbsXDays,1), '.', 'LineStyle', 'none', 'Color', [0 0 0])
hold on
errorbar([1:19], nanmedian(CRAmpsXDays), mad(CRAmpsXDays,1), '.', 'LineStyle', 'none', 'Color', [1 0 1])
xlabel('Sessions')
ylabel('CR Probability')

% plot experimental group animals' performance across extinction
CRProbsXDays = [];
CRAmpsXDays = [];
for m = 1:length(mice)
    thisMouse = mice(m,1);
    idx = [];
    for i = 1:length(daily.mouse)
        if strcmpi(thisMouse, daily.mouse{i,1})
            idx(end+1)=i;
        end
    end
    temp = find(daily.sesstype==1);
    temp2 = ismember(idx,temp);
    testidx = idx(temp2);
    crprobvals = daily.crprob(testidx(1)-2:testidx(1)+9,1)';
    crampvals = daily.cradjampHit(testidx(1)-2:testidx(1)+9,1)';
    if isempty(CRProbsXDays)
        CRProbsXDays = crprobvals;
        CRAmpsXDays = crampvals;
    elseif size(CRProbsXDays,2) > size(crprobvals,2)
        while size(CRProbsXDays,2) > size(crprobvals,2)
            crprobvals = [crprobvals,NaN];
            crampvals = [crampvals,NaN];
        end
        CRProbsXDays = [CRProbsXDays; crprobvals];
        CRAmpsXDays = [CRAmpsXDays; crampvals];
    elseif size(CRProbsXDays,2)<size(crprobvals,2)
        while size(CRProbsXDays,2)<size(crprobvals,2)
            CRProbsXDays = [CRProbsXDays,nan(size(CRProbsXDays,1),1)];
            CRAmpsXDays = [CRAmpsXDays,nan(size(CRAmpsXDays,1),1)];
        end
        CRProbsXDays = [CRProbsXDays; crprobvals];
        CRAmpsXDays = [CRAmpsXDays; crampvals];
    else
        CRProbsXDays = [CRProbsXDays; crprobvals];
        CRAmpsXDays = [CRAmpsXDays; crampvals];
    end
end

colordef white
figure
errorbar([1:12], nanmedian(CRProbsXDays), mad(CRProbsXDays,1), '.', 'LineStyle', 'none', 'Color', [0 0 0])
hold on
errorbar([1:12], nanmedian(CRAmpsXDays), mad(CRAmpsXDays,1), '.', 'LineStyle', 'none', 'Color', [1 0 1])
xlabel('Sessions')
ylabel('CR Probability')

%% Get performance data for the control group animals
basedir = 'E:\pcp2ChR2 data\rebound';
mice = {'OK234';'OK235';'OK236';'OK237',;'OK238';'OK239';'OK240';'OK241'};

daily.mouse = {};
daily.crprob = [];
daily.cradjamp = [];
daily.cradjampHit = [];
daily.sesstype = []; % 0 = training, 1 = extinction, 2 = laser, 3 = laser extinction trial
for m = 1:length(mice)
    mousedir = [basedir,'\',mice{m,1}];
    cd(mousedir)
    
    days = dir('19*');
    for d = 1:length(days)
        
        daydir = [mousedir,'\',days(d,1).name];
        cd(daydir)
        
        % check if there is a trialdata file present
        loaded = false;
        if exist('newTrialdata.mat','file')==2
            load('newTrialdata.mat')
            loaded = true;
        elseif exist('trialdata.mat','file')==2
            load('trialdata.mat')
            loaded = true;
        end
        
        if loaded == true
            % check if this day had CS trials or if it was a laser test day
            numCSTrials = sum(trials.c_csdur>0);
            
            if numCSTrials > 0
                CSUSidx = find(trials.c_csdur>0 & trials.c_usdur>0);
                CSOnlyidx = find(trials.c_csdur>0 & trials.c_usdur==0);
                laserIdx = find(trials.laser.dur>0);               
                if ~isempty(CSUSidx)
                    checkidx = CSUSidx;
                    if length(laserIdx)>length(CSUSidx)                        
                        daily.sesstype(end+1,1) = 3;
                    else
                        daily.sesstype(end+1,1) = 0;
                    end
                else
                    checkidx = CSOnlyidx;
                    daily.sesstype(end+1,1) = 1;
                end
                
                baseline = nan(length(checkidx),1);
                cradjamp = nan(length(checkidx),1);
                stable = nan(length(checkidx),1);
                %eyeadj = nan(length(checkidx),size(trials.eyelidpos,2));
                for t = 1:length(checkidx)
                    baseline(t,1) = mean(trials.eyelidpos(checkidx(t,1),1:40));
                    stable(t,1) = max(trials.eyelidpos(checkidx(t,1),1:40))<0.3;
                    cradjamp(t,1) = max(trials.eyelidpos(checkidx(t,1),76:85))-baseline(t,1);
                    %eyeadj(t,:)=trials.eyelidpos(checkidx(t,1),:)-baseline(t,1);
                end
                
                meancradjamp = median(cradjamp(stable==1,1));
                meancradjampHit = median(cradjamp(stable==1 & cradjamp>=0.1));
                crprob = sum(cradjamp(stable==1,1)>=0.1)./sum(stable==1);
                
%                 figure
%                 plot(eyeadj(stable==1,:)')
%                 hold on
%                 plot(mean(eyeadj(stable==1,:)),'Color',[0 0 0],'LineWidth',3)
%                 title(num2str(meancradjamp))
%                 pause
                
                daily.mouse{end+1,1} = mice{m,1};
                daily.cradjamp(end+1,1) = meancradjamp;
                daily.cradjampHit(end+1,1) = meancradjampHit;
                daily.crprob(end+1,1) = crprob;
                
                clear trials meancradjamp crprob baseline cradjamp stable CSUSidx...
                    checkidx CSOnlyidx
            else
                daily.mouse{end+1,1} = mice{m,1};
                daily.sesstype(end+1,1) = 2;
                daily.cradjamp(end+1,1) = nan;
                daily.crprob(end+1,1) = nan;
                daily.cradjampHit(end+1,1) = nan;
            end
            
            clear numCSTrials
            
        end
        
        clear loaded daydir
    end
    
    clear mousedir
end

% plot control group animals' performance across training
CRProbsXDays = [];
CRAmpsXDays = [];
for m = 1:length(mice)
    thisMouse = mice(m,1);
    idx = [];
    for i = 1:length(daily.mouse)
        if strcmpi(thisMouse, daily.mouse{i,1})
            idx(end+1)=i;
        end
    end
    temp = find(daily.sesstype==1);
    temp2 = ismember(idx,temp);
    testidx = idx(temp2);
    pretestidx = testidx(1);
    posttestidx = testidx(2);
    crprobvals = daily.crprob(pretestidx+1:posttestidx-1,1)';
    crampvals = daily.cradjampHit(pretestidx+1:posttestidx-1,1)';
    if isempty(CRProbsXDays)
        CRProbsXDays = crprobvals;
        CRAmpsXDays = crampvals;
    elseif size(CRProbsXDays,2) > size(crprobvals,2)
        while size(CRProbsXDays,2) > size(crprobvals,2)
            crprobvals = [crprobvals,NaN];
            crampvals = [crampvals,NaN];
        end
        CRProbsXDays = [CRProbsXDays; crprobvals];
        CRAmpsXDays = [CRAmpsXDays; crampvals];
    elseif size(CRProbsXDays,2)<size(crprobvals,2)
        while size(CRProbsXDays,2)<size(crprobvals,2)
            CRProbsXDays = [CRProbsXDays,nan(size(CRProbsXDays,1),1)];
            CRAmpsXDays = [CRAmpsXDays,nan(size(CRAmpsXDays,1),1)];
        end
        CRProbsXDays = [CRProbsXDays; crprobvals];
        CRAmpsXDays = [CRAmpsXDays; crampvals];
    else
        CRProbsXDays = [CRProbsXDays; crprobvals];
        CRAmpsXDays = [CRAmpsXDays; crampvals];
    end
end

colordef white
figure
errorbar([1:19], nanmedian(CRProbsXDays), mad(CRProbsXDays,1), '.', 'LineStyle', 'none', 'Color', [0 0 0])
hold on
errorbar([1:19], nanmedian(CRAmpsXDays), mad(CRAmpsXDays,1), '.', 'LineStyle', 'none', 'Color', [1 0 1])
xlabel('Sessions')
ylabel('CR Probability')
ylim([0 1])