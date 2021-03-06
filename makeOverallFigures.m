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
load('overallData_200504.mat')



rbstats.mouse = [];
rbstats.phase = {};
rbstats.amp = [];
rbstats.lat = [];
rbstats.lasamp = [];
rbstats.laspow = [];
rbstats.winstart = [];

crstats.mouse = [];
crstats.phase = {};
crstats.amp = [];

daycrstats.mouse = [];
daycrstats.phase = {};
daycrstats.amp = [];
daycrstats.prob = [];
daycrstats.meantr = [];

daystats.mouse = [];
daystats.phase = {};
daystats.rb.amp = [];
daystats.rb.hitamp = [];
daystats.rb.lat = [];
daystats.rb.prob = [];
daystats.lasamp = [];
daystats.laspow = [];
daystats.meanRBTr = [];
daystats.meanRBTrHit = [];
daystats.winstart = [];


timeVector = 1:size(postPairedData.eyelidpos,2);
timeVector = timeVector * 0.00488372;
timeVector = timeVector - 0.2;

mice = unique(postPairedData.mouse);

for m = 1:length(mice)
    [rbstats, crstats, daystats, daycrstats] = getrbprops(preUnpairedData,'preUnpaired', rbstats, daystats,...
         crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(lastUnpairedData,'lastUnpaired', rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(postUnpairedData,'postUnpaired', rbstats, daystats,...
         crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(prePairedData,'prePaired', rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(lastAcqData,'lastAcq', rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(postPairedData,'postPaired', rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(lastExtData,'lastExt', rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    [rbstats, crstats, daystats, daycrstats] = getrbprops(postExtData,'postExt', rbstats, daystats,...
        crstats, daycrstats, m, mice, timeVector, 0);
    %plot(daystats.meanRBTr')
end

%% TO DO: get test inhibition data into some usable format

%% compare animals performance within group across light intensities
% NOTE: took out pre-paired data from OK235 because it was breaking the
% signrank function since he had missing data (missed the post paired test for health reasons)
[signrankp_paired, ttestp_paired] = makeSubplotSummary(daystats, 211, 218, 'Paired');
[signrankp_unpaired, ttestp_unpaired] = makeSubplotSummary(daystats, 234, 241, 'Unpaired');

%% need to do something with these to show that can combine the groups. basically no differences between groups so it's ok...?
[prevals, postvals, premad, postmad]=makePlots_RBxPower(daystats, daystats.rb.amp, [15,30,60], ...
    [211,213,214,215,216,217,218],...
    1, 'prePaired', 'postPaired');
xlabel('Laser Power (mW)')
ylabel('Rebound Amplitude (FEC), n=7')
title('initially got paired training')
ylim([0 1])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')
checkShapiroWilk = prevals(:,1)-postvals(:,1); % R shapiro.test: W = 0.85763, p-value = 0.1442
checkShapiroWilk = prevals(:,2)-postvals(:,2); % R shapiro.test: W = 0.86774, p-value = 0.1773
checkShapiroWilk = prevals(:,3)-postvals(:,3); % R shapiro.test: W = 0.91333, p-value = 0.4194
[h, p, ci, stats] = ttest(prevals(:,1),postvals(:,1)) % p = 0.0461
[h, p, ci, stats] = ttest(prevals(:,2),postvals(:,2)) % p = 0.0135
[h, p, ci, stats] = ttest(prevals(:,3),postvals(:,3)) % p = 0.0198

[prevals, postvals2, premad, postmad]=makePlots_RBxPower(daystats, daystats.rb.amp, [15,30,60], ...
    [234,236,237,238,239,240,241],...
    1, 'prePaired', 'postPaired');
xlabel('Laser Power (mW)')
ylabel('Rebound Amplitude (FEC), n=7')
title('previously got unpaired training')
ylim([0 1])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')
checkShapiroWilk = prevals(:,1)-postvals(:,1); % R shapiro.test: W = 0.85763, p-value = 0.1442
checkShapiroWilk = prevals(:,2)-postvals(:,2); % R shapiro.test: W = 0.86774, p-value = 0.1773
checkShapiroWilk = prevals(:,3)-postvals(:,3); % R shapiro.test: W = 0.91333, p-value = 0.4194
[h, p, ci, stats] = ttest(prevals(:,1),postvals(:,1)) % p = 0.0461
[h, p, ci, stats] = ttest(prevals(:,2),postvals(:,2)) % p = 0.0135
[h, p, ci, stats] = ttest(prevals(:,3),postvals(:,3)) % p = 0.0198


[pretestprobs, trainprobs, pretestprobmad, trainprobmad]=makePlots_RBxPower(daystats, daystats.rb.prob, [15,30,60], ...
    [211,213,214,215,216,217,218],...
    1, 'prePaired', 'postPaired');
xlabel('Laser Power (mW)')
ylabel('Rebound Probability (FEC), n = 7')
title('initially got paired training')
ylim([0 1])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')

[pretestprobs, trainprobs, pretestprobmad, trainprobmad]=makePlots_RBxPower(daystats, daystats.rb.prob, [15,30,60], ...
    [234,236,237,238,239,240,241],...
    1, 'prePaired', 'postPaired');
xlabel('Laser Power (mW)')
ylabel('Rebound Probability (FEC), n = 7')
title('initially got unpaired training')
ylim([0 1])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')

%% get data for all mice pre and post training
[prevals, postvals, premad, postmad]=makePlots_RBxPower(daystats, daystats.rb.amp, [15,30,60], ...
    [211,213,214,215,216,217,218,234,236,237,238,239,240,241],...
    1, 'prePaired', 'postPaired');
xlabel('Laser Power (mW)')
ylabel('Rebound Amplitude (FEC), n=14')
title('groups mixed')
ylim([0 0.6])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')
[p,h,stats] = signrank(prevals(:,1),postvals(:,1)) % p = 0.02, alpha is 0.0125
[p,h,stats] = signrank(prevals(:,2),postvals(:,2)) % p = 0.0012
[p,h,stats] = signrank(prevals(:,3),postvals(:,3)) % p = 0.00061035


[pretestprobs, trainprobs, pretestprobmad, trainprobmad]=makePlots_RBxPower(daystats, daystats.rb.prob, [15,30,60], ...
    [211,213,214,215,216,217,218,234,236,237,238,239,240,241],...
    1, 'prePaired', 'postPaired');
xlabel('Laser Power (mW)')
ylabel('Rebound Probability (FEC), n = 7')
title('groups mixed')
ylim([0 1])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')
[p,h,stats] = signrank(pretestprobs(:,1), trainprobs(:,1)) % p = 0.0781, alpha is 0.0125
[p,h,stats] = signrank(pretestprobs(:,2), trainprobs(:,2)) % p = 0.00024414
[p,h,stats] = signrank(pretestprobs(:,3), trainprobs(:,3)) % p = 0.00024414

% plot boxplot, RB amplitudes across all intensities
predata = nanmedian(prevals,2);
postdata = nanmedian(postvals,2);
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
ylim([0 0.6])
for i = 1:14
    plot([1 2],[predata(i) postdata(i)], 'LineStyle', ':', 'Color', [0 0 0])
end
scatter(ones(14,1), predata, 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[0 0 1])
scatter(ones(14,1)*2, postdata, 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[1 0 0])
ylabel('RB amplitude across all intensities')

% plot boxplot, rb amplitudes at highest intensity
prequant = quantile(prevals(:,3),3);
postquant = quantile(postvals(:,3),3);
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
for i = 1:14
    plot([1 2],[prevals(i,3) postvals(i,3)], 'LineStyle', ':', 'Color', [0 0 0])
end
scatter(ones(14,1), prevals(:,3), 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[0 0 1])
scatter(ones(14,1)*2, postvals(:,3), 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[1 0 0])
ylabel('RB amplitude at 60 mW')


% plot boxplot, RB probability across all intensities
    % is this the right way to compute this
predata = nanmedian(pretestprobs,2);
postdata = nanmedian(trainprobs,2);
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
for i = 1:14
    plot([1 2],[predata(i) postdata(i)], 'LineStyle', ':', 'Color', [0 0 0])
end
scatter(ones(14,1), predata, 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[0 0 1])
scatter(ones(14,1)*2, postdata, 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[1 0 0])
ylabel('RB probability across all intensities')

% plot boxplot, rb probability at highest intensity
prequant = quantile(pretestprobs(:,3),3);
postquant = quantile(trainprobs(:,3),3);
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
for i = 1:14
    plot([1 2],[pretestprobs(i,3) trainprobs(i,3)], 'LineStyle', ':', 'Color', [0 0 0])
end
scatter(ones(14,1), pretestprobs(:,3), 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[0 0 1])
scatter(ones(14,1)*2, trainprobs(:,3), 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[1 0 0])
ylabel('RB probability at 60 mW')

% plot boxplot, rb probability at 30 mW
prequant = quantile(pretestprobs(:,2),3);
postquant = quantile(trainprobs(:,2),3);
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
for i = 1:14
    plot([1 2],[pretestprobs(i,2) trainprobs(i,2)], 'LineStyle', ':', 'Color', [0 0 0])
end
scatter(ones(14,1), pretestprobs(:,2), 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[0 0 1])
scatter(ones(14,1)*2, trainprobs(:,2), 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[1 0 0])
ylabel('RB probability at 30 mW')


pause

%% STOPPED CHECKING THINGS HERE

[postTrainHitAmps, postExtHitAmps, ~, ~]=makePlots_RBxPower(daystats, daystats.rb.hitamp, [15,30,60], [211,213,214,215,216,217,218], 1, 2, 3);
xlabel('Laser Power (mW)')
ylabel('Rebound Amplitude (FEC)')
ylim([0 1])
xlim([0.5 3.5])
legend('Post-Test','Post-Test','Post-Ext','Post-Ext', 'Location', 'NorthWest')
[pretestHitAmps, ~, ~, ~]=makePlots_RBxPower(daystats, daystats.rb.hitamp, [15,30,60], [211,213,214,215,216,217,218], 1, 1, nan);


figure
quants30 = quantile(postTrainHitAmps(:,2),3);
quants60 = quantile(postTrainHitAmps(:,3),3);
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
     
     % select intensity with the highest probability of rebound for each
     % mouse
     idxT = nan(size(trainprobs,1),1);
     idxP = nan(size(trainprobs,1),1);
     idxE = nan(size(trainprobs,1),1);
     idxToUse = nan(size(trainprobs,1),1);
     for i = 1:size(trainprobs,1)
         [~,idxT(i,1)]=max(trainprobs(i,:));
         [~,idxP(i,1)]=max(pretestprobs(i,:));
         [~,idxE(i,1)]=max(extprobs(i,:));
         idxToUse(i,1) = max([idxT(i,1),idxP(i,1),idxE(i,1)]);
     end
     
% for experimental group extinction: looking at rebound amplitudes
[trainamps, extamps, trainampmad, extampmad]=makePlots_RBxPower(daystats, daystats.rb.amp, [15,30,60], [211,213,214,215,216,217,218], 1, 2, 3);
xlabel('Laser Power (mW)')
ylabel('Rebound Amplitude (FEC)')
ylim([0 1])
xlim([0.5 3.5])
[pretestamps, ~, pretestampmad, ~]=makePlots_RBxPower(daystats, daystats.rb.amp, [15,30,60], [211,213,214,215,216,217,218], 0, 1, nan);
legend('Post-Training','Post-Training','Post-Extinction','Post-Extinction','Pre-Test','Pre-Test', 'Location', 'NorthWest')

% write a CSV with all the rebound probability and amplitude data
% let col1 be experiment phase, col2 be mouse, col3 be laser power, and col
% 4 be RB probability
theseMice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218'};
powers = [15;30;60];
bigcell = cell(63,6); % 7*3*3 rows (n*phases*powers), 6 columns
iter = 1;
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'pretest';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = pretestprobs(i,c);
        bigcell{iter,5} = pretestamps(i,c);
        bigcell{iter, 6} = pretestHitAmps(i,c);
        iter = iter+1;
    end
end
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'post_training';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = trainprobs(i,c);
        bigcell{iter,5} = trainamps(i,c);
        bigcell{iter, 6} = postTrainHitAmps(i,c);
        iter = iter+1;
    end
end
for i = 1:7
    for c = 1:3
        bigcell{iter,1} = 'post_extinction';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = extprobs(i,c);
        bigcell{iter,5} = extamps(i,c);
        bigcell{iter, 6} = postExtHitAmps(i,c);
        iter = iter+1;
    end
end
bigtable = cell2table(bigcell,'VariableNames',{'phase','mouse','power_mW','RB_prob','RB_amp','RB_hitAmp'});
writetable(bigtable, 'RBData_20200108.csv')

% make a CSV with restricted RB probabilities based on the maximal RB
% evoked at training
theseMice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218'};
powers = [15;30;60];
bigcell = cell(21,4); % 7*3 rows (n*phases*powers), 4 columns
iter = 1;
for i = 1:7
    c = idxT(i,1);
    bigcell{iter,1} = 'pretest';
    bigcell{iter,2} = theseMice(i,1);
    bigcell{iter,3} = powers(c,1);
    bigcell{iter,4} = pretestprobs(i,c);
    bigcell{iter,5} = pretestamps(i,c);
    bigcell{iter, 6} = pretestHitAmps(i,c);
    iter = iter+1;
end
for i = 1:7
    c = idxT(i,1);
    bigcell{iter,1} = 'post_training';
    bigcell{iter,2} = theseMice(i,1);
    bigcell{iter,3} = powers(c,1);
    bigcell{iter,4} = trainprobs(i,c);
    bigcell{iter,5} = trainamps(i,c);
    bigcell{iter, 6} = postTrainHitAmps(i,c);
    iter = iter+1;
end
for i = 1:7
    c = idxT(i,1);
    bigcell{iter,1} = 'post_extinction';
    bigcell{iter,2} = theseMice(i,1);
    bigcell{iter,3} = powers(c,1);
    bigcell{iter,4} = extprobs(i,c);
    bigcell{iter,5} = extamps(i,c);
    bigcell{iter, 6} = postExtHitAmps(i,c);
    iter = iter+1;
end
bigtable = cell2table(bigcell,'VariableNames',{'phase','mouse','power_mW','RB_prob','RB_amp','RB_hitAmp'});
writetable(bigtable, 'RBProbs_max_20200108.csv')

% make figure 4 panel d: pre vs post for training and extinction
pretraini = [];
posttraini = [];
postunpi = [];
for i = 1:length(bigtable.phase)
    if strcmpi('pretest',bigtable.phase{i,1})
        pretraini = [pretraini;i];
    elseif strcmpi('post_training',bigtable.phase{i,1})
        posttraini = [posttraini; i];
    elseif strcmpi('post_extinction',bigtable.phase{i,1})
        postunpi = [postunpi; i];
    end
end
prepaired = bigtable.RB_prob(pretraini,1);
postpaired = bigtable.RB_prob(posttraini,1);
postunpaired = bigtable.RB_prob(postunpi,1);
medianVals = [median(prepaired), median(postpaired), median(postunpaired)];
tempquantPre = quantile(prepaired,3);
tempquantPost = quantile(postpaired,3);
tempquantPostE = quantile(postunpaired,3);
pairedPosErrBars = [tempquantPre(3)-tempquantPre(2), tempquantPost(3)-tempquantPost(2), tempquantPostE(3)-tempquantPostE(2)];
pairedNegErrBars = [tempquantPre(2)-tempquantPre(1), tempquantPost(2)-tempquantPost(1), tempquantPostE(2)-tempquantPostE(1)];
figure
errorbar([1,2,3], medianVals, pairedNegErrBars, pairedPosErrBars, 'Color', [0 0 0])
hold on
scatter([1,2,3],medianVals,10,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
xlim([0.5 3.5])
ylabel('RB Probability')


% make boxplot for the probabilities
tempquantPre = quantile(prepaired,3);
tempquantPost = quantile(postpaired,3);
tempquantPostE = quantile(postunpaired,3);
figure
%plot([1 2 3], [median(prepaired) median(postpaired) median(postpaired)], 'Color', [0 0 0])
hold on
plot([0.75 1.25], [nanmedian(prepaired) nanmedian(prepaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([0.75 1.25], [tempquantPre(3) tempquantPre(3)], 'Color', [0.5 0.5 0.5])
plot([0.75 1.25], [tempquantPre(1) tempquantPre(1)], 'Color', [0.5 0.5 0.5])
plot([1.75 2.25], [nanmedian(postpaired) nanmedian(postpaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([1.75 2.25], [tempquantPost(3) tempquantPost(3)], 'Color', [0.5 0.5 0.5])
plot([1.75 2.25], [tempquantPost(1) tempquantPost(1)], 'Color', [0.5 0.5 0.5])
plot([2.75 3.25], [nanmedian(postunpaired) nanmedian(postunpaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([2.75 3.25], [tempquantPostE(3) tempquantPostE(3)], 'Color', [1 0 1])
plot([2.75 3.25], [tempquantPostE(1) tempquantPostE(1)], 'Color', [1 0 1])
for i = 1:7
    if ~isnan(postunpaired(i,1))
        plot([1 2 3], [prepaired(i) postpaired(i) postunpaired(i)], 'Color', [0 0 1], 'LineStyle', ':')
    end
end
scatter(ones(7,1), prepaired, 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
scatter(ones(7,1)*2, postpaired, 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
scatter(ones(7,1)*3, postunpaired, 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0])
ylabel('Rebound Probability (FEC)')

% make boxplot for the amplitudes
prepaired = bigtable.RB_hitAmp(pretraini,1);
postpaired = bigtable.RB_hitAmp(posttraini,1);
postunpaired = bigtable.RB_hitAmp(postunpi,1);
tempquantPre = quantile(prepaired,3);
tempquantPost = quantile(postpaired,3);
tempquantPostE = quantile(postunpaired,3);
figure
hold on
plot([0.75 1.25], [nanmedian(prepaired) nanmedian(prepaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([0.75 1.25], [tempquantPre(3) tempquantPre(3)], 'Color', [0.5 0.5 0.5])
plot([0.75 1.25], [tempquantPre(1) tempquantPre(1)], 'Color', [0.5 0.5 0.5])
plot([1.75 2.25], [nanmedian(postpaired) nanmedian(postpaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([1.75 2.25], [tempquantPost(3) tempquantPost(3)], 'Color', [0.5 0.5 0.5])
plot([1.75 2.25], [tempquantPost(1) tempquantPost(1)], 'Color', [0.5 0.5 0.5])
plot([2.75 3.25], [nanmedian(postunpaired) nanmedian(postunpaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([2.75 3.25], [tempquantPostE(3) tempquantPostE(3)], 'Color', [1 0 1])
plot([2.75 3.25], [tempquantPostE(1) tempquantPostE(1)], 'Color', [1 0 1])
for i = 1:7
    if ~isnan(postunpaired(i,1))
        plot([1 2 3], [prepaired(i) postpaired(i) postunpaired(i)], 'Color', [0 0 1], 'LineStyle', ':')
    end
end
scatter(ones(7,1), prepaired, 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
scatter(ones(7,1)*2, postpaired, 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
scatter(ones(7,1)*3, postunpaired, 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0])
ylabel('Rebound Amplitude (FEC)')

% make boxplot for the amplitudes, include all trials
prepaired = bigtable.RB_amp(pretraini,1);
postpaired = bigtable.RB_amp(posttraini,1);
postunpaired = bigtable.RB_amp(postunpi,1);
tempquantPre = quantile(prepaired,3);
tempquantPost = quantile(postpaired,3);
tempquantPostE = quantile(postunpaired,3);
figure
hold on
plot([0.75 1.25], [nanmedian(prepaired) nanmedian(prepaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([0.75 1.25], [tempquantPre(3) tempquantPre(3)], 'Color', [0.5 0.5 0.5])
plot([0.75 1.25], [tempquantPre(1) tempquantPre(1)], 'Color', [0.5 0.5 0.5])
plot([1.75 2.25], [nanmedian(postpaired) nanmedian(postpaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([1.75 2.25], [tempquantPost(3) tempquantPost(3)], 'Color', [0.5 0.5 0.5])
plot([1.75 2.25], [tempquantPost(1) tempquantPost(1)], 'Color', [0.5 0.5 0.5])
plot([2.75 3.25], [nanmedian(postunpaired) nanmedian(postunpaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([2.75 3.25], [tempquantPostE(3) tempquantPostE(3)], 'Color', [1 0 1])
plot([2.75 3.25], [tempquantPostE(1) tempquantPostE(1)], 'Color', [1 0 1])
for i = 1:7
    if ~isnan(postunpaired(i,1))
        plot([1 2 3], [prepaired(i) postpaired(i) postunpaired(i)], 'Color', [0 0 1], 'LineStyle', ':')
    end
end
scatter(ones(7,1), prepaired, 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
scatter(ones(7,1)*2, postpaired, 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
scatter(ones(7,1)*3, postunpaired, 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0])
ylabel('Rebound Amplitude (FEC)')

% for control mice
[preunpprob, postunpprob, preunpmad, postunpmad]=makePlots_RBxPower(daystats, daystats.rb.prob, [15,30,60], [234,235,236,237,238,239,240,241], 1, 1, 2);
xlabel('Laser Power (mW)')
ylabel('Rebound Probability (FEC)')
ylim([0 1])
xlim([0.5 3.5])
legend('Pre-Test','Pre-Test','Post-Test','Post-Test', 'Location', 'NorthWest')

idxT_cont = nan(size(postunpprob,1),1);
idxP_cont = nan(size(postunpprob,1),1);
for i = 1:size(postunpprob,1)
    [~,idxT_cont(i,1)]=max(postunpprob(i,:));
    [~,idxP_cont(i,1)]=max(preunpprob(i,:));
end


% make a CSV with restricted RB probabilities based on the maximal RB
% evoked at the posttest
theseMice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218'};
powers = [15;30;60];
bigcell = cell(28,4); % 7*3 rows (n*phases*powers), 4 columns
iter = 1;
for i = 1:7
    c = idxT(i,1);
    bigcell{iter,1} = 'pre_training';
    bigcell{iter,2} = theseMice(i,1);
    bigcell{iter,3} = powers(c,1);
    bigcell{iter,4} = pretestprobs(i,c);
    iter = iter+1;
end
for i = 1:7
    c = idxT(i,1);
    bigcell{iter,1} = 'post_training';
    bigcell{iter,2} = theseMice(i,1);
    bigcell{iter,3} = powers(c,1);
    bigcell{iter,4} = trainprobs(i,c);
    iter = iter+1;
end
theseMice = {'OK234';'OK235';'OK236';'OK237';'OK238';'OK239';'OK241';'OK242'};
for i = 1:7
    c = idxT_cont(i,1);
    bigcell{iter,1} = 'pre_unpaired';
    bigcell{iter,2} = theseMice(i,1);
    bigcell{iter,3} = powers(c,1);
    bigcell{iter,4} = preunpprob(i,c);
    iter = iter+1;
end
for i = 1:7
    c = idxT_cont(i,1);
    bigcell{iter,1} = 'post_unpaired';
    bigcell{iter,2} = theseMice(i,1);
    bigcell{iter,3} = powers(c,1);
    bigcell{iter,4} = postunpprob(i,c);
    iter = iter+1;
end
bigtable = cell2table(bigcell,'VariableNames',{'phase','mouse','power_mW','RB_prob'});
writetable(bigtable, 'RBProbs_expVsCont_forFig2_20200108.csv')

% make figure 3 panel c: pre vs post for paired and unpaired training
% find different phase indices
pretraini = [];
posttraini = [];
preunpi = [];
postunpi = [];
for i = 1:length(bigtable.phase)
    if strcmpi('pre_training',bigtable.phase{i,1})
        pretraini = [pretraini;i];
    elseif strcmpi('post_training',bigtable.phase{i,1})
        posttraini = [posttraini; i];
    elseif strcmpi('pre_unpaired',bigtable.phase{i,1})
        preunpi = [preunpi; i];
    elseif strcmpi('post_unpaired',bigtable.phase{i,1})
        postunpi = [postunpi; i];
    end
end
prepaired = bigtable.RB_prob(pretraini,1);
postpaired = bigtable.RB_prob(posttraini,1);
preunpaired = bigtable.RB_prob(preunpi,1);
postunpaired = bigtable.RB_prob(postunpi,1);
pairedVals = [median(prepaired), median(postpaired)];
tempquantPre = quantile(prepaired,3);
tempquantPost = quantile(postpaired,3);
pairedPosErrBars = [tempquantPre(3)-tempquantPre(2), tempquantPost(3)-tempquantPost(2)];
pairedNegErrBars = [tempquantPre(2)-tempquantPre(1), tempquantPost(2)-tempquantPost(1)];
unpairedVals = [median(preunpaired), median(postunpaired)];
tempquantPre = quantile(preunpaired,3);
tempquantPost = quantile(postunpaired,3);
unpairedPosErrBars = [tempquantPre(3)-tempquantPre(2), tempquantPost(3)-tempquantPost(2)];
unpairedNegErrBars = [tempquantPre(2)-tempquantPre(1), tempquantPost(2)-tempquantPost(1)];
figure
errorbar([1,2], pairedVals, pairedNegErrBars, pairedPosErrBars, 'Color', [0 0 0])
hold on
scatter([1,2],pairedVals,10,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
errorbar([1,2], unpairedVals,  unpairedNegErrBars, unpairedPosErrBars, 'Color', [1 0 0])
scatter([1,2],unpairedVals,10,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0])
xlim([0.5 2.5])
ylabel('RB Probability')

% boxplot data for just the unpaired training animals
tempquantPre = quantile(preunpaired,3);
tempquantPost = quantile(postunpaired,3);
figure
%plot([1 2 3], [median(prepaired) median(postpaired) median(postpaired)], 'Color', [0 0 0])
hold on
plot([0.75 1.25], [nanmedian(preunpaired) nanmedian(preunpaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([0.75 1.25], [tempquantPre(3) tempquantPre(3)], 'Color', [0.5 0.5 0.5])
plot([0.75 1.25], [tempquantPre(1) tempquantPre(1)], 'Color', [0.5 0.5 0.5])
plot([1.75 2.25], [nanmedian(postunpaired) nanmedian(postunpaired)], 'Color', [1 0 0], 'LineWidth', 2)
plot([1.75 2.25], [tempquantPost(3) tempquantPost(3)], 'Color', [0.5 0.5 0.5])
plot([1.75 2.25], [tempquantPost(1) tempquantPost(1)], 'Color', [0.5 0.5 0.5])
for i = 1:7
    if ~isnan(postunpaired(i,1))
        plot([1 2], [preunpaired(i) postunpaired(i)], 'Color', [0 0 1], 'LineStyle', ':')
    end
end
scatter(ones(7,1), preunpaired, 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
scatter(ones(7,1)*2, postunpaired, 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
ylabel('Rebound Probability (FEC)')
ylim([0 1])

% plot the average pre and post eyelid traces for the compared
% intensities/times
theseMice = [211;213;214;215;216;217;218];
pretrain_RBTr = [];
pretrain_RBTrHit = [];
posttrain_RBTr = [];
posttrain_RBTrHit = [];
postext_RBTr = [];
postext_RBTrHit = [];
for m = 1:length(theseMice)
    arrayidx = daystats.mouse==theseMice(m,1);
    temppows = daystats.laspow(arrayidx);
    tempphases = daystats.phase(arrayidx);
    temptracesHit = daystats.meanRBTrHit(arrayidx,:);
    temptracesAll = daystats.meanRBTr(arrayidx,:);
    
    if idxT(m,1)==3
        thisPow = 60;
    elseif idxT(m,1)==2
        thisPow = 30;
    elseif idxT(m,1)==2
        thisPow = 15;
    end
    
    % get pre traces
    preidx = tempphases==1 & temppows == thisPow;
    pretrain_RBTr = [pretrain_RBTr; temptracesAll(preidx,:)];
    pretrain_RBTrHit = [pretrain_RBTrHit; temptracesHit(preidx,:)];
    
    % get post traces
    postidx = tempphases==2 & temppows==thisPow;
    posttrain_RBTr = [posttrain_RBTr; temptracesAll(postidx,:)];
    posttrain_RBTrHit = [posttrain_RBTrHit; temptracesHit(postidx,:)];
    
    % get ext traces
    postEidx = tempphases==3 & temppows==thisPow;
    postext_RBTr = [postext_RBTr; temptracesAll(postEidx,:)];
    postext_RBTrHit = [postext_RBTrHit; temptracesHit(postEidx,:)];
    
    clear thisPow
end
figure
subplot(2,2,1)
shadedErrorBar(timeVector, median(posttrain_RBTr), mad(posttrain_RBTr,1))
hold on
shadedErrorBar(timeVector, median(postext_RBTr), mad(postext_RBTr,1), '-b')
shadedErrorBar(timeVector, median(pretrain_RBTr), mad(pretrain_RBTr,1),'-r')
plot([0 0],[0 1],'LineStyle',':','Color',[0 0 0])
plot([0.85 0.85],[0 1], 'LineStyle',':','Color',[0 0 0])
ylabel('Eyelid Position (FEC')
xlabel('Time from Laser Onset (s)')
ylim([0 1])
xlim([-0.1 1.450])
title('median +/- mad, all trials')

subplot(2,2,2)
shadedErrorBar(timeVector, mean(posttrain_RBTr), std(posttrain_RBTr)./sqrt(size(posttrain_RBTr,1)))
hold on
shadedErrorBar(timeVector, mean(postext_RBTr), std(postext_RBTr)./sqrt(size(posttrain_RBTr,1)), '-b')
shadedErrorBar(timeVector, mean(pretrain_RBTr), std(pretrain_RBTr)./sqrt(size(posttrain_RBTr,1)),'-r')
plot([0 0],[0 1],'LineStyle',':','Color',[0 0 0])
plot([0.85 0.85],[0 1], 'LineStyle',':','Color',[0 0 0])
ylabel('Eyelid Position (FEC')
xlabel('Time from Laser Onset (s)')
ylim([0 1])
xlim([-0.1 1.450])
title('mean +/- std, all trials')

subplot(2,2,3)
shadedErrorBar(timeVector, median(posttrain_RBTrHit), mad(posttrain_RBTrHit,1))
hold on
shadedErrorBar(timeVector, nanmedian(postext_RBTrHit), mad(postext_RBTrHit,1), '-b')
shadedErrorBar(timeVector, nanmedian(pretrain_RBTrHit), mad(pretrain_RBTrHit,1),'-r')
plot([0 0],[0 1],'LineStyle',':','Color',[0 0 0])
plot([0.85 0.85],[0 1], 'LineStyle',':','Color',[0 0 0])
ylabel('Eyelid Position (FEC')
xlabel('Time from Laser Onset (s)')
ylim([0 1])
xlim([-0.1 1.450])
title('median +/- mad, hit trials')

subplot(2,2,4)
shadedErrorBar(timeVector, mean(posttrain_RBTrHit), std(posttrain_RBTrHit)./sqrt(sum(~isnan(posttrain_RBTrHit))))
hold on
shadedErrorBar(timeVector, nanmean(postext_RBTrHit), nanstd(postext_RBTrHit)./sqrt(sum(~isnan(posttrain_RBTrHit))), '-b')
shadedErrorBar(timeVector, nanmean(pretrain_RBTrHit), nanstd(pretrain_RBTrHit)./sqrt(sum(~isnan(posttrain_RBTrHit))),'-r')
plot([0 0],[0 1],'LineStyle',':','Color',[0 0 0])
plot([0.85 0.85],[0 1], 'LineStyle',':','Color',[0 0 0])
ylabel('Eyelid Position (FEC')
xlabel('Time from Laser Onset (s)')
ylim([0 1])
xlim([-0.1 1.450])
title('mean +/- std, hit trials')


% now for controls
theseMice = [234;235;236;237;238;239;240;241];
pretrain_RBTr = [];
pretrain_RBTrHit = [];
posttrain_RBTr = [];
posttrain_RBTrHit = [];
for m = 1:length(theseMice)
    arrayidx = daystats.mouse==theseMice(m,1);
    temppows = daystats.laspow(arrayidx);
    tempphases = daystats.phase(arrayidx);
    temptracesHit = daystats.meanRBTrHit(arrayidx,:);
    temptracesAll = daystats.meanRBTr(arrayidx,:);
    
    if idxT_cont(m,1)==3
        thisPow = 60;
    elseif idxT_cont(m,1)==2
        thisPow = 30;
    elseif idxT_cont(m,1)==1
        thisPow = 15;
    end
    
    % get pre traces
    preidx = tempphases==1 & temppows == thisPow;
    pretrain_RBTr = [pretrain_RBTr; temptracesAll(preidx,:)];
    pretrain_RBTrHit = [pretrain_RBTrHit; temptracesHit(preidx,:)];
    
    % get post traces
    postidx = tempphases==2 & temppows==thisPow;
    posttrain_RBTr = [posttrain_RBTr; temptracesAll(postidx,:)];
    posttrain_RBTrHit = [posttrain_RBTrHit; temptracesHit(postidx,:)];
    
    clear thisPow
end
figure
shadedErrorBar(timeVector, median(posttrain_RBTr), mad(posttrain_RBTr,1))
hold on
shadedErrorBar(timeVector, median(pretrain_RBTr), mad(pretrain_RBTr,1),'-r')
plot([0 0],[0 1],'LineStyle',':','Color',[0 0 0])
plot([0.85 0.85],[0 1], 'LineStyle',':','Color',[0 0 0])
ylabel('Eyelid Position (FEC')
xlabel('Time from Laser Onset (s)')
ylim([0 1])

% and now for only successful rebounds at 30 and 60 mw post training
theseMice = [211;213;214;215;216;217;218];
RBTrHit30 = [];
RBTrHit60 = [];
for m = 1:length(theseMice)
    arrayidx = daystats.mouse==theseMice(m,1);
    temppows = daystats.laspow(arrayidx);
    tempphases = daystats.phase(arrayidx);
    temptracesHit = daystats.meanRBTrHit(arrayidx,:);
    temptracesAll = daystats.meanRBTr(arrayidx,:);
    
    % get pre traces
    idx30 = tempphases==2 & temppows == 30;
    RBTrHit30 = [RBTrHit30; temptracesHit(idx30,:)];
    
    % get post traces
    idx60 = tempphases==2 & temppows== 60;
    RBTrHit60 = [RBTrHit60; temptracesHit(idx60,:)];
    
    clear thisPow
end
figure
shadedErrorBar(timeVector, median(RBTrHit60), mad(RBTrHit60,1))
hold on
shadedErrorBar(timeVector, median(RBTrHit30), mad(RBTrHit30,1),'-r')
plot([0 0],[0 1],'LineStyle',':','Color',[0 0 0])
plot([0.85 0.85],[0 1], 'LineStyle',':','Color',[0 0 0])
ylabel('Eyelid Position (FEC')
xlabel('Time from Laser Onset (s)')
ylim([0 1])


% plot the average pre and post CRs for training
theseMice = [211;213;214;215;216;217;218];
pretrain_Tr = [];
posttrain_Tr = [];
for m = 1:length(theseMice)
    arrayidx = daycrstats.mouse==theseMice(m,1);
    tempphases = daycrstats.phase(arrayidx);
    temptracesAll = daycrstats.meantr(arrayidx,:);
    
    % get pre traces
    goHere = ['E:\pcp2ChR2 data\rebound\OK', num2str(theseMice(m,1))];
    cd(goHere)
    folders = dir;
    goHere = [goHere,'\',folders(4,1).name];
    cd(goHere)
    load('trialdata.mat')
    csusidx = trials.c_csdur>0 & trials.c_usdur>0;
    thisTrace = mean(trials.eyelidpos(csusidx,:));
    if size(thisTrace,2)<340
        addnan = nan(1,340-size(thisTrace,2));
        thisTrace = [thisTrace, addnan];
    end
    pretrain_Tr = [pretrain_Tr; thisTrace];
    
    % get post traces
    postidx = tempphases==1.5;
    posttrain_Tr = [posttrain_Tr; temptracesAll(postidx,:)];
    
    clear thisTrace
end
figure
shadedErrorBar(timeVector, median(pretrain_Tr), mad(pretrain_Tr,1))
hold on
shadedErrorBar(timeVector, median(posttrain_Tr), mad(posttrain_Tr,1),'-r')
plot([0 0],[0 1],'LineStyle',':','Color',[0 0 0])
plot([0.85 0.85],[0 1], 'LineStyle',':','Color',[0 0 0])
ylabel('Eyelid Position (FEC')
xlabel('Time from Laser Onset (s)')
ylim([0 1])


% plot the average pre and post CRs for unpairedtraining
theseMice = [234;235;236;237;238;239;240;241];
pretrain_Tr = [];
posttrain_Tr = [];
for m = 1:length(theseMice)
    arrayidx = daycrstats.mouse==theseMice(m,1);
    tempphases = daycrstats.phase(arrayidx);
    temptracesAll = daycrstats.meantr(arrayidx,:);
    
    % get pre traces
    goHere = ['E:\pcp2ChR2 data\rebound\OK', num2str(theseMice(m,1))];
    cd(goHere)
    folders = dir;
    goHere = [goHere,'\',folders(4,1).name];
    cd(goHere)
    try
        load('newtrialdata.mat')
    catch ME
        load('trialdata.mat')
    end
    csusidx = trials.c_csdur>0 & trials.c_usdur==0;
    thisTrace = mean(trials.eyelidpos(csusidx,:));
    if size(thisTrace,2)<340
        addnan = nan(1,340-size(thisTrace,2));
        thisTrace = [thisTrace, addnan];
    end
    pretrain_Tr = [pretrain_Tr; thisTrace];
    
    % get post traces
    postidx = tempphases==1.5;
    posttrain_Tr = [posttrain_Tr; temptracesAll(postidx,:)];
    
    clear thisTrace
end
figure
shadedErrorBar(timeVector, median(pretrain_Tr), mad(pretrain_Tr,1))
hold on
shadedErrorBar(timeVector, median(posttrain_Tr), mad(posttrain_Tr,1),'-r')
plot([0 0],[0 1],'LineStyle',':','Color',[0 0 0])
plot([0.85 0.85],[0 1], 'LineStyle',':','Color',[0 0 0])
ylabel('Eyelid Position (FEC')
xlabel('Time from Laser Onset (s)')
ylim([0 1])

% save a csv for comparing rebounds after training in the experimental and
% the control animals
theseMice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218'};
powers = [15;30;60];
bigcell = cell(28,4); % 7*2*3 rows (n*phases*powers), 4 columns
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
        bigcell{iter,1} = 'pre_unpaired';
        bigcell{iter,2} = theseMice(i,1);
        bigcell{iter,3} = powers(c,1);
        bigcell{iter,4} = preunpprob(i,c);
        iter = iter+1;
    end
end
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


% plot the average pre and post extinction CR traces
theseMice = [211;213;214;215;216;217;218];
pretrain_Tr = [];
posttrain_Tr = [];
for m = 1:length(theseMice)
    arrayidx = daycrstats.mouse==theseMice(m,1);
    tempphases = daycrstats.phase(arrayidx);
    temptracesAll = daycrstats.meantr(arrayidx,:);
    
    % get pre traces
    preidx = tempphases==1.5;
    pretrain_Tr = [pretrain_Tr; temptracesAll(preidx,:)];
    
    % get post traces
    postidx = tempphases==2.5;
    posttrain_Tr = [posttrain_Tr; temptracesAll(postidx,:)];
    
    clear thisTrace
end
figure
shadedErrorBar(timeVector, median(pretrain_Tr), mad(pretrain_Tr,1))
hold on
shadedErrorBar(timeVector, median(posttrain_Tr), mad(posttrain_Tr,1),'-r')
plot([0 0],[0 1],'LineStyle',':','Color',[0 0 0])
plot([0.85 0.85],[0 1], 'LineStyle',':','Color',[0 0 0])
ylabel('Eyelid Position (FEC')
xlabel('Time from Laser Onset (s)')
ylim([0 1])

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
daily.eyetrace = [];
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
                eyeadj = nan(length(checkidx),size(trials.eyelidpos,2));
                for t = 1:length(checkidx)
                    baseline(t,1) = mean(trials.eyelidpos(checkidx(t,1),1:40));
                    stable(t,1) = max(trials.eyelidpos(checkidx(t,1),1:40))<0.3;
                    cradjamp(t,1) = max(trials.eyelidpos(checkidx(t,1),76:85))-baseline(t,1);
                    eyeadj(t,:)=trials.eyelidpos(checkidx(t,1),:)-baseline(t,1);
                end
                
                meancradjamp = median(cradjamp(stable==1,1));
                meancradjampHit = median(cradjamp(stable==1 & cradjamp>=0.1));
                crprob = sum(cradjamp(stable==1,1)>=0.1)./sum(stable==1);
                eyetrace = mean(eyeadj(stable==1,:));
                
                if length(eyetrace)<440
                    diffr = 440 - length(eyetrace);
                    eyetrace = [eyetrace, nan(1,diffr)];
                end
                
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
                daily.eyetrace = [daily.eyetrace; eyetrace];
                
                
                clear trials meancradjamp crprob baseline cradjamp stable CSUSidx...
                    checkidx CSOnlyidx
            else
                daily.mouse{end+1,1} = mice{m,1};
                daily.sesstype(end+1,1) = 2;
                daily.cradjamp(end+1,1) = nan;
                daily.crprob(end+1,1) = nan;
                daily.cradjampHit(end+1,1) = nan;
                daily.eyetrace = [daily.eyetrace; nan(1, 440)];
            end
            
            clear numCSTrials
            
        end
        
        clear loaded daydir
    end
    
    clear mousedir
end

traceData.tday1 = [];
traceData.tday2 = [];
traceData.tday4 = [];
traceData.tday6 = [];
traceData.tday8 = [];
traceData.tday10 = [];
traceData.tday12 = [];
traceData.tdayLast = [];

traceData.tdayEBaseline = [];
traceData.eday1 = [];
traceData.eday2 = [];
traceData.eday4 = [];
traceData.eday6 = [];
traceData.eday8 = [];
traceData.eday10 = [];
for m = 1:length(mice)
    thisMouse = mice(m,1);
    idx = [];
    for i = 1:length(daily.mouse)
        if strcmpi(thisMouse, daily.mouse{i,1})
            idx(end+1)=i;
        end
    end
    thisSesstype = daily.sesstype(idx,1);
    thisEyetrace = daily.eyetrace(idx,:);
    
    tempidx = find(thisSesstype==0);
    tempdiffs = diff(tempidx);
    temptempidx = find(tempdiffs>1);
    stopidx = temptempidx(1);
    trainEyetraces = thisEyetrace(tempidx(1:stopidx),:);
    traceData.tday1 = [traceData.tday1; trainEyetraces(1,:)];
    traceData.tday2 = [traceData.tday2; trainEyetraces(2,:)];
    traceData.tday4 = [traceData.tday4; trainEyetraces(4,:)];
    traceData.tday6 = [traceData.tday6; trainEyetraces(6,:)];
    traceData.tday8 = [traceData.tday8; trainEyetraces(8,:)];
    traceData.tday10 = [traceData.tday10; trainEyetraces(10,:)];
    traceData.tday12 = [traceData.tday12; trainEyetraces(12,:)];
    traceData.tdayLast = [traceData.tdayLast; trainEyetraces(end,:)];
    
    traceData.tdayEBaseline = [traceData.tdayEBaseline; thisEyetrace(tempidx(temptempidx(2)),:)];
    
    tempidx = find(thisSesstype==1);
    extEyetraces = thisEyetrace(tempidx,:);
    traceData.eday1 = [traceData.eday1; extEyetraces(1,:)];
    traceData.eday2 = [traceData.eday2; extEyetraces(2,:)];
    traceData.eday4 = [traceData.eday4; extEyetraces(4,:)];
    traceData.eday6 = [traceData.eday6; extEyetraces(6,:)];
    traceData.eday8 = [traceData.eday8; extEyetraces(8,:)];
    traceData.eday10 = [traceData.eday10; extEyetraces(10,:)];
    clear idx
end


figure
plot(timeVector, nanmean(traceData.tday1))
hold on
plot(timeVector, nanmean(traceData.tday2))
plot(timeVector, nanmean(traceData.tday4))
plot(timeVector, nanmean(traceData.tday6))
plot(timeVector, nanmean(traceData.tday8))
plot(timeVector, nanmean(traceData.tday10))
plot(timeVector, nanmean(traceData.tday12))
plot(timeVector, nanmean(traceData.tdayLast))
xlim([-0.05 0.3])
legend('day1','2','4','6','8','10','12','last','Location','NorthWest')
ylabel('FEC-baseline')
xlabel('time from tone (s)')
title('training traces')

figure
plot(timeVector, nanmean(traceData.tdayEBaseline))
hold on
plot(timeVector, nanmean(traceData.eday1))
plot(timeVector, nanmean(traceData.eday2))
plot(timeVector, nanmean(traceData.eday4))
plot(timeVector, nanmean(traceData.eday6))
plot(timeVector, nanmean(traceData.eday8))
plot(timeVector, nanmean(traceData.eday10))
xlim([-0.05 0.3])
ylim([0 1])
legend('baseline','day1','2','4','6','8','10','Location','NorthWest')
ylabel('FEC')
xlabel('time from tone (s)')
title('extinction traces')


% let performance threshold be >80% CRs and >50% median amp of hit trials
threshProb = 0.8;
threshAmp = 0.5; %0.5 if use col 4, this is more accurate for asymptotes though
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
        if tempArr(iteridx,2)>=threshProb && tempArr(iteridx,4)>=threshAmp
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
        if tempArr(iteridx,2)>=threshProb && tempArr(iteridx,4)>=threshAmp
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

figure
scatter(daystats.rb.prob(rbidx), numdays.reacq)
xlabel('RB probability after extinction')
ylabel('days to reacq')
ylim([0 5.5])
xlim([0 1])
temprb = daystats.rb.prob(rbidx);
[r, p] = corr(temprb, numdays.reacq);
lsline
titlestring = ['Pearsons r = ', num2str(r), '; p = ', num2str(p)];
title(titlestring)

figure
scatter(numdays.reacq, daystats.rb.amp(rbidx))
ylabel('RB amplitude after extinction')
xlabel('days to reacq')
xlim([0 5.5])
ylim([0 0.5])
temprb = daystats.rb.amp(rbidx);
[r, p] = corr(numdays.reacq, temprb);
lsline
titlestring = ['Pearsons r = ', num2str(r), '; p = ', num2str(p)];
title(titlestring)
% should probably try trials to savings instead of days

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
daily.eyetrace = [];
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
                eyeadj = nan(length(checkidx),size(trials.eyelidpos,2));
                for t = 1:length(checkidx)
                    baseline(t,1) = mean(trials.eyelidpos(checkidx(t,1),1:40));
                    stable(t,1) = max(trials.eyelidpos(checkidx(t,1),1:40))<0.3;
                    cradjamp(t,1) = max(trials.eyelidpos(checkidx(t,1),76:85))-baseline(t,1);
                    eyeadj(t,:)=trials.eyelidpos(checkidx(t,1),:)-baseline(t,1);
                end
                
                meancradjamp = median(cradjamp(stable==1,1));
                meancradjampHit = median(cradjamp(stable==1 & cradjamp>=0.1));
                crprob = sum(cradjamp(stable==1,1)>=0.1)./sum(stable==1);
                eyetrace = mean(eyeadj(stable==1,:));
                
                if length(eyetrace)<440
                    diffr = 440 - length(eyetrace);
                    eyetrace = [eyetrace, nan(1,diffr)];
                end
                
                daily.mouse{end+1,1} = mice{m,1};
                daily.cradjamp(end+1,1) = meancradjamp;
                daily.cradjampHit(end+1,1) = meancradjampHit;
                daily.crprob(end+1,1) = crprob;
                daily.eyetrace = [daily.eyetrace; eyetrace];
                
                clear trials meancradjamp crprob baseline cradjamp stable CSUSidx...
                    checkidx CSOnlyidx
            else
                daily.mouse{end+1,1} = mice{m,1};
                daily.sesstype(end+1,1) = 2;
                daily.cradjamp(end+1,1) = nan;
                daily.crprob(end+1,1) = nan;
                daily.cradjampHit(end+1,1) = nan;
                daily.eyetrace = [daily.eyetrace; nan(1,440)];
            end
            
            clear numCSTrials
            
        end
        
        clear loaded daydir
    end
    
    clear mousedir
end

traceData.tday1 = [];
traceData.tday2 = [];
traceData.tday4 = [];
traceData.tday6 = [];
traceData.tday8 = [];
traceData.tday10 = [];
traceData.tdat12 = [];
traceData.tdayLast = [];
for m = 1:length(mice)
    thisMouse = mice(m,1);
    idx = [];
    for i = 1:length(daily.mouse)
        if strcmpi(thisMouse, daily.mouse{i,1})
            idx(end+1)=i;
        end
    end
    thisSesstype = daily.sesstype(idx,1);
    thisEyetrace = daily.eyetrace(idx,:);
    
    tempidx = find(thisSesstype==1);
    trainEyetraces = thisEyetrace(tempidx(1:stopidx),:);
    traceData.tday1 = [traceData.tday1; trainEyetraces(1,:)];
    traceData.tday2 = [traceData.tday2; trainEyetraces(2,:)];
    traceData.tday4 = [traceData.tday4; trainEyetraces(4,:)];
    traceData.tday6 = [traceData.tday6; trainEyetraces(6,:)];
    traceData.tday8 = [traceData.tday8; trainEyetraces(8,:)];
    traceData.tday10 = [traceData.tday10; trainEyetraces(10,:)];
    traceData.tday12 = [traceData.tday10; trainEyetraces(12,:)];
    traceData.tdayLast = [traceData.tdayLast; trainEyetraces(end,:)];
    clear idx
end


figure
plot(timeVector, nanmean(traceData.tday1))
hold on
plot(timeVector, nanmean(traceData.tday2))
plot(timeVector, nanmean(traceData.tday4))
plot(timeVector, nanmean(traceData.tday6))
plot(timeVector, nanmean(traceData.tday8))
plot(timeVector, nanmean(traceData.tday10))
plot(timeVector, nanmean(traceData.tday12))
plot(timeVector, nanmean(traceData.tdayLast))
xlim([-0.05 0.3])
legend('day1','2','4','6','8','10','12','last','Location','NorthWest')
ylabel('FEC-baseline')
xlabel('time from tone (s)')
title('unpaired training traces')
ylim([0 1])

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