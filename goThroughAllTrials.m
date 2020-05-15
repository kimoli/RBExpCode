clear all
close all

%machine = 'COMPUPITAR';
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

% load the spreadsheet specifying which days belong to which manipulation
% sessions
[~, ~, raw] = xlsread('RBExptDates.xlsx');
xlsdata = raw(1:20,1:10);

% put the dates into an orderly structure
dates.mouse = {};
dates.testInhibDate = nan(19,1); % 19 total subjects
dates.preunpDate = nan(19,1);
dates.lastunpDate = nan(19,1);
dates.postunpDate = nan(19,1);
dates.prepairDate = nan(19,1); % same idx for posttraining, just 1 column over
dates.lastacqDate = nan(19,1);
dates.postpairDate = nan(19,1);
dates.lastextDate = nan(19,1);
dates.postextDate = nan(19,1); % the first group did not get extinction, remove 4
for m = 2:20 % cycle through each mouse (each row of spreadsheet)
    dates.mouse{m-1,1}=xlsdata{m,1};
    for c = 2:10
        switch xlsdata{1,c}
            case 'Inhib Block CR'
                dates.testInhibDate(m-1,1)=cell2mat(xlsdata(m,c));
            case 'Pre-Unpaired'
                dates.preunpDate(m-1,1)=cell2mat(xlsdata(m,c));
            case 'Last Unp'
                dates.lastunpDate(m-1,1)=cell2mat(xlsdata(m,c));
            case 'Post-Unpaired'
                dates.postunpDate(m-1,1)=cell2mat(xlsdata(m,c));
            case 'Pre-Paired'
                dates.prepairDate(m-1,1)=cell2mat(xlsdata(m,c));
            case 'Last Acq'
                dates.lastacqDate(m-1,1)=cell2mat(xlsdata(m,c));
            case 'Post-Paired'
                dates.postpairDate(m-1,1)=cell2mat(xlsdata(m,c));
            case 'Last Ext'
                dates.lastextDate(m-1,1)=cell2mat(xlsdata(m,c));
            case 'Post-Extinction'
                dates.postextDate(m-1,1)=cell2mat(xlsdata(m,c));
        end
    end
end

% user specifies which mice to include in the dataset
mice = {'OK213';'OK211';'OK214';'OK215';'OK216';'OK217';'OK218';'OK234';...
    'OK235';'OK236';'OK237';'OK238';'OK239';'OK240';'OK241'};

% load all the trialdata
load('RBExpt_allTrials.mat')

% make corrected time vector
timeVector = 1:size(data.eyelidpos,2);
timeVector = timeVector * 0.00488372;
timeVector = timeVector - 0.2;


paired.mouse = {};
paired.session = [];
paired.crprob = [];
paired.cradjamp = [];
paired.rbamp = [];
paired.rbprob = [];
paired.rbtrace = [];
paired.eyelidposadj = [];

unpaired.mouse = {};
unpaired.session = [];
unpaired.crprob = [];
unpaired.cradjamp = [];
unpaired.rbamp = [];
unpaired.rbprob = [];
unpaired.rbtrace = [];
unpaired.eyelidposadj = [];

extinction.mouse = {};
extinction.session = [];
extinction.crprob = [];
extinction.cradjamp = [];
extinction.rbamp = [];
extinction.rbprob = [];
extinction.rbtrace = [];
extinction.eyelidposadj = [];

rebound.prePaired.mouse = {};
rebound.prePaired.laserint = [];
rebound.prePaired.rbamp = [];
rebound.prePaired.rbprob = [];
rebound.prePaired.rbtrace = [];

rebound.postPaired.mouse = {};
rebound.postPaired.laserint = [];
rebound.postPaired.rbamp = [];
rebound.postPaired.rbprob = [];
rebound.postPaired.rbtrace = [];

rebound.preUnp.mouse = {};
rebound.preUnp.laserint = [];
rebound.preUnp.rbamp = [];
rebound.preUnp.rbprob = [];
rebound.preUnp.rbtrace = [];

rebound.postUnp.mouse = {};
rebound.postUnp.laserint = [];
rebound.postUnp.rbamp = [];
rebound.postUnp.rbprob = [];
rebound.postUnp.rbtrace = [];

rebound.postExt.mouse = {};
rebound.postExt.laserint = [];
rebound.postExt.rbamp = [];
rebound.postExt.rbprob = [];
rebound.postExt.rbtrace = [];


for m = 1:length(mice)
    
    thisMouse = mice{m,1};
    mouseidx = find(strcmpi(thisMouse, dates.mouse));
    
    
    
    % get group acquisition curves
    startTraining = dates.prepairDate(mouseidx,1)+1;
    if strcmpi(thisMouse,'OK237') || strcmpi(thisMouse,'OK240')
        startTraining = 200106;
    end
    stopTraining = dates.lastacqDate(mouseidx,1);
    if ~isnan(startTraining)
        [trainingSessions] = getDateVector(startTraining, stopTraining);
        [paired]=updateCRRBProbStruct(trainingSessions, thisMouse, data, paired, timeVector, 'training');
    end
    
    % get group unpaired training curves
    startTraining = dates.preunpDate(mouseidx,1)+1;
    stopTraining = dates.lastunpDate(mouseidx,1);
    if ~isnan(startTraining)
        [trainingSessions] = getDateVector(startTraining, stopTraining);
        [unpaired]=updateCRRBProbStruct(trainingSessions, thisMouse, data, unpaired, timeVector, 'unpaired');
    end
    
    % get group extinction curves
    startTraining = dates.postpairDate(mouseidx,1)+1;
    stopTraining = dates.lastextDate(mouseidx,1);
    if ~isnan(startTraining)
        [trainingSessions] = getDateVector(startTraining, stopTraining);
        [extinction]=updateCRRBProbStruct(trainingSessions, thisMouse, data, extinction, timeVector, 'extinction');
    end
    
    % get rebound traces/amps/probabilities before training
    trainingSession = dates.prepairDate(mouseidx,1);
    if isnan(trainingSession)
        trainingSession = dates.postunpDate(mouseidx,1);
    end
    [rebound.prePaired]=updateRBStruct(trainingSession, thisMouse, data, rebound.prePaired, timeVector);
    
    % rebound features after training
    trainingSession = dates.postpairDate(mouseidx,1);
    if ~isnan(trainingSession)
        [rebound.postPaired]=updateRBStruct(trainingSession, thisMouse, data, rebound.postPaired, timeVector);
    end

    % rebound features before & after unpaired training
    trainingSession = dates.preunpDate(mouseidx,1);
    if ~isnan(trainingSession)
        [rebound.preUnp]=updateRBStruct(trainingSession, thisMouse, data, rebound.preUnp, timeVector);
    end
    trainingSession = dates.postunpDate(mouseidx,1);
    if ~isnan(trainingSession)
        [rebound.postUnp]=updateRBStruct(trainingSession, thisMouse, data, rebound.postUnp, timeVector);
    end
    
    % rebound features after extinction
    trainingSession = dates.postextDate(mouseidx,1);
    if ~isnan(trainingSession)
        [rebound.postExt]=updateRBStruct(trainingSession, thisMouse, data, rebound.postExt, timeVector);
    end

    
    clear thisMouse mouseidx midx eyelidpos sessdate trialtype csdur usdur...
        isi laserdur laserint laserdelay startTraining stopTraining
end


% find best laser intensity for rebounds for each mouse
anims = unique(rebound.postPaired.mouse);
rbprobs = nan(length(anims),3);
rbamps = nan(length(anims),3);
possibints = [15,30,60];
bestInt_prob = nan(length(anims),1);
bestInt_amp = nan(length(anims),1);
for a = 1:length(anims)
    aidx = find(strcmpi(rebound.postPaired.mouse,anims{a,1}));
    temp = rebound.postPaired.rbprob(aidx)';
    rbprobs(a,:) = temp(1:3);
    [~,maxidx] = max(temp);
    bestInt_prob(a,1) = possibints(maxidx);
    
    temp = rebound.postPaired.rbamp(aidx)';
    rbamps(a,:) = temp(1:3);
    [~,maxidx] = max(temp);
    bestInt_amp(a,1) = possibints(maxidx);
end

% make group figure for task acquisition
pairedphasecrprob = nan(20,2);
pairedphasecradjamp = nan(20,2);
unpairedphasecrprob = nan(20,2);
unpairedphasecradjamp = nan(20,2);
pairedphaseeyelidpos = nan(20,440);
unpairedphaseeyelidpos = nan(20,440);
for s = 1:20
    pairedphasecrprob(s,1) = nanmedian(paired.crprob(paired.session==s,1));
    pairedphasecrprob(s,2) = mad(paired.crprob(paired.session==s,1),1);
    
    unpairedphasecrprob(s,1) = nanmedian(unpaired.crprob(unpaired.session==s,1));
    unpairedphasecrprob(s,2) = mad(unpaired.crprob(unpaired.session==s,1),1);
    
    pairedphasecradjamp(s,1) = nanmedian(paired.cradjamp(paired.session==s,1));
    pairedphasecradjamp(s,2) = mad(paired.cradjamp(paired.session==s,1),1);
    
    unpairedphasecradjamp(s,1) = nanmedian(unpaired.cradjamp(unpaired.session==s,1));
    unpairedphasecradjamp(s,2) = mad(unpaired.cradjamp(unpaired.session==s,1),1);
    
    pairedphaseeyelidpos(s,:) = nanmean(paired.eyelidposadj(paired.session==s,:));
    unpairedphaseeyelidpos(s,:) = nanmean(unpaired.eyelidposadj(unpaired.session==s,:));
end
pairedmice = unique(paired.mouse);
lasteyelidtrace = nan(length(pairedmice),440);
for m = 1:length(pairedmice)
    sessions = unique(paired.session(strcmpi(paired.mouse,pairedmice{m,1})));
    lastsession = max(sessions);
    if lastsession>2
        lasteyelidtrace(m,:) = paired.eyelidposadj(...
            strcmpi(paired.mouse,pairedmice{m,1}) ...
            & paired.session==lastsession,:);
    end
end
unpairedmice = unique(unpaired.mouse);
unpairedlasteyelidtrace = nan(length(unpairedmice),440);
for m = 1:length(unpairedmice)
    sessions = unique(unpaired.session(strcmpi(unpaired.mouse,unpairedmice{m,1})));
    lastsession = max(sessions);
    unpairedlasteyelidtrace(m,:) = unpaired.eyelidposadj(...
        strcmpi(unpaired.mouse,unpairedmice{m,1}) ...
        & unpaired.session==lastsession,:);
end
figure
subplot(2,6,[1:6])
a=errorbar(1:length(pairedphasecrprob), pairedphasecrprob(:,1), pairedphasecrprob(:,2));
set(a,'CapSize',0, 'Marker', 'o')
hold on
b=errorbar(1:length(unpairedphasecrprob), unpairedphasecrprob(:,1), unpairedphasecrprob(:,2));
set(b,'CapSize',0, 'Marker', 'o')
ylabel('CR Probability')
legend('paired training', 'unpaired training', 'Location', 'NorthWest')
xlabel('session')
text(14, 0.4, ['paired n = ', num2str(length(unique(paired.mouse)))])
text(14, 0.3, ['unpaired n = ', num2str(length(unique(unpaired.mouse)))])
subplot(2,6,[7,8])
hold on
plot(timeVector, pairedphaseeyelidpos(1,:)')
plot(timeVector, nanmean(lasteyelidtrace))
plot(timeVector, unpairedphaseeyelidpos(1,:)')
plot(timeVector, nanmean(unpairedlasteyelidtrace))
xlim([0 0.3])
ylim([0 1])
legend('day 1, P', 'last day, P', 'day 1, U', 'last day, U', 'Location', 'NorthWest')
subplot(2,6,[9,10]) % for this one, just plot the RB's at the best intensity
rbtrace.prepaired = nan(length(anims),440);
rbtrace.postpaired = nan(length(anims),440);
rbtrace.preunpaired = nan(length(anims),440);
rbtrace.postunpaired = nan(length(anims),440);
for a = 1:length(anims)
    if ~strcmpi(anims{m,1},'OK235')
        [rbtrace.prepaired] = getTrace(rebound.prePaired, rbtrace.prepaired, a, bestInt_prob, anims{a,1});
        [rbtrace.postpaired] = getTrace(rebound.postPaired, rbtrace.postpaired, a, bestInt_prob, anims{a,1});
    end
    [rbtrace.preunpaired] = getTrace(rebound.preUnp, rbtrace.preunpaired, a, bestInt_prob, anims{a,1});
    [rbtrace.postunpaired] = getTrace(rebound.postUnp, rbtrace.postunpaired, a, bestInt_prob, anims{a,1});
end
hold on
plot(timeVector, nanmean(rbtrace.prepaired)')
plot(timeVector, nanmean(rbtrace.postpaired)')
plot(timeVector, nanmean(rbtrace.preunpaired)')
plot(timeVector, nanmean(rbtrace.postunpaired)')
xlim([0.85 1.4])
ylim([0 1])
subplot(2,6,11) % for this one, just plot the RB's at the best intensity 
rbprob.prepaired = nan(length(anims),1);
rbprob.postpaired = nan(length(anims),1);
rbprob.preunpaired = nan(length(anims),1);
rbprob.postunpaired = nan(length(anims),1);
for a = 1:length(anims)
    if ~strcmpi(anims{m,1},'OK235')
        [rbprob.prepaired] = getProb(rebound.prePaired, rbprob.prepaired, a, bestInt_prob, anims{a,1});
        [rbprob.postpaired] = getProb(rebound.postPaired, rbprob.postpaired, a, bestInt_prob, anims{a,1});
    end
    [rbprob.preunpaired] = getProb(rebound.preUnp, rbprob.preunpaired, a, bestInt_prob, anims{a,1});
    [rbprob.postunpaired] = getProb(rebound.postUnp, rbprob.postunpaired, a, bestInt_prob, anims{a,1});
end
hold on
plotMedianBoxplot(rbprob.prepaired, 1, 0.25, [0 0 1])
plotMedianBoxplot(rbprob.postpaired, 2, 0.25, [0 0 1])
plotMedianBoxplot(rbprob.preunpaired, 4, 0.25, [1 0 0])
plotMedianBoxplot(rbprob.postunpaired, 5, 0.25, [1 0 0])
text(4, 0.9, ['n = ',num2str(sum(~isnan(rbprob.prepaired)))], 'Color', [0 0 1])
text(4, 0.8, ['n = ',num2str(sum(~isnan(rbprob.preunpaired)))], 'Color', [1 0 0])
plot([1 1.9], [1 1], 'Color', [0 0 1])
text(1.5, 1.05, '**', 'HorizontalAlignment', 'center')
plot([2.1 5], [1 1], 'Color', [1 0 1])
text(3.5, 1.05, '**', 'HorizontalAlignment', 'center')
xlim([0.5 5.5])
ylim([0 1.1])
ylabel('Rebound Probability')
xticks([1,2,4,5])
xticklabels({'BP','AP','BU','BP'})
subplot(2,6,12) % for this one, just plot the RB's at the best intensity 
rbamp.prepaired = nan(length(anims),1);
rbamp.postpaired = nan(length(anims),1);
rbamp.preunpaired = nan(length(anims),1);
rbamp.postunpaired = nan(length(anims),1);
for a = 1:length(anims)
    [rbamp.prepaired] = getAmp(rebound.prePaired, rbamp.prepaired, a, bestInt_prob, anims{a,1});
    [rbamp.postpaired] = getAmp(rebound.postPaired, rbamp.postpaired, a, bestInt_prob, anims{a,1});
    [rbamp.preunpaired] = getAmp(rebound.preUnp, rbamp.preunpaired, a, bestInt_prob, anims{a,1});
    [rbamp.postunpaired] = getAmp(rebound.postUnp, rbamp.postunpaired, a, bestInt_prob, anims{a,1});
end
hold on
plotMedianBoxplot(rbamp.prepaired, 1, 0.25, [0 0 1])
plotMedianBoxplot(rbamp.postpaired, 2, 0.25, [0 0 1])
plotMedianBoxplot(rbamp.preunpaired, 4, 0.25, [1 0 0])
plotMedianBoxplot(rbamp.postunpaired, 5, 0.25, [1 0 0])
plot([1 1.9], [0.7 0.7], 'Color', [0 0 1])
text(1.5, 0.75, '**', 'HorizontalAlignment', 'center')
plot([2.1 5], [0.7 0.7], 'Color', [1 0 1])
text(3.5, 0.75, '**', 'HorizontalAlignment', 'center')
xlim([0.5 5.5])
ylim([0 1])
ylabel('Rebound Amplitude (FEC)')
xticks([1,2,4,5])
xticklabels({'BP','AP','BU','BP'})

%% write rebound probability and amplitude information to a table for R to look at
values = [rbprob.prepaired, rbprob.postpaired];
headers = {'PrePairProb','PostPairProb'};
tempcsv = [headers;num2cell(values)];
tempcsv = cell2table(tempcsv(2:end,:),'VariableNames',tempcsv(1,:));
writetable(tempcsv,'PreVsPostPairedTraining_RBProb_forPairedComp.csv')

values = [rbprob.preunpaired, rbprob.postunpaired];
headers = {'PreUnpairProb','PostUnpairProb'};
tempcsv = [headers;num2cell(values)];
tempcsv = cell2table(tempcsv(2:end,:),'VariableNames',tempcsv(1,:));
writetable(tempcsv,'PreVsPostUnpairedTraining_RBProb_forPairedComp.csv')

values = [rbprob.postpaired, rbprob.postunpaired];
headers = {'PostPairProb','PostUnpairProb'};
tempcsv = [headers;num2cell(values)];
tempcsv = cell2table(tempcsv(2:end,:),'VariableNames',tempcsv(1,:));
writetable(tempcsv,'PostPairVsUnpairedTraining_RBProb_forPairedComp.csv')

values = [rbprob.prepaired, rbprob.preunpaired];
headers = {'PrePairProb','PreUnpairProb'};
tempcsv = [headers;num2cell(values)];
tempcsv = cell2table(tempcsv(2:end,:),'VariableNames',tempcsv(1,:));
writetable(tempcsv,'PrePairVsUnpairedTraining_RBProb_forPairedComp.csv')

values = [rbamp.prepaired, rbamp.postpaired];
headers = {'PrePairAmp','PostPairAmp'};
tempcsv = [headers;num2cell(values)];
tempcsv = cell2table(tempcsv(2:end,:),'VariableNames',tempcsv(1,:));
writetable(tempcsv,'PreVsPostPairedTraining_RBAmp_forPairedComp.csv')

values = [rbamp.preunpaired, rbamp.postunpaired];
headers = {'PreUnpairAmp','PostUnpairAmp'};
tempcsv = [headers;num2cell(values)];
tempcsv = cell2table(tempcsv(2:end,:),'VariableNames',tempcsv(1,:));
writetable(tempcsv,'PreVsPostUnpairedTraining_RBAmp_forPairedComp.csv')

values = [rbamp.postpaired, rbamp.postunpaired];
headers = {'PostPairAmp','PostUnpairAmp'};
tempcsv = [headers;num2cell(values)];
tempcsv = cell2table(tempcsv(2:end,:),'VariableNames',tempcsv(1,:));
writetable(tempcsv,'PostPairVsUnpairedTraining_RBAmp_forPairedComp.csv')

values = [rbamp.prepaired, rbamp.preunpaired];
headers = {'PrePairAmp','PreUnpairAmp'};
tempcsv = [headers;num2cell(values)];
tempcsv = cell2table(tempcsv(2:end,:),'VariableNames',tempcsv(1,:));
writetable(tempcsv,'PrePairVsUnpairedTraining_RBAmp_forPairedComp.csv')

% subplot(2,4,[5,6])
% a=errorbar(1:length(pairedphasecradjamp), pairedphasecradjamp(:,1), pairedphasecradjamp(:,2));
% set(a,'CapSize',0, 'Marker', 'o')
% hold on
% b=errorbar(1:length(unpairedphasecradjamp), unpairedphasecradjamp(:,1), unpairedphasecradjamp(:,2));
% set(b,'CapSize',0, 'Marker', 'o')
% ylabel('CR Amplitude')
% legend('paired training', 'unpaired training', 'Location', 'NorthWest')
% xlabel('session')
% text(14, 0.3, ['paired n = ', num2str(length(unique(paired.mouse)))])
% text(14, 0.2, ['unpaired n = ', num2str(length(unique(unpaired.mouse)))])