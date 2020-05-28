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
paired.cradjampTm = [];
paired.rbamp = [];
paired.rbprob = [];
paired.rbtrace = [];
paired.rbtraceHit = [];
paired.eyelidposadj = [];
paired.eyelidposadjHit = [];

unpaired.mouse = {};
unpaired.session = [];
unpaired.crprob = [];
unpaired.cradjamp = [];
unpaired.cradjampTm = [];
unpaired.rbamp = [];
unpaired.rbprob = [];
unpaired.rbtrace = [];
unpaired.rbtraceHit = [];
unpaired.eyelidposadj = [];
unpaired.eyelidposadjHit = [];

extinction.mouse = {};
extinction.session = [];
extinction.crprob = [];
extinction.cradjamp = [];
extinction.cradjampTm = [];
extinction.rbamp = [];
extinction.rbprob = [];
extinction.rbtrace = [];
extinction.rbtraceHit = [];
extinction.eyelidposadj = [];
extinction.eyelidposadjHit = [];

savings.mouse = {};
savings.session = [];
savings.crprob = [];
savings.cradjamp = [];
savings.cradjampTm = [];
savings.rbamp = [];
savings.rbprob = [];
savings.rbtrace = [];
savings.rbtraceHit = [];
savings.eyelidposadj = [];
savings.eyelidposadjHit = [];

rebound.prePaired.mouse = {};
rebound.prePaired.laserint = [];
rebound.prePaired.rbamp = [];
rebound.prePaired.rbprob = [];
rebound.prePaired.rbtrace = [];
rebound.prePaired.rbtraceHit = [];

rebound.postPaired.mouse = {};
rebound.postPaired.laserint = [];
rebound.postPaired.rbamp = [];
rebound.postPaired.rbprob = [];
rebound.postPaired.rbtrace = [];
rebound.postPaired.rbtraceHit = [];

rebound.preUnp.mouse = {};
rebound.preUnp.laserint = [];
rebound.preUnp.rbamp = [];
rebound.preUnp.rbprob = [];
rebound.preUnp.rbtrace = [];
rebound.preUnp.rbtraceHit = [];

rebound.postUnp.mouse = {};
rebound.postUnp.laserint = [];
rebound.postUnp.rbamp = [];
rebound.postUnp.rbprob = [];
rebound.postUnp.rbtrace = [];
rebound.postUnp.rbtraceHit = [];

rebound.postExt.mouse = {};
rebound.postExt.laserint = [];
rebound.postExt.rbamp = [];
rebound.postExt.rbprob = [];
rebound.postExt.rbtrace = [];
rebound.postExt.rbtraceHit = [];


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
    checkdate = startTraining;
    monthRollover = 0;
    %pause
    if ~isnan(checkdate)
        while sum(data.date==checkdate & strcmpi(thisMouse, data.mouse))<150 % while not an extinction session
            if sum(data.date==checkdate)==0
                startstring = num2str(checkdate);
                monthnum = str2double(startstring(3:4));
                if monthnum==9 || monthnum==4 || monthnum==6 || monthnum==11
                    % month has 30 days
                    lastDayInMonth = [startstring(1:4),'30'];
                elseif monthnum==2
                    % month days depends on if year is a leap year (if it was 2019 or
                    % 2020 when I did the experiment, for simplification)
                    if strcmpi(startstring(1:2),'19')
                        lastDayInMonth = [startstring(1:4),'27'];
                    elseif strcmpi(startstring(1:2),'20')
                        lastDayInMonth = [startstring(1:4),'28'];
                    else
                        disp('script not designed to accommodate this year')
                        pause
                    end
                else
                    % month has 31 days
                    lastDayInMonth = [startstring(1:4),'31'];
                end
                if checkdate>=str2double(lastDayInMonth)
                    if monthnum==12 % need to increment year + month
                        pause
                        year = num2str(str2double(startstring(1:2))+1);
                        month = '01';
                        day = '01';
                        checkdate = str2double([year,month,day]);
                    else % need to increment month only
                        month = num2str(str2double(startstring(3:4))+1);
                        if length(month)==1
                            month = ['0',month];
                        end
                        day = '01';
                        checkdate = str2double([startstring(1:2),month,day]);
                    end
                    monthRollover = 2;
                else % just a missing training session
                    checkdate = checkdate + 1;
                    if monthRollover > 0
                        monthRollover = monthRollover -1;
                    end
                end
            else
                checkdate= checkdate+1;
                if monthRollover > 0
                    monthRollover = monthRollover -1;
                end
            end
        end
        if monthRollover>0
            startTraining = str2double(lastDayInMonth)-2+monthRollover;
        else
            startTraining = checkdate-2;
        end
%         sum(data.date==startTraining & strcmpi(thisMouse, data.mouse))
%         sum(data.date==checkdate-1 & strcmpi(thisMouse, data.mouse))
%         sum(data.date==checkdate & strcmpi(thisMouse, data.mouse))
%         pause
    end
    if ~isnan(startTraining)
        [trainingSessions] = getDateVector(startTraining, stopTraining);
        [extinction]=updateCRRBProbStruct(trainingSessions, thisMouse, data, extinction, timeVector, 'extinction');
    end
    
    % savings
    mousedates = unique(data.date(strcmpi(data.mouse,thisMouse),1));
    extposttestidx = find(mousedates==dates.postextDate(mouseidx,1));
    startTraining = mousedates(extposttestidx+1,1);
    stopTraining = mousedates(end,1);
%     dates.postextDate(mouseidx,1)
%     startTraining
%     stopTraining
%     pause
    if ~isnan(startTraining)
        [trainingSessions] = getDateVector(startTraining, stopTraining);
        [savings]=updateCRRBProbStruct(trainingSessions, thisMouse, data, savings, timeVector, 'savings');
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

%% make a version of the figure that would allow you to test the order effects in a normal ANOVA
% is this not going to work because I have  different number of days for
% acquisition and extinction/unpaired training?

% make group figure for task acquisition
PFMice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218'};
[pairedphasecrprob_PF, pairedphasecradjamp_PF, unpairedphasecrprob_PF,...
    unpairedphasecradjamp_PF, pairedphaseeyelidpos_PF, unpairedphaseeyelidpos_PF]=...
    getGroupOfMiceData(PFMice, paired, extinction);
UFMice = {'OK234';'OK236';'OK237';'OK238';'OK239';'OK240';'OK241'};
[pairedphasecrprob_UF, pairedphasecradjamp_UF, unpairedphasecrprob_UF,...
    unpairedphasecradjamp_UF, pairedphaseeyelidpos_UF, unpairedphaseeyelidpos_UF]=...
    getGroupOfMiceData(UFMice, paired, unpaired);

%% make extinction figure
extnctionmice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218';...
    'OK234';'OK236';'OK237';'OK238';'OK239';'OK240';'OK241'};
[extphasecrprob, extphasecradjamp, savphasecrprob,...
    savphasecradjamp, extphaseeyelidpos, savphaseeyelidpos]=...
    getGroupOfMiceData(extnctionmice, extinction, savings);
figure
subplot(3,6,[1:6])
b=errorbar(1:length(extphasecrprob), extphasecrprob(:,1), extphasecrprob(:,2));
set(b,'CapSize',0, 'Marker', 'o', 'LineStyle', 'none')
hold on
c=errorbar(13:12+length(savphasecrprob), savphasecrprob(:,1), savphasecrprob(:,2));
set(c,'CapSize',0, 'Marker', 'o', 'LineStyle', 'none')
plot([2.5 2.5], [0 1], 'Color', [0 0 0], 'LineStyle', '--')
plot([12.5 12.5], [0 1], 'Color', [0 0 0], 'LineStyle', '--')
xlim([0.5 17.5])
ylabel('CR Probability')
xlabel('session')
text(10, 0.9, ['n = ', num2str(length(unique(unpaired.mouse)))])
subplot(3,6,[7,8])
hold on
for i = 2:12
plot(timeVector, extphaseeyelidpos(i,:)')
end
xlim([0 0.3])
ylim([0 1])
legend('B','1','2','3','4','5','6','7','8','9','10', 'Location', 'NorthWest')
subplot(3,6,[9,10]) % for this one, just plot the RB's at the best intensity
rbtrace.preunpaired = nan(length(anims),440);
rbtrace.postunpaired = nan(length(anims),440);
for a = 1:length(anims)
    [rbtrace.preunpaired] = getTrace(rebound.postPaired, rbtrace.preunpaired, a, bestInt_prob, anims{a,1});
    [rbtrace.postunpaired] = getTrace(rebound.postExt, rbtrace.postunpaired, a, bestInt_prob, anims{a,1});
end
hold on
plot(timeVector, nanmean(rbtrace.preunpaired)')
plot(timeVector, nanmean(rbtrace.postunpaired)')
xlim([0 1.4])
ylim([0 1])
subplot(3,6,11) % for this one, just plot the RB's at the best intensity 
rbprob.preunpaired = nan(length(anims),1);
rbprob.postunpaired = nan(length(anims),1);
for a = 1:length(anims)
    [rbprob.preunpaired] = getProb(rebound.postPaired, rbprob.preunpaired, a, bestInt_prob, anims{a,1});
    [rbprob.postunpaired] = getProb(rebound.postExt, rbprob.postunpaired, a, bestInt_prob, anims{a,1});
end
hold on
plotMedianBoxplot(rbprob.preunpaired, 1, 0.25, [1 0 0])
plotMedianBoxplot(rbprob.postunpaired, 2, 0.25, [1 0 0])
xlim([0.5 2.5])
ylim([0 1.1])
ylabel('Rebound Probability')
subplot(3,6,12) % for this one, just plot the RB's at the best intensity 
rbamp.preunpaired = nan(length(anims),1);
rbamp.postunpaired = nan(length(anims),1);
for a = 1:length(anims)
    [rbamp.preunpaired] = getAmp(rebound.postPaired, rbamp.preunpaired, a, bestInt_prob, anims{a,1});
    [rbamp.postunpaired] = getAmp(rebound.postExt, rbamp.postunpaired, a, bestInt_prob, anims{a,1});
end
hold on
plotMedianBoxplot(rbamp.preunpaired, 1, 0.25, [1 0 0])
plotMedianBoxplot(rbamp.postunpaired, 2, 0.25, [1 0 0])
xlim([0.5 2.5])
ylim([0 1])
ylabel('Rebound Amplitude (FEC)')
sessionsToShow10Pct_acq = nan(length(extnctionmice), 1);
sessionsToShow10Pct_sav = nan(length(extnctionmice), 1);
reachThreshEyepos = nan(length(extnctionmice),440);
reachThreshEyepos_Sav = nan(length(extnctionmice),440);
lasteyelidtraceHit = nan(length(extnctionmice),440);
thresh = [0.6, 0.8, 0.6, 0.6, 0.6, 0.55, 0.6, 0.6, 0.7, 0.5, 0.6, 0.5, 0.6, 0.7];
for i=1:length(extnctionmice)    
    temp = paired.cradjamp(strcmpi(paired.mouse, extnctionmice{i,1}),:)>=thresh(i);
    temp = find(temp);
    if i==7 || i==10 || i==12 % big startles on first day of training
        sessionsToShow10Pct_acq(i,1) = temp(2);
        
        reachThreshEyepos(i,:) = paired.eyelidposadjHit(temp(2),:);
    else
        sessionsToShow10Pct_acq(i,1) = temp(1);
        
        reachThreshEyepos(i,:) = paired.eyelidposadjHit(temp(1),:);
    end
    temp = strcmpi(paired.mouse, extnctionmice{i,1});
    temp = find(temp);
    lasteyelidtraceHit(i,:) = paired.eyelidposadjHit(temp(end),:);
    
    temp = savings.cradjamp(strcmpi(savings.mouse, extnctionmice{i,1}),:)>=thresh(i);
    temp = find(temp);
    sessionsToShow10Pct_sav(i,1) = temp(1);
    reachThreshEyepos_Sav(i,:) = savings.eyelidposadjHit(temp(1),:);
end
savsessions = sessionsToShow10Pct_acq-sessionsToShow10Pct_sav;
subplot(3,6,[13,14])
a = scatter(sessionsToShow10Pct_sav, rbprob.postunpaired, 4);
set(a, 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1])
hold on
lsline
[r,p]=corr(sessionsToShow10Pct_sav, rbprob.postunpaired, 'Type', 'Spearman');
text(7, 0.8, ['r =',num2str(r)])
text(7, 0.7, ['p=',num2str(p)])
xlim([0.5 10.5])
ylim([0 1])
ylabel('RB Prob after Ext')
xlabel('Sessions to Reacquire')
subplot(3,6,15)
hold on
plotMedianBoxplot(sessionsToShow10Pct_acq, 1, 0.25, [1 0 0])
plotMedianBoxplot(sessionsToShow10Pct_sav, 2, 0.25, [1 0 0])
ylabel('Sessions')
xticks([1,2])
xticklabels({'Acq','Sav'})
xlim([0.5 2.5])
subplot(3,6,[16 17])
plot(timeVector, nanmean(reachThreshEyepos)')
hold on
plot(timeVector, nanmean(reachThreshEyepos_Sav)')
plot(timeVector, nanmean(lasteyelidtraceHit))
xlim([0 0.3])
legend('Reacquired', 'Last Training', 'Location', 'NorthWest')


sessionsToShow10Pct_ext = nan(length(extnctionmice), 1);
for i=1:length(extnctionmice)
    temp = extinction.crprob(strcmpi(extinction.mouse, extnctionmice{i,1}),:)<=0.15;
    temp = find(temp);
    sessionsToShow10Pct_ext(i,1) = temp(1);
end
figure
scatter(sessionsToShow10Pct_ext, rbprob.postpaired)
hold on
lsline
xlim([0 11])
ylim([0 1])
[r,p]= corr(sessionsToShow10Pct_ext, rbamp.postpaired, 'Type', 'Spearman');
text(2, 0.8, ['r =',num2str(r)])
text(2, 0.7, ['p=',num2str(p)])
xlabel('Sessions to Extinction')
ylabel('Rebound Probability after Training')
pause
figure
hold on
plot(timeVector, extphaseeyelidpos(12,:))
for s = 1:5
    plot(timeVector, savphaseeyelidpos(s,:))
end
xlim([0 0.3])
ylim([0 1])
legend('B','1','2','3','4','5')

disp('GOT TO PART WANT TO PAY ATTENTION TO')
pause

%% figure illustrating relationship between task acquisition and rebound acquisition
rbprobs_withinTraining = nan(20,1);
rbtraces_withinTraining = nan(20,440);
pairedphaseeyelidposThruTrain = nan(20,440);
crprobs_withinTraining = nan(20,1);
rbamps_withinTraining = nan(20,1);
cramps_withinTraining = nan(20,1);
rbprobs_extinction = nan(20,1);
rbamps_extinction = nan(20,1);
rbtraces_extinction = nan(20,440);
pairedphaseeyelidposThruExt = nan(20,440);
crprobs_ext = nan(20,1);
cramps_ext = nan(20,1);
outsem = [];
outsemC = [];
for s = 1:20
    idx = paired.session==s;
    rbprobs_withinTraining(s,1) = nanmean(paired.rbprob(idx));
    rbamps_withinTraining(s,1) = nanmean(paired.rba0000000mp(idx));
    temp = paired.rbtraceHit(idx,:);
    rbtraces_withinTraining(s,:) = nanmean(temp(~isnan(paired.rbprob(idx)),:));
    temp = paired.eyelidposadj(idx,:);
    pairedphaseeyelidposThruTrain(s,:) = nanmean(temp(~isnan(paired.rbprob(idx)),:));
    temp = paired.crprob(idx,:);
    crprobs_withinTraining(s,:) = nanmean(temp(~isnan(paired.rbprob(idx)),:));
    temp = paired.cradjampTm(idx,:);
    cramps_withinTraining(s,:) = nanmean(temp(~isnan(paired.rbprob(idx)),:));
        
    idx = extinction.session==s;
    rbprobs_extinction(s,1) = nanmean(extinction.rbprob(idx));
    rbamps_extinction(s,1) = nanmean(extinction.rbamp(idx));
    outsem(s,1) = nanstd(extinction.rbamp(idx))./sqrt(sum(~isnan(extinction.rbprob(idx))));
    temp = extinction.rbtraceHit(idx,:);
    rbtraces_extinction(s,:) = nanmean(temp(~isnan(extinction.rbprob(idx)),:));
    temp = extinction.eyelidposadj(idx,:);
    pairedphaseeyelidposThruExt(s,:) = nanmean(temp(~isnan(extinction.rbprob(idx)),:));
    temp = extinction.crprob(idx,:);
    crprobs_ext(s,:) = nanmedian(temp(~isnan(extinction.rbprob(idx)),:));
    temp = extinction.cradjampTm(idx,:);
    cramps_ext(s,:) = nanmedian(temp(~isnan(extinction.rbprob(idx)),:));
    outsemC(s,1) = nanstd(temp(~isnan(extinction.rbprob(idx))))./sqrt(sum(~isnan(extinction.rbprob(idx))));
end
[r,p]=corr(rbprobs_withinTraining(1:19), crprobs_withinTraining(1:19,1),'Type','Spearman'); % is the right way to do these comparisons to limit the number of days?
[r,p]=corr(rbamps_withinTraining(1:19), cramps_withinTraining(1:19,1),'Type','Spearman');
[r,p]=corr(rbprobs_extinction(2:12), crprobs_ext(2:12,1),'Type','Spearman');
[r,p]=corr(rbamps_extinction(2:12), cramps_ext(2:12,1),'Type','Spearman');
figure
subplot(2,4,1)
a = scatter(crprobs_withinTraining(1:19,1),rbprobs_withinTraining(1:19),  4);
set(a, 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1])
hold on
lsline
xlim([0 1])
ylim([0 1])
xlabel('CR Probability')
ylabel('Rebound Probability')
title('acquisition')
subplot(2,4,2)
a = scatter(rbamps_withinTraining(1:19), cramps_withinTraining(1:19,1), 4);
set(a, 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1])
hold on
lsline
xlim([0 1])
ylim([0 1])
xlabel('CR Amplitude')
ylabel('Rebound Amplitude')
title('acquisition')
subplot(2,4,3)
daysToPlot = [1;5;10;12;17;18;19];
hold on
for d = 1:length(daysToPlot)
    plot(timeVector, pairedphaseeyelidposThruTrain(daysToPlot(d),:))
end
xlim([0 0.55])
ylim([-0.025 1])
xlabel('Time from Tone (s)')
ylabel('Eyelid Position (FEC)')
subplot(2,4,4)
daysToPlot = [1;5;10;12;17;18;19];
hold on
for d = 1:length(daysToPlot)
    plot(timeVector-0.85, rbtraces_withinTraining(daysToPlot(d),:))
end
xlim([0 0.55])
ylim([-0.025 1])
xlabel('Time from Tone (s)')
ylabel('Eyelid position')
subplot(2,4,5)
a = scatter(crprobs_ext(2:12,1),rbprobs_extinction(2:12),   4);
set(a, 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1])
hold on
lsline
xlim([0 1])
ylim([0 1])
xlabel('CR Probability')
ylabel('Rebound Probability')
title('extinction')
subplot(2,4,6)
a = scatter(cramps_ext(2:12,1), rbamps_extinction(2:12),  4);
set(a, 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1])
hold on
lsline
xlim([0 1])
ylim([0 1])
xlabel('CR Amplitude')
ylabel('Rebound Amplitude')
title('extinction')
subplot(2,4,7)
daysToPlot = [1;3;4;5;8;9;12];
hold on
for d = 1:length(daysToPlot)
    plot(timeVector, pairedphaseeyelidposThruExt(daysToPlot(d),:))
end
xlim([0 0.55])
ylim([-0.025 1])
xlabel('Time from Tone (s)')
ylabel('Eyelid Position (FEC)')
subplot(2,4,8)
daysToPlot = [1;3;4;5;8;9;12];
hold on
for d = 1:length(daysToPlot)
    plot(timeVector-0.85, rbtraces_extinction(daysToPlot(d),:))
end
xlim([0 0.55])
ylim([-0.025 1])
xlabel('Time from Tone (s)')
ylabel('Eyelid position')

colordef black
figure
hold on
e = errorbar(1:11, rbamps_extinction(2:12), outsem(2:12), 'LineStyle', 'none', 'Color', [0 1 1]);
set(e, 'CapSize', 0)
a = scatter(1:11, rbamps_extinction(2:12), 20, 'MarkerFaceColor', [0 1 1], 'MarkerEdgeColor', [0 1 1]);
hold on
b = scatter(1:11, cramps_ext(2:12,1),  20, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [1 1 1]);
e = errorbar(1:11, cramps_ext(2:12), outsemC(2:12), 'LineStyle', 'none', 'Color', [1 1 1]);
set(e, 'CapSize', 0)
ylabel('Amplitude')
xlabel('Sessions')
legend([a, b], 'RB Amp', 'CR Amp', 'Location', 'NorthOutside')
ylim([0 0.6])


%% plot data for the inhibition test trials

inhibData.CSUStrace = nan(7,440);
inhibData.CSlaserTrace = nan(7,440);
for m = 1:length(dates.mouse)
    inhibSession = dates.testInhibDate(m,1);
    thisMouse = dates.mouse{m,1};
    idx = find(strcmpi(data.mouse,thisMouse) & data.date==inhibSession & data.type==1);
    
    eyelidpos = data.eyelidpos(idx,:);
    csdur = data.csdur(idx,:);
    usdur = data.usdur(idx,:);
    laserdur = data.laserdur(idx,:);
    laserdel = data.laserdelay(idx,:);
    
    laseridx = laserdur>0;
    trainingidx = laserdur==0;
    
    figure
    plot(nanmean(eyelidpos(laseridx,:))')
    hold on
    plot(nanmean(eyelidpos(trainingidx,:))')
    xlim([40 90])
    pause

    if sum(laseridx)>0
    end
end
% I think I want different data for this

% find sessions where I tried to block the CR, 160 ms latency laser and 35
% ms duration
inhibTestDirs = {};
basedir = 'E:\pcp2ChR2 data\rebound';
cd(basedir)
mice = dir('OK*');
for m = 1:length(mice)
    goHere = [basedir, '\', mice(m,1).name];
    cd(goHere)
    days = dir('19*');
    for d = 1:length(days)
        daydir = [goHere, '\', days(d,1).name];
        cd(daydir)
        if exist('newtrialdata.mat','file')==2
            load('newtrialdata.mat')
            checkday = 1;
        elseif exist('trialdata.mat','file')==2
            load('trialdata.mat')
            checkday = 1;
        else
            checkday = 0;
        end
        if checkday==1
            laserPlusCS = trials.laser.dur>0 & trials.c_csdur>0 & trials.c_usdur>0;
            if sum(laserPlusCS)>0
                trials.laser.delay(laserPlusCS)
                pause
                inhibTestDirs = [inhibTestDirs; daydir];
                break
            end
        end
    end
end

inhibCR.mouse = cell(4,1);
inhibCR.inhibTrace = nan(4,200);
inhibCR.baselineTrace = nan(4,200);
inhibCR.inhibCRAmp = nan(4,1);
inhibCR.baselineCRAmp = nan(4,1);
inhibCR.inhibCRProb = nan(4,1);
inhibCR.baselineCRProb = nan(4,1);

for i = 1:length(inhibTestDirs)
    cd(inhibTestDirs{i,1})
    load('trialdata.mat')
    laserPlusCS = find(trials.laser.dur>0 & trials.c_csdur>0 & trials.laser.delay==160);
    baselineTrials = find(trials.laser.dur==0 & trials.c_csdur>0 & trials.c_usdur>0);
    
    inhibCR.mouse{i,1} = inhibTestDirs{i,1}(29:33);
    [inhibCR.inhibCRAmp(i,1), inhibCR.inhibCRProb(i,1), inhibCR.inhibTrace(i,:)]=getCRProbAdjampTrace(trials, laserPlusCS);
    
%     day = str2double(inhibTestDirs{i,1}(end-5:end));
%     prevday = day-1;
%     goTo = strcat(inhibTestDirs{i,1}(1:end-6), num2str(prevday));
%     cd(goTo)
%     load('trialdata.mat')
%    baselineTrials = find(trials.laser.dur==0 & trials.c_csdur>0 & trials.c_usdur>0);
    [inhibCR.baselineCRAmp(i,1), inhibCR.baselineCRProb(i,1), inhibCR.baselineTrace(i,:)]=getCRProbAdjampTrace(trials, baselineTrials);
end

figure
subplot(1,3,[1,2])
a=shadedErrorBar(timeVector(1:200), nanmean(inhibCR.baselineTrace)', nanstd(inhibCR.baselineTrace)./sqrt(4), '-r', 1);
hold on
b=shadedErrorBar(timeVector(1:200), nanmean(inhibCR.inhibTrace)', nanstd(inhibCR.inhibTrace)./sqrt(4), '-b', 1);
plot([0.16 0.16], [0 1], 'Color', [0 0 0], 'LineStyle', '--')
plot([0.16+0.035 0.16+0.035], [0 1], 'Color', [0 0 0], 'LineStyle', '--')
xlim([0 0.28])
ylim([0 1])
ylabel('Eyelid Position (FEC)')
xlabel('time from tone (s)')
legend([a.mainLine, b.mainLine], 'CS + US', 'CS + laser + US', 'Location', 'NorthWest')
text(0.025, 0.7, 'n = 4')
subplot(1,3,3)
hold on
plotMedianBoxplot(inhibCR.baselineCRAmp, 1, 0.25, [1 0 0])
plotMedianBoxplot(inhibCR.inhibCRAmp, 2, 0.25, [0 0 1])
ylim([0 1])
xlim([0.5 2.5])
ylabel('Eyelid Closure Before Puff')
xticks([1 2])
xticklabels({'CS + US', 'CS + laser + US'})
text(0.75, 0.75, 't(3) = -8.22, p = 0.004')

% for right this minute, just use the data from the first mouse?

%% Below here is saved junk for exploratory savings analysis
% how many days until mouse showed > 10% CRs in acquisition
sessionsToShow10Pct_acq = nan(length(extnctionmice), 1);
sessionsToShow10Pct_sav = nan(length(extnctionmice), 1);
for i=1:length(extnctionmice)
    temp = paired.crprob(strcmpi(paired.mouse, extnctionmice{i,1}),:)>=0.1;
    temp = find(temp);
    sessionsToShow10Pct_acq(i,1) = temp(1);
    
    
    temp = savings.crprob(strcmpi(savings.mouse, extnctionmice{i,1}),:)>=0.1;
    temp = find(temp);
    sessionsToShow10Pct_sav(i,1) = temp(1);
end
[sessionsToShow10Pct_acq, sessionsToShow10Pct_sav]


sessionsToShow10Pct_acq = nan(length(extnctionmice), 1);
sessionsToShow10Pct_sav = nan(length(extnctionmice), 1);
for i=1:length(extnctionmice)
    temp = paired.crprob(strcmpi(paired.mouse, extnctionmice{i,1}),:)>=0.5;
    temp = find(temp);
    sessionsToShow10Pct_acq(i,1) = temp(1);
    
    
    temp = savings.crprob(strcmpi(savings.mouse, extnctionmice{i,1}),:)>=0.5;
    temp = find(temp);
    sessionsToShow10Pct_sav(i,1) = temp(1);
end
hist(sessionsToShow10Pct_acq-sessionsToShow10Pct_sav)

sessionsToShow10Pct_acq = nan(length(extnctionmice), 1);
sessionsToShow10Pct_sav = nan(length(extnctionmice), 1);
for i=1:length(extnctionmice)
    temp = paired.crprob(strcmpi(paired.mouse, extnctionmice{i,1}),:)>=0.65;
    temp = find(temp);
    sessionsToShow10Pct_acq(i,1) = temp(1);
    
    
    temp = savings.crprob(strcmpi(savings.mouse, extnctionmice{i,1}),:)>=0.65;
    temp = find(temp);
    sessionsToShow10Pct_sav(i,1) = temp(1);
end
figure
hist(sessionsToShow10Pct_acq-sessionsToShow10Pct_sav)
savsessions = sessionsToShow10Pct_acq-sessionsToShow10Pct_sav;

% this one works but it would be sort of weird to write it up
figure
scatter(savsessions, rbprob.preunpaired)
[r,p]=corr(savsessions, rbprob.preunpaired, 'Type', 'Spearman')

figure
scatter(savsessions, rbprob.postunpaired)
[r,p]=corr(savsessions, rbprob.postunpaired, 'Type', 'Spearman')


figure
scatter(sessionsToShow10Pct_sav, rbprob.postunpaired)
[r,p]=corr(sessionsToShow10Pct_sav, rbprob.postunpaired, 'Type', 'Kendall')
[r,p]=corr(sessionsToShow10Pct_sav, rbprob.postunpaired, 'Type', 'Spearman')
[r,p]=corr(log(sessionsToShow10Pct_sav), log(rbprob.postunpaired), 'Type', 'Spearman')



figure
scatter(savsessions, rbprob.postunpaired)
[r,p]=corr(sessionsToShow10Pct_acq-sessionsToShow10Pct_sav, rbprob.postunpaired, 'Type', 'Spearman')
[r,p]=corr(log(sessionsToShow10Pct_acq-sessionsToShow10Pct_sav), log(rbprob.postunpaired), 'Type', 'Spearman')

% try excluding the 3 mice with no rebounds after extinction
figure
scatter(savsessions, rbprob.postunpaired)
hold on
scatter(savsessions(rbprob.postunpaired>0), rbprob.postunpaired(rbprob.postunpaired>0))
[r,p]=corr(savsessions(rbprob.postunpaired>0), rbprob.postunpaired(rbprob.postunpaired>0), 'Type', 'Spearman')

[rbprob.postunpaired(rbprob.postunpaired==0) rbprob.postpaired(rbprob.postunpaired==0)]


% try excluding the 2 mice with a lot of savings and no rebounds
figure
scatter(savsessions, rbprob.postunpaired)
hold on
scatter(savsessions(savsessions~=8), rbprob.postunpaired(savsessions~=8))
[r,p]=corr(savsessions(savsessions~=8), rbprob.postunpaired(savsessions~=8), 'Type', 'Spearman')


figure
scatter(savsessions, rbamp.postunpaired)
[r,p]=corr(sessionsToShow10Pct_acq-sessionsToShow10Pct_sav, rbamp.postunpaired, 'Type', 'Spearman')
[r,p]=corr(log(sessionsToShow10Pct_acq-sessionsToShow10Pct_sav), log(rbamp.postunpaired), 'Type', 'Spearman')

figure
scatter(sessionsToShow10Pct_sav, rbprob.postunpaired)
[r,p]=corr(sessionsToShow10Pct_sav, rbprob.postunpaired, 'Type', 'Spearman')
[r,p]=corr(log(sessionsToShow10Pct_sav), log(rbprob.postunpaired), 'Type', 'Spearman')

figure
scatter(sessionsToShow10Pct_sav, rbamp.postunpaired)
[r,p]=corr(sessionsToShow10Pct_sav, rbamp.postunpaired, 'Type', 'Spearman')
[r,p]=corr(log(sessionsToShow10Pct_sav), log(rbamp.postunpaired), 'Type', 'Spearman')


figure
scatter(sessionsToShow10Pct_sav(sessionsToShow10Pct_sav<10), rbprob.postunpaired(sessionsToShow10Pct_sav<10)-rbprob.preunpaired(sessionsToShow10Pct_sav<10))
[r,p]=corr(sessionsToShow10Pct_sav(sessionsToShow10Pct_sav<10), rbprob.postunpaired(sessionsToShow10Pct_sav<10)-rbprob.preunpaired(sessionsToShow10Pct_sav<10), 'Type', 'Spearman')

figure
scatter(savsessions, rbprob.preunpaired)
[r,p]=corr(savsessions, rbprob.preunpaired, 'Type', 'Spearman')


figure
scatter(savsessions(sessionsToShow10Pct_sav<10), rbamp.postunpaired(sessionsToShow10Pct_sav<10)-rbamp.preunpaired(sessionsToShow10Pct_sav<10))
[r,p]=corr(savsessions(sessionsToShow10Pct_sav<10), rbamp.postunpaired(sessionsToShow10Pct_sav<10)-rbamp.preunpaired(sessionsToShow10Pct_sav<10), 'Type', 'Spearman')


%% try looking at savings in terms of CR Amplitude

sessionsToShow10Pct_acq = nan(length(extnctionmice), 1);
sessionsToShow10Pct_sav = nan(length(extnctionmice), 1);
thresh = [0.6, 0.8, 0.6, 0.6, 0.6, 0.55, 0.6, 0.6, 0.7, 0.5, 0.6, 0.5, 0.6, 0.7];
for i=1:length(extnctionmice)
%     figure
%     plot(paired.cradjamp(strcmpi(paired.mouse, extnctionmice{i,1}),:))
%     hold on
%     plot(savings.cradjamp(strcmpi(savings.mouse, extnctionmice{i,1}),:))
%     pause
%     close all
   
    
    temp = paired.cradjamp(strcmpi(paired.mouse, extnctionmice{i,1}),:)>=thresh(i);
    temp = find(temp);
    if i==10 || i==12 % big startles on first day of training
        sessionsToShow10Pct_acq(i,1) = temp(2);
    else
        sessionsToShow10Pct_acq(i,1) = temp(1);
    end
    
    
    temp = savings.cradjamp(strcmpi(savings.mouse, extnctionmice{i,1}),:)>=thresh(i);
    temp = find(temp);
    sessionsToShow10Pct_sav(i,1) = temp(1);
end
savsessions = sessionsToShow10Pct_acq-sessionsToShow10Pct_sav;
figure
scatter(sessionsToShow10Pct_sav, rbprob.postunpaired)
[r,p]=corr(sessionsToShow10Pct_sav, rbprob.postunpaired, 'Type', 'Spearman')
% THE ONE ABOVE WORKS AND IT WOULD BE EASY TO WRITE UP

% this one works but it would be sort of weird to write it up
figure
scatter(savsessions, rbprob.preunpaired)
[r,p]=corr(savsessions, rbprob.preunpaired, 'Type', 'Spearman')

figure
scatter(savsessions, rbprob.postunpaired)
[r,p]=corr(savsessions, rbprob.postunpaired, 'Type', 'Spearman')

figure
scatter(savsessions, rbamp.preunpaired)
[r,p]=corr(savsessions, rbamp.preunpaired, 'Type', 'Spearman')

figure
scatter(savsessions, rbamp.postunpaired)
[r,p]=corr(savsessions, rbamp.postunpaired, 'Type', 'Spearman')

figure
scatter(sessionsToShow10Pct_sav, rbamp.postunpaired)
[r,p]=corr(sessionsToShow10Pct_sav, rbamp.postunpaired, 'Type', 'Spearman')
[r,p]=corr(sessionsToShow10Pct_sav, rbamp.preunpaired, 'Type', 'Spearman')
[r,p]=corr(sessionsToShow10Pct_sav, rbprob.postunpaired, 'Type', 'Spearman')
[r,p]=corr(sessionsToShow10Pct_sav, rbprob.preunpaired, 'Type', 'Spearman')

% 
% sessionsToShow10Pct_acq = nan(length(extnctionmice), 1);
% sessionsToShow10Pct_sav = nan(length(extnctionmice), 1);
% thresh = [ones(1,7)*0.7,0.65,0.7, 0.65, 0.65, 0.7, 0.65, 0.7];
% for i=1:length(extnctionmice)
%     temp = paired.crprob(strcmpi(paired.mouse, extnctionmice{i,1}),:)>=thresh(i);
%     temp = find(temp);
%     sessionsToShow10Pct_acq(i,1) = temp(1);
%     
%     
%     temp = savings.crprob(strcmpi(savings.mouse, extnctionmice{i,1}),:)>=thresh(i);
%     temp = find(temp);
%     sessionsToShow10Pct_sav(i,1) = temp(1);
% end
% figure
% hist(sessionsToShow10Pct_acq-sessionsToShow10Pct_sav)



% try only looking at relationships at 30 mW
rbprob.preunpaired = nan(length(anims),1);
rbprob.postunpaired = nan(length(anims),1);
for a = 1:length(anims)
    [rbprob.preunpaired] = getProb(rebound.postPaired, rbprob.preunpaired, a, ones(length(bestInt_prob),1)*60, anims{a,1});
    [rbprob.postunpaired] = getProb(rebound.postExt, rbprob.postunpaired, a, ones(length(bestInt_prob),1)*60, anims{a,1});
end
rbamp.preunpaired = nan(length(anims),1);
rbamp.postunpaired = nan(length(anims),1);
for a = 1:length(anims)
    [rbamp.preunpaired] = getAmp(rebound.postPaired, rbamp.preunpaired, a, ones(length(bestInt_prob),1)*60, anims{a,1});
    [rbamp.postunpaired] = getAmp(rebound.postExt, rbamp.postunpaired, a, ones(length(bestInt_prob),1)*60, anims{a,1});
end

figure
savsessions = sessionsToShow10Pct_acq-sessionsToShow10Pct_sav;
scatter(savsessions, rbprob.postunpaired)
hold on
scatter(savsessions(savsessions~=8), rbprob.postunpaired(savsessions~=8))
[r,p]=corr(savsessions(savsessions~=8), rbprob.postunpaired(savsessions~=8), 'Type', 'Spearman')


% TO DO add savings analysis
% TO DO add field for RB hit trials only traces/amplitudes

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