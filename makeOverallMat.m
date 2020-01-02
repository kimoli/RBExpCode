%%% make a big file with all the trialdata for the overall pre and post
%%% analyses

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
[num, txt, raw] = xlsread('prePostDays.xlsx');
xlsdata = raw(1:53,1:4);

mousename = {};
testInhibIdx = nan(11,1); % 11 total subjects
pretestIdx = nan(19,1); % same idx for posttraining, just 1 column over
lastacqIdx = nan(7,1);
postextIdx = nan(7,1); % the first group did not get extinction, remove 4
iter_testInhib = 1;
iter_pretest = 1;
iter_postext = 1;
iter_lastacq = 1;
for i = 2:53
    if strcmpi(xlsdata(i,2), 'test inhib')
        testInhibIdx(iter_testInhib,1) = i;
        mousename = [mousename; xlsdata{i,1}];
        iter_testInhib = iter_testInhib + 1;
    elseif strcmpi(xlsdata(i,2), 'RB post acquisition')
        pretestIdx(iter_pretest, 1) = i;
        iter_pretest = iter_pretest + 1;
    elseif strcmpi(xlsdata(i,2), 'RB post extinction')
        postextIdx(iter_postext,1) = i;
        iter_postext = iter_postext + 1;
    elseif strcmpi(xlsdata(i,2), 'last acquisition')
        lastacqIdx(iter_lastacq,1) = i;
        iter_lastacq = iter_lastacq + 1;
    else
        disp('bad experiment label in data')
        pause
    end
end
clear iter_testInhib iter_pretest iter_postext

% fetch data
testInhibData = concatData(basedir, xlsdata, 4, testInhibIdx); % posttests in col 4
pretestData = concatData(basedir, xlsdata, 3, pretestIdx); % pretest in col 3
posttestData = concatData(basedir, xlsdata, 4, pretestIdx);
postextData = concatData(basedir, xlsdata, 4, postextIdx);
lastacqData = concatData(basedir, xlsdata, 4, lastacqIdx);
lastextData = concatData(basedir, xlsdata, 3, postextIdx); % last extinction day is in column 3

% save data
if strcmpi(machine, 'OREK')
    savedir = 'E:\pcp2ChR2 data\rebound';
elseif strcmpi(machine, 'COMPUPITAR')
    savedir = 'D:\pcp2ChR2 data\rebound';
end
cd(savedir)
save('overallData_191219.mat', 'testInhibData', 'pretestData', 'posttestData', ...
    'postextData', 'lastacqData', 'lastextData')