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
[num, txt, raw] = xlsread('RBExptDates.xlsx');
xlsdata = raw(1:20,1:10);

dates.mouse = {};
dates.testInhibDate = nan(19,1); % 11 total subjects
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
%% stopped editing here
% fetch data
testInhibData = concatData_20200504(basedir, dates.mouse, dates.testInhibDate);
preUnpairedData = concatData_20200504(basedir, dates.mouse, dates.preunpDate);
lastUnpairedData = concatData_20200504(basedir, dates.mouse, dates.lastunpDate);
postUnpairedData = concatData_20200504(basedir, dates.mouse, dates.postunpDate);
prePairedData = concatData_20200504(basedir, dates.mouse, dates.prepairDate);
lastAcqData = concatData_20200504(basedir, dates.mouse, dates.lastacqDate);
postPairedData = concatData_20200504(basedir, dates.mouse, dates.postpairDate);
lastExtData = concatData_20200504(basedir, dates.mouse, dates.lastextDate);
postExtData = concatData_20200504(basedir, dates.mouse, dates.postextDate);

% save data
if strcmpi(machine, 'OREK')
    savedir = 'E:\pcp2ChR2 data\rebound';
elseif strcmpi(machine, 'COMPUPITAR')
    savedir = 'D:\pcp2ChR2 data\rebound';
end
cd(savedir)
save('overallData_200504.mat', 'testInhibData', 'preUnpairedData', ...
    'lastUnpairedData', 'postUnpairedData', 'prePairedData', 'lastAcqData',...
    'postPairedData', 'lastExtData', 'postExtData', 'dates')