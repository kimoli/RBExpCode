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

% % saving this part of the code because I think it will be useful for the
% % next level of processing, but for now it is limiting the main data mat in
% % a way that I don't like
% [num, txt, raw] = xlsread('RBExptDates.xlsx');
% xlsdata = raw(1:20,1:10);
% 
% dates.mouse = {};
% dates.testInhibDate = nan(19,1); % 11 total subjects
% dates.preunpDate = nan(19,1);
% dates.lastunpDate = nan(19,1);
% dates.postunpDate = nan(19,1);
% dates.prepairDate = nan(19,1); % same idx for posttraining, just 1 column over
% dates.lastacqDate = nan(19,1);
% dates.postpairDate = nan(19,1);
% dates.lastextDate = nan(19,1);
% dates.postextDate = nan(19,1); % the first group did not get extinction, remove 4
% for m = 2:20 % cycle through each mouse (each row of spreadsheet)
%     dates.mouse{m-1,1}=xlsdata{m,1};
%     for c = 2:10
%         switch xlsdata{1,c}
%             case 'Inhib Block CR'
%                 dates.testInhibDate(m-1,1)=cell2mat(xlsdata(m,c));
%             case 'Pre-Unpaired'
%                 dates.preunpDate(m-1,1)=cell2mat(xlsdata(m,c));
%             case 'Last Unp'
%                 dates.lastunpDate(m-1,1)=cell2mat(xlsdata(m,c));
%             case 'Post-Unpaired'
%                 dates.postunpDate(m-1,1)=cell2mat(xlsdata(m,c));
%             case 'Pre-Paired'
%                 dates.prepairDate(m-1,1)=cell2mat(xlsdata(m,c));
%             case 'Last Acq'
%                 dates.lastacqDate(m-1,1)=cell2mat(xlsdata(m,c));
%             case 'Post-Paired'
%                 dates.postpairDate(m-1,1)=cell2mat(xlsdata(m,c));
%             case 'Last Ext'
%                 dates.lastextDate(m-1,1)=cell2mat(xlsdata(m,c));
%             case 'Post-Extinction'
%                 dates.postextDate(m-1,1)=cell2mat(xlsdata(m,c));
%         end
%     end
% end
% %% stopped editing here
% % fetch data
% testInhibData = concatData_20200504(basedir, dates.mouse, dates.testInhibDate);
% preUnpairedData = concatData_20200504(basedir, dates.mouse, dates.preunpDate);
% lastUnpairedData = concatData_20200504(basedir, dates.mouse, dates.lastunpDate);
% postUnpairedData = concatData_20200504(basedir, dates.mouse, dates.postunpDate);
% prePairedData = concatData_20200504(basedir, dates.mouse, dates.prepairDate);
% lastAcqData = concatData_20200504(basedir, dates.mouse, dates.lastacqDate);
% postPairedData = concatData_20200504(basedir, dates.mouse, dates.postpairDate);
% lastExtData = concatData_20200504(basedir, dates.mouse, dates.lastextDate);
% postExtData = concatData_20200504(basedir, dates.mouse, dates.postextDate);

mice = {'OK213';'OK211';'OK214';'OK215';'OK216';'OK217';'OK218';'OK234';...
    'OK235';'OK236';'OK237';'OK238';'OK239';'OK240';'OK241'};
data.mouse = {};
data.eyelidpos = [];
data.date = [];
data.vidscorr = [];
data.type = [];
data.csdur = [];
data.usdur = [];
data.isi = [];
data.laserdur = [];
data.laserint = [];
data.laserdelay = [];
for m = 1:length(mice)
    mouseDir = [basedir,'\',mice{m,1}];
    cd(mouseDir)
    temp19 = dir('19*');
    temp20 = dir('20*');
    days = [temp19;temp20];
    clear temp19 temp20
    for d = 1:length(days)
        dayDir = [mouseDir,'\',days(d,1).name];
        cd(dayDir)
        if exist('newTrialdata.mat','file')==2
            vidscorr = 1;
            load('newTrialdata.mat')
            pulldata = 1;
        elseif exist('trialdata.mat','file')==2
            vidscorr = 0;
            load('trialdata.mat')
            pulldata = 1;
        else
            pulldata = 0;
        end
        
        if pulldata==1
            addmouse = cell(length(trials.c_csdur),1);
            [addmouse{1:end}] = deal(mice{m,1});
            addeyelidpos = nan(length(trials.c_csdur),440);
            for t = 1:length(trials.c_csdur)
                addeyelidpos(t,1:size(trials.eyelidpos,2)) = trials.eyelidpos(t,:);
            end
            
            
            data.mouse = [data.mouse;addmouse];
            data.eyelidpos = [data.eyelidpos;addeyelidpos];
            data.date = [data.date;ones(length(trials.c_csdur),1)*str2double(days(d,1).name)];
            data.vidscorr = [data.vidscorr;ones(length(trials.c_csdur),1)*vidscorr];
            data.type = [data.type;strcmpi(trials.type,'Conditioning')];
            data.csdur = [data.csdur;trials.c_csdur];
            data.usdur = [data.usdur; trials.c_usdur];
            data.isi = [data.isi; trials.c_isi];
            data.laserdur = [data.laserdur;trials.laser.dur];
            data.laserint = [data.laserint;trials.laser.amp];
            data.laserdelay = [data.laserdelay;trials.laser.delay];
            
            clear addmouse addeyelidpos
        end
    end
end

% save data
if strcmpi(machine, 'OREK')
    savedir = 'E:\pcp2ChR2 data\rebound';
elseif strcmpi(machine, 'COMPUPITAR')
    savedir = 'D:\pcp2ChR2 data\rebound';
end
cd(savedir)
save('RBExpt_allTrials.mat', 'data')