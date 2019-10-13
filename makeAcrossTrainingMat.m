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

mice = {'OK214'; 'OK216'; 'OK217'};

xtraindata.eyelidpos = [];
xtraindata.mouse = [];
xtraindata.day = [];
xtraindata.csdur = [];
xtraindata.lasdur = [];
xtraindata.usdur = [];
xtraindata.lasdel = [];
xtraindata.lasamp = [];

for m = 1:length(mice)
    mouseFolder = [basedir, '\', mice{m,1}];
    cd(mouseFolder)
    days = dir('19*');
    for d = 1:length(days)
        dayFolder = [mouseFolder, '\', days(d,1).name];
        cd(dayFolder)
        if exist('newTrialdata.mat','file')==2
            load('newTrialdata.mat')
            try
                xtraindata.eyelidpos = [xtraindata.eyelidpos;trials.eyelidpos];
            catch ME
                addCols = size(xtraindata.eyelidpos,2) - size(trials.eyelidpos,2);
                if addCols <= 1 % trials is bigger than xtrain
                    addCols = addCols*-1;
                    addMe = nan(size(xtraindata.eyelidpos,1),addCols);
                    xtraindata.eyelidpos = [xtraindata.eyelidpos, addMe];
                elseif addCols >= 1 % xtraindata is bigger than trials                    
                    addMe = nan(size(trials.eyelidpos,1),addCols);
                    trials.eyelidpos = [trials.eyelidpos, addMe];
                end
                xtraindata.eyelidpos = [xtraindata.eyelidpos;trials.eyelidpos];
            end
            xtraindata.csdur = [xtraindata.csdur; trials.c_csdur];
            xtraindata.lasdur = [xtraindata.lasdur; trials.laser.dur];
            xtraindata.usdur = [xtraindata.usdur; trials.c_usdur];
            xtraindata.lasdel = [xtraindata.lasdel; trials.laser.delay];
            xtraindata.lasamp = [xtraindata.lasamp; trials.laser.amp];
            
            addOnes = ones(length(trials.c_csdur),1);
            mouse = str2double(mice{m,1}(3:end));
            day = str2double(days(d,1).name);
            xtraindata.mouse = [xtraindata.mouse; addOnes*mouse];
            xtraindata.day = [xtraindata.day; addOnes*day];
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
save('acrossTrainingData.mat', 'xtraindata')