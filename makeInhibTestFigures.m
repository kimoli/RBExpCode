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
load('overallData.mat')

mice = [211,213,214,215,216,217,218];
timeVector = 1:size(testInhibData.eyelidpos,2);
timeVector = timeVector * 0.00488372;
timeVector = timeVector - 0.2;


for m = 1:length(mice)
    thisMouse = mice(m);
    mouseIdx = find(testInhibData.mouse==thisMouse);
    thisDay = testInhibData.day(mouseIdx);
    if length(unique(thisDay))>1
        disp('bad data?')
        pause
    end
    
    if unique(testInhibData.lasdur(mouseIdx)) > 2
        disp('mouse has more then 2 testes laser durs')
        pause
    end
    
    csLasTrials = testInhibData.mouse==thisMouse & testInhibData.lasdur>0 & ...
        testInhibData.csdur>0;
    csNoLasNoPuffTrials = testInhibData.mouse==thisMouse & testInhibData.lasdur==0 & ...
        testInhibData.csdur>0 & testInhibData.usdur==0;
    
    csLasData = testInhibData.eyelidpos(csLasTrials,:);
    lasdel = testInhibData.lasdel(csLasTrials,:)./1000;
    lasdur = testInhibData.lasdur(csLasTrials,:)./1000;
    csNoLasNoPuffData = testInhibData.eyelidpos(csNoLasNoPuffTrials,:);
    
    for i = 1:size(csLasData,1)
        bl = mean(csLasData(i, 1:40));
        csLasData(i,:) = csLasData(i,:) - bl;
    end
    for i = 1:size(csNoLasNoPuffData,1)
        bl = mean(csNoLasNoPuffData(i, 1:40));
        csNoLasNoPuffData(i,:) = csNoLasNoPuffData(i,:) - bl;
    end
    
    figure
    subplot(2,1,1)
    plot(timeVector, csLasData)
    hold on
    plot(timeVector, mean(csLasData), 'LineWidth', 3)
    plot([lasdel(1), lasdel(1)], [0 1])
    plot([lasdel(1)+lasdur(1), lasdel(1)+lasdur(1)], [0 1])
    title(num2str(thisMouse))
    xlabel('Time from Tone (ms)')
    xlim([-0.2 1.4])
    ylim([-0.1 1])
    subplot(2,1,2)
    plot(timeVector, csNoLasNoPuffData)
    hold on
    plot(timeVector, mean(csNoLasNoPuffData), 'LineWidth', 3)
    xlabel('Time from Tone (ms)')
    ylabel('Eyelid Position (FEC)')
    xlim([-0.2 1.4])
    ylim([-0.1 1])
    
end
