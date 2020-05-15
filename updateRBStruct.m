function [phase]=updateRBStruct(trainingSession, thisMouse, data, phase, timeVector)
    % get all the data from this mouse
    midx = strcmpi(data.mouse, thisMouse);
    eyelidpos = data.eyelidpos(midx,:);
    sessdate = data.date(midx,:);
    trialtype = data.type(midx,:);
    csdur = data.csdur(midx,:);
    laserdur = data.laserdur(midx,:);
    laserint = data.laserint(midx,:);
    
    didx = sessdate==trainingSession & trialtype==1 & csdur==0 & laserdur>0;
    
    temp = find(didx);
    laserOff = laserdur(temp(1),1)./1000;
    laserOffBin = find(timeVector>laserOff-0.00488372 & timeVector<laserOff+0.00488372);
    laserOffBin = max(laserOffBin);

    
    laserints = unique(laserint(didx,1));
    if laserints(1)<laserints(2) && laserints(2)<laserints(3)
        %disp('laser intensities in order')
    else
        disp('ASSUMPTION ABOUT LASER INTENSITY ORDER FALSIFIED')
        pause
    end
    for i = 1:length(laserints)
        idx = find(laserint==laserints(i,1) & didx==1);
        
        baseline = nan(length(idx),1);
        rbadjamp = nan(length(idx),1);
        rbtrace = nan(length(idx),440);
        for t = 1:length(idx)
            baseline(t,1) = mean(eyelidpos(idx(t,1),laserOffBin-10:laserOffBin));
            rbadjamp(t,1) = max(eyelidpos(idx(t,1),laserOffBin+1:end))-baseline(t,1);
            rbtrace(t,:) = eyelidpos(idx(t,1),:)-baseline(t,1);
        end
        
        phase.mouse = [phase.mouse; thisMouse];
        switch i
            case 1
                phase.laserint = [phase.laserint; 15];
            case 2
                phase.laserint = [phase.laserint; 30];
            case 3
                phase.laserint = [phase.laserint; 60];
            case 4
                phase.laserint = [phase.laserint; NaN];
        end
        phase.rbamp = [phase.rbamp;nanmean(rbadjamp)];
        phase.rbprob = [phase.rbprob;sum(rbadjamp>0.1)/length(idx)];
        phase.rbtrace = [phase.rbtrace;nanmean(rbtrace)];
        
    end
    
end