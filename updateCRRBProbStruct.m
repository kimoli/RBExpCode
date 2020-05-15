function [phase]=updateCRRBProbStruct(trainingSessions, thisMouse, data, phase, timeVector, phasename)
    % get all the data from this mouse
    midx = strcmpi(data.mouse, thisMouse);
    eyelidpos = data.eyelidpos(midx,:);
    sessdate = data.date(midx,:);
    trialtype = data.type(midx,:);
    csdur = data.csdur(midx,:);
    laserdur = data.laserdur(midx,:);
    if strcmp(phasename, 'extinction')
        usdur = data.usdur(midx,:);
        startidx = 1;
        checkidx = find(sessdate==trainingSessions(startidx,1) & trialtype==1);
        while sum(csdur(checkidx,1)>0 & usdur(checkidx,1)>0)>0 % while there are still paired trials
            startidx = startidx + 1;
            checkidx = find(sessdate==trainingSessions(startidx,1) & trialtype==1);
        end
        if startidx > 1
            startidx = startidx - 1; % keep one baseline session
        end
    else
        startidx = 1;
    end
    sessioniter = 0;
    for d = startidx:length(trainingSessions)
        if sum(sessdate==trainingSessions(d,1))>0 % to detect if session is there
            didx = find(sessdate==trainingSessions(d,1) & trialtype==1 & csdur>0 & laserdur==0);
            sessioniter = sessioniter + 1;
            
            baseline = nan(length(didx),1);
            cradjamp = nan(length(didx),1);
            crvel = nan(length(didx),1);
            eyelidposadj = nan(length(didx),440);
            stable = nan(length(didx),1);
            for t = 1:length(didx)
                baseline(t,1) = mean(eyelidpos(didx(t,1),1:40));
                stable(t,1) = max(eyelidpos(didx(t,1),1:40))<0.3;
                cradjamp(t,1) = max(eyelidpos(didx(t,1),70:82))-baseline(t,1);
                eyelidposadj(t,:) = eyelidpos(didx(t,1),:)-baseline(t,1);
                vel = diff(eyelidpos(didx(t,1),:))./0.005;
                crvel(t,1) = median(vel(63:76));
            end
            phase.mouse = [phase.mouse; thisMouse];
            phase.session = [phase.session;sessioniter];
            phase.crprob = [phase.crprob;sum(stable==1 & cradjamp>=0.1 & crvel>1)./sum(stable==1)];
            phase.cradjamp = [phase.cradjamp;nanmean(cradjamp)];
            phase.eyelidposadj = [phase.eyelidposadj;nanmean(eyelidposadj)];
            
            clear baseline cradjamp eyelidposadj stable didx
            
            didx = find(sessdate==trainingSessions(d,1) & trialtype==1 & laserdur>0);
            % assumes that the laser occurs at the same time on all trials
            % where it is present, since that is how the experiment was
            % designed. I am also not tracking laser intensity because the
            % experiment was designed so that the laser presentations
            % throughout training were all 30 mW. Also assumes that the
            % laser is on 0 delay.
            if ~isempty(didx)
                laserOff = laserdur(didx(1),1)./1000;
                laserOffBin = find(timeVector>laserOff-0.00488372 & timeVector<laserOff+0.00488372);
                laserOffBin = max(laserOffBin);
                baseline = nan(length(didx),1);
                rbadjamp = nan(length(didx),1);
                rbtrace = nan(length(didx),440);
                for t = 1:length(didx)
                    baseline(t,1) = mean(eyelidpos(didx(t,1),laserOffBin-10:laserOffBin));
                    rbadjamp(t,1) = max(eyelidpos(didx(t,1),laserOffBin+1:end))-baseline(t,1);
                    rbtrace(t,:) = eyelidpos(didx(t,1),:)-baseline(t,1);
                end
                phase.rbamp = [phase.rbamp;nanmean(rbadjamp)];
                phase.rbprob = [phase.rbprob;sum(rbadjamp>0.1)/length(didx)];
                phase.rbtrace = [phase.rbtrace;nanmean(rbtrace)];
                clear laserOff laserOn laserOffBin baseline rbadjamp
            else
                phase.rbamp = [phase.rbamp; NaN];
                phase.rbprob = [phase.rbprob; NaN];
                phase.rbtrace = [phase.rbtrace; nan(1,440)];
            end
            clear didx
        end
    end
end