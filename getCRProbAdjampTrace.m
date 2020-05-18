function [meancradjamp, crprob, eyelidtrace]=getCRProbAdjampTrace(trials, idx)
baseline = nan(length(idx),1);
cradjamp = nan(length(idx),1);
stable = nan(length(idx),1);
eyelidposadj = nan(length(idx),size(trials.eyelidpos,2));
crvel = nan(length(idx),1);
for t = 1:length(idx)
    baseline(t,1) = mean(trials.eyelidpos(idx(t,1),1:40));
    cradjamp(t,1) = max(trials.eyelidpos(idx(t,1),85))-baseline(t,1);
    vel = diff(trials.eyelidpos(idx(t,1),:))./0.005;
    crvel(t,1) = median(vel(63:76));
    eyelidposadj(t,:) = trials.eyelidpos(idx(t,1),:)-baseline(t,1);
    stable(t,1) = max(trials.eyelidpos(idx(t,1),1:40))<0.3;
end

crprob = sum(stable==1 & cradjamp>0.1 & crvel>=1)./sum(stable==1);
meancradjamp = nanmean(cradjamp(stable==1,1));
eyelidtrace = nanmean(eyelidposadj(stable==1,:));
end