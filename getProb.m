function [rbprobphase] = getProb(phase, rbprobphase, a, bestInt_prob, thisMouse)

idx = find(phase.laserint==bestInt_prob(a,1) & strcmpi(phase.mouse, thisMouse));
if ~isempty(idx)
    rbprobphase(a,:) = phase.rbprob(idx,:);
end

end