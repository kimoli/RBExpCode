function [rbampphase] = getAmp(phase, rbampphase, a, bestInt_prob, thisMouse)

idx = find(phase.laserint==bestInt_prob(a,1) & strcmpi(phase.mouse, thisMouse));
if ~isempty(idx)
    rbampphase(a,:) = phase.rbamp(idx,:);
end

end