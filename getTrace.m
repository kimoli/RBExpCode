function [rbtracephase] = getTrace(phase, rbtracephase, a, bestInt_prob, thisMouse)

idx = find(phase.laserint==bestInt_prob(a,1) & strcmpi(phase.mouse, thisMouse));
if ~isempty(idx)
    rbtracephase(a,:) = phase.rbtrace(idx,:);
end

end