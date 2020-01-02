function [prevals, postvals, premad, postmad]=makePlots_RBxPower(daystats, datafield, powers, mice, makefig)

prevals = nan(length(mice),length(powers));
postvals = nan(length(mice),length(powers));
for p = 1:length(powers)
    for m = 1:length(mice)
        preidx = find(daystats.phase == 1 & daystats.laspow == powers(p) & daystats.mouse==mice(m)); % pretest
        prevals(m,p) = mean(datafield(preidx));
        postidx = find(daystats.phase == 2 & daystats.laspow == powers(p) & daystats.mouse==mice(m)); % post training
        postvals(m,p) = mean(datafield(postidx));
    end
end

if makefig == 1
    figure
end

premad = mad(prevals,1);
% MAD Version
errorbar(1:length(powers), median(prevals), premad,'-b')
hold on
scatter(1:length(powers), median(prevals),10,...
    'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1])
postmad = mad(postvals,1);
errorbar(1:length(powers), median(postvals), postmad,'-r')
hold on
scatter(1:length(powers), median(postvals),10,...
    'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0])


% IQR Version
% 
% tempquants = quantile(prevals,3);
% iqr1 = tempquants(2,:)-tempquants(1,:);
% iqr2 = tempquants(3,:)-tempquants(2,:);
% errorbar(1:length(powers), median(prevals), iqr1,iqr2,'-b')
% hold on
% scatter(1:length(powers), median(prevals),10,...
%     'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1])
% 
% tempquants = quantile(postvals,3);
% iqr1 = tempquants(2,:)-tempquants(1,:);
% iqr2 = tempquants(3,:)-tempquants(2,:);
% errorbar(1:length(powers), median(postvals), iqr1,iqr2,'-r')
% hold on
% scatter(1:length(powers), median(postvals),10,...
%     'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0])
end