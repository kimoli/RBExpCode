function [signrankp, ttestp]=makePlots_RBxTrain(daystats, datafield, mW, mouseStart,...
    mouseEnd, group, makefig)


preidx = find(daystats.phase == 1 & daystats.laspow == mW & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd); % pretest
postidx = find(daystats.phase == 2 & daystats.laspow == mW & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd); % post training
signrankp.pre_post = signrank(datafield(preidx), datafield(postidx)); % ns with correction, but <0.05
[~,ttestp.pre_post]=ttest(datafield(preidx), datafield(postidx)); % signif
if makefig == 1
    figure
end
if strcmpi(group, 'paired')
    postextidx = find(daystats.phase == 3 & daystats.laspow == mW & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd); % post extinction
    signrankp.pre_postext = signrank(datafield(preidx), datafield(postextidx));
    signrankp.post_postext = signrank(datafield(postidx), datafield(postextidx));
    [~,ttestp.pre_postext]=ttest(datafield(preidx), datafield(postextidx));
    [~,ttestp.post_postext]=ttest(datafield(postidx), datafield(postextidx));
    for i = 1:length(preidx)
        plot(1:3, [datafield(preidx(i)), datafield(postidx(i)),...
            datafield(postextidx(i))],'LineStyle','--', 'Color', [0 0 0])
        hold on
        scatter(1:3, [datafield(preidx(i)), datafield(postidx(i)),...
            datafield(postextidx(i))],10,...
            'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 1 1])
    end
else
    for i = 1:length(preidx)
        plot(1:2, [datafield(preidx(i)), datafield(postidx(i))],'LineStyle','--', 'Color', [0 0 0])
        hold on
        scatter(1:2, [datafield(preidx(i)), datafield(postidx(i))],10,...
            'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 1 1])
    end
end
end