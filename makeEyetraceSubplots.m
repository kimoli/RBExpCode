function makeEyetraceSubplots(daystats, daycrstats, mouseStart, mouseEnd, group,timeVector)

if strcmpi(group,'paired')
    numcols=3;
else
    numcols = 2;
end
panel = 1;


%% plot the eyelid traces centered on event of interest
% plot for lowest intensity
pretestidx = daystats.phase == 1 & daystats.laspow==15 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
posttestidx = daystats.phase == 2  & daystats.laspow==15 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;

figure
subplot(4,numcols,panel)
plot(timeVector-0.85, daystats.meanRBTr(pretestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector-0.85, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
ylabel('15 mW')
title('pre training')
xlim([-0.1 0.6])
panel = panel+1;

subplot(4,numcols,panel)
panel = panel+1;
plot(timeVector-0.85, daystats.meanRBTr(posttestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector-0.85, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
title('post training')
xlim([-0.1 0.6])

if strcmpi(group,'paired')
    postextidx = daystats.phase == 3  & daystats.laspow==15 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
    
    subplot(4,numcols,panel)
    panel = panel+1;
    plot(timeVector-0.85, daystats.meanRBTr(postextidx,:), 'Color', [0.5 0.5 0.5])
    hold on
    plot(timeVector-0.85, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 0 0],...
        'LineWidth', 2)
    ylim([0 1])
    title('post extinction')
    xlim([-0.1 0.6])
end

% plot for middle intensity
pretestidx = daystats.phase == 1 & daystats.laspow==30& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
posttestidx = daystats.phase == 2  & daystats.laspow==30& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;

subplot(4,numcols,panel)
panel = panel+1;
plot(timeVector-0.85, daystats.meanRBTr(pretestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector-0.85, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([-0.1 0.6])
ylabel('30 mW')

subplot(4,numcols,panel)
panel = panel+1;
plot(timeVector-0.85, daystats.meanRBTr(posttestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector-0.85, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([-0.1 0.6])

if strcmpi(group,'paired')
    postextidx = daystats.phase == 3  & daystats.laspow==30& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
    
    subplot(4,numcols,panel)
    panel = panel+1;
    plot(timeVector-0.85, daystats.meanRBTr(postextidx,:), 'Color', [0.5 0.5 0.5])
    hold on
    plot(timeVector-0.85, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 0 0],...
        'LineWidth', 2)
    ylim([0 1])
    xlim([-0.1 0.6])
end

% plot for highest intensity
pretestidx = daystats.phase == 1 & daystats.laspow==60& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
posttestidx = daystats.phase == 2  & daystats.laspow==60& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;

subplot(4,numcols,panel)
panel = panel+1;
plot(timeVector-0.85, daystats.meanRBTr(pretestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector-0.85, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([-0.1 0.6])
ylabel('60 mW')

subplot(4,numcols,panel)
panel = panel+1;
plot(timeVector-0.85, daystats.meanRBTr(posttestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector-0.85, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([-0.1 0.6])

if strcmpi(group,'paired')
    postextidx = daystats.phase == 3  & daystats.laspow==60& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
    
    subplot(4,numcols,panel)
    panel = panel+1;
    plot(timeVector-0.85, daystats.meanRBTr(postextidx,:), 'Color', [0.5 0.5 0.5])
    hold on
    plot(timeVector-0.85, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 0 0],...
        'LineWidth', 2)
    ylim([0 1])
    xlim([-0.1 0.6])
end

% plot crs
a = find(daycrstats.phase==1);
b = find(daycrstats.mouse>=mouseStart);
c = find(daycrstats.mouse<=mouseEnd);
temp = intersect(a,b);
pretestidx = intersect(temp,c);

a = find(daycrstats.phase==1.5);
temp = intersect(a,b);
posttestidx = intersect(temp,c);

a = find(daycrstats.phase==2.5);
temp = intersect(a,b);
postextidx = intersect(temp,c);

subplot(4,numcols,panel)
panel = panel+1;
plot(timeVector, daycrstats.meantr(pretestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daycrstats.meantr(pretestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([-0.1 0.6])
ylabel('CS')

subplot(4,numcols,panel)
panel = panel+1;
plot(timeVector, daycrstats.meantr(posttestidx,:), 'Color', [0.5 0.5 0.5])
hold on
plot(timeVector, nanmean(daycrstats.meantr(posttestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)
ylim([0 1])
xlim([-0.1 0.6])

if strcmpi(group,'paired')
    a = find(daycrstats.phase==2.5);
    temp = intersect(a,b);
    postextidx = intersect(temp,c);
    
    subplot(4,numcols,panel)
    panel = panel+1;
    plot(timeVector, daycrstats.meantr(postextidx,:), 'Color', [0.5 0.5 0.5])
    hold on
    plot(timeVector, nanmean(daycrstats.meantr(postextidx,:)), 'Color', [0 0 0],...
        'LineWidth', 2)
    ylim([0 1])
    xlim([-0.1 0.6])
end

end