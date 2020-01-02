function makeEyetraceSubplots_collapsed(daystats, daycrstats, mouseStart, mouseEnd, group,timeVector)

if strcmpi(group,'paired')
    numcols=3;
    panel = [1,2,3;4,5,6];
else
    numcols = 2;
    panel = [1,2;3,4];
end


%% plot the eyelid traces centered on event of interest
% plot for lowest intensity
pretestidx = daystats.phase == 1 & daystats.laspow==15 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
posttestidx = daystats.phase == 2  & daystats.laspow==15 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;

figure
subplot(2,numcols,panel(1,1))
plot(timeVector, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 1 1],...
    'LineWidth', 1)
title('pre training')

subplot(2,numcols,panel(1,2))
plot(timeVector, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 1 1],...
    'LineWidth', 1)
title('post training')

if strcmpi(group,'paired')
    postextidx = daystats.phase == 3  & daystats.laspow==15 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
    
    subplot(2,numcols,panel(1,3))
    plot(timeVector, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 1 1],...
        'LineWidth', 1)
    title('post extinction')
end

% plot for middle intensity
pretestidx = daystats.phase == 1 & daystats.laspow==30& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
posttestidx = daystats.phase == 2  & daystats.laspow==30& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;

subplot(2,numcols,panel(1,1))
hold on
plot(timeVector, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [1 0 0],...
    'LineWidth', 1)

subplot(2,numcols,panel(1,2))
hold on
plot(timeVector, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [1 0 0],...
    'LineWidth', 1)

if strcmpi(group,'paired')
    postextidx = daystats.phase == 3  & daystats.laspow==30& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
    
    subplot(2,numcols,panel(1,3))
    hold on
    plot(timeVector, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [1 0 0],...
        'LineWidth', 1)
end

% plot for highest intensity
pretestidx = daystats.phase == 1 & daystats.laspow==60& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
posttestidx = daystats.phase == 2  & daystats.laspow==60& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;

subplot(2,numcols,panel(1,1))
plot(timeVector, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 0 1],...
    'LineWidth', 1)

subplot(2,numcols,panel(1,2))
plot(timeVector, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 0 1],...
    'LineWidth', 1)

if strcmpi(group,'paired')
    postextidx = daystats.phase == 3  & daystats.laspow==60& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
    
    subplot(2,numcols,panel(1,3))
    plot(timeVector, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 0 1],...
        'LineWidth', 1)
end

% plot all intensity means for all trials
pretestidx = daystats.phase == 1 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
posttestidx = daystats.phase == 2  & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
subplot(2,numcols,panel(1,1))
plot(timeVector, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)

subplot(2,numcols,panel(1,2))
plot(timeVector, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 2)

if strcmpi(group,'paired')
    postextidx = daystats.phase == 3  & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
    
    subplot(2,numcols,panel(1,3))
    plot(timeVector, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 0 0],...
        'LineWidth', 2)
end

% plot all intensity means for hit trials only
pretestidx = daystats.phase == 1 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
posttestidx = daystats.phase == 2  & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
subplot(2,numcols,panel(1,1))
plot(timeVector, nanmean(daystats.meanRBTrHit(pretestidx,:)), 'Color', [1 0 1],...
    'LineWidth', 2)
ylim([0 1])
xlim([-0.2 1.4])
ylabel('Eyelid Position (FEC)')

subplot(2,numcols,panel(1,2))
plot(timeVector, nanmean(daystats.meanRBTrHit(posttestidx,:)), 'Color', [1 0 1],...
    'LineWidth', 2)
ylim([0 1])
xlim([-0.2 1.4])

if strcmpi(group,'paired')
    postextidx = daystats.phase == 3  & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
    
    subplot(2,numcols,panel(1,3))
    plot(timeVector, nanmean(daystats.meanRBTrHit(postextidx,:)), 'Color', [1 0 1],...
        'LineWidth', 2)
    ylim([0 1])
    xlim([-0.2 1.4])
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

subplot(2,numcols,panel(2,1))
plot(timeVector, nanmean(daycrstats.meantr(pretestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 1)
ylim([0 1])
xlim([-0.2 1.4])
ylabel('CS')

subplot(2,numcols,panel(2,2))
plot(timeVector, nanmean(daycrstats.meantr(posttestidx,:)), 'Color', [0 0 0],...
    'LineWidth', 1)
ylim([0 1])
xlim([-0.2 1.4])

if strcmpi(group,'paired')
    a = find(daycrstats.phase==2.5);
    temp = intersect(a,b);
    postextidx = intersect(temp,c);
    
    subplot(2,numcols,panel(2,3))
    plot(timeVector, nanmean(daycrstats.meantr(postextidx,:)), 'Color', [0 0 0],...
        'LineWidth', 1)
    ylim([0 1])
    xlim([-0.2 1.4])
end

end