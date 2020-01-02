function makeEyetraceSubplots_collapsed(daystats, daycrstats, mouseStart, mouseEnd, group,timeVector)

if strcmpi(group,'paired')
    numcols=3;
    panel = [1,2,3;4,5,6];
else
    numcols = 2;
    panel = [1,2;3,4];
end


%% plot the eyelid traces centered on event of interest
% % plot for lowest intensity
% pretestidx = daystats.phase == 1 & daystats.laspow==15 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
% posttestidx = daystats.phase == 2  & daystats.laspow==15 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
% 
% figure
% subplot(2,numcols,panel(1,1))
% plot(timeVector, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 1 1],...
%     'LineWidth', 1)
% title('pre training')
% 
% subplot(2,numcols,panel(1,2))
% plot(timeVector, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 1 1],...
%     'LineWidth', 1)
% title('post training')
% 
% if strcmpi(group,'paired')
%     postextidx = daystats.phase == 3  & daystats.laspow==15 & daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
%     
%     subplot(2,numcols,panel(1,3))
%     plot(timeVector, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 1 1],...
%         'LineWidth', 1)
%     title('post extinction')
% end
% 
% % plot for middle intensity
% pretestidx = daystats.phase == 1 & daystats.laspow==30& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
% posttestidx = daystats.phase == 2  & daystats.laspow==30& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
% 
% subplot(2,numcols,panel(1,1))
% hold on
% plot(timeVector, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [1 0 0],...
%     'LineWidth', 1)
% 
% subplot(2,numcols,panel(1,2))
% hold on
% plot(timeVector, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [1 0 0],...
%     'LineWidth', 1)
% 
% if strcmpi(group,'paired')
%     postextidx = daystats.phase == 3  & daystats.laspow==30& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
%     
%     subplot(2,numcols,panel(1,3))
%     hold on
%     plot(timeVector, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [1 0 0],...
%         'LineWidth', 1)
% end
% 
% % plot for highest intensity
% pretestidx = daystats.phase == 1 & daystats.laspow==60& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
% posttestidx = daystats.phase == 2  & daystats.laspow==60& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
% 
% subplot(2,numcols,panel(1,1))
% plot(timeVector, nanmean(daystats.meanRBTr(pretestidx,:)), 'Color', [0 0 1],...
%     'LineWidth', 1)
% 
% subplot(2,numcols,panel(1,2))
% plot(timeVector, nanmean(daystats.meanRBTr(posttestidx,:)), 'Color', [0 0 1],...
%     'LineWidth', 1)
% 
% if strcmpi(group,'paired')
%     postextidx = daystats.phase == 3  & daystats.laspow==60& daystats.mouse>=mouseStart & daystats.mouse<=mouseEnd;
%     
%     subplot(2,numcols,panel(1,3))
%     plot(timeVector, nanmean(daystats.meanRBTr(postextidx,:)), 'Color', [0 0 1],...
%         'LineWidth', 1)
% end

% plot all intensity means for all trials
thisMouse = mouseStart;
pretestTraces = [];
posttestTraces = [];
while thisMouse <= mouseEnd
    pretestidx = daystats.phase == 1 & daystats.mouse==thisMouse;
    if sum(pretestidx)>0
        posttestidx = daystats.phase == 2 & daystats.mouse==thisMouse;
        pretestTraces = [pretestTraces; nanmean(daystats.meanRBTr(pretestidx,:))];
        posttestTraces = [posttestTraces; nanmean(daystats.meanRBTr(posttestidx,:))];
    end
    thisMouse = thisMouse + 1;
end
subplot(2,numcols,panel(1,1))
shadedErrorBar(timeVector, nanmean(pretestTraces), ...
    nanstd(pretestTraces)./sqrt(size(pretestTraces,1)),'-k')

subplot(2,numcols,panel(1,2))
shadedErrorBar(timeVector, nanmean(posttestTraces), ...
    nanstd(posttestTraces)./sqrt(size(posttestTraces,1)),'-k')

if strcmpi(group,'paired')
    postextTraces = [];
    thisMouse = mouseStart;
    while thisMouse <= mouseEnd
        postextidx = daystats.phase == 3 & daystats.mouse==thisMouse;
        if sum(postextidx)>0
            postextTraces = [postextTraces; nanmean(daystats.meanRBTr(postextidx,:))];
        end
        thisMouse = thisMouse + 1;
    end
    subplot(2,numcols,panel(1,3))
    shadedErrorBar(timeVector, nanmean(postextTraces), ...
        nanstd(postextTraces)./sqrt(size(postextTraces,1)),'-k')
end

% plot all intensity means for hit trials only
thisMouse = mouseStart;
pretestTraces = [];
posttestTraces = [];
while thisMouse <= mouseEnd
    pretestidx = daystats.phase == 1 & daystats.mouse==thisMouse;
    if sum(pretestidx)>0
        posttestidx = daystats.phase == 2 & daystats.mouse==thisMouse;
        pretestTraces = [pretestTraces; nanmean(daystats.meanRBTrHit(pretestidx,:))];
        posttestTraces = [posttestTraces; nanmean(daystats.meanRBTrHit(posttestidx,:))];
    end
    thisMouse = thisMouse + 1;
end
subplot(2,numcols,panel(1,1))
hold on
shadedErrorBar(timeVector, nanmean(pretestTraces), ...
    nanstd(pretestTraces)./sqrt(size(pretestTraces,1)),'-r')
ylim([0 1])
xlim([-0.2 1.4])
ylabel('Eyelid Position (FEC')

subplot(2,numcols,panel(1,2))
hold on
shadedErrorBar(timeVector, nanmean(posttestTraces), ...
    nanstd(posttestTraces)./sqrt(size(posttestTraces,1)),'-r')
ylim([0 1])
xlim([-0.2 1.4])

if strcmpi(group,'paired')
    postextTraces = [];
    thisMouse = mouseStart;
    while thisMouse <= mouseEnd
        postextidx = daystats.phase == 3 & daystats.mouse==thisMouse;
        if sum(postextidx)>0
            postextTraces = [postextTraces; nanmean(daystats.meanRBTrHit(postextidx,:))];
        end
        thisMouse = thisMouse + 1;
    end
    subplot(2,numcols,panel(1,3))
    hold on
    shadedErrorBar(timeVector, nanmean(postextTraces), ...
        nanstd(postextTraces)./sqrt(size(postextTraces,1)),'-r')
end
ylim([0 1])
xlim([-0.2 1.4])

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