function [signrankp, ttestp]=makeSubplotSummary(daystats, startMouse, endMouse, group)

if strcmpi(group,'paired')
    xmax = 3.25;
else
    xmax = 2.25;
end

figure
subplot(2,3,1)
[signrankp.int60.amp, ttestp.int60.amp] = makePlots_RBxTrain(daystats, daystats.rb.amp, 60, startMouse,...
    endMouse, group,0);
ylabel('RB amp (FEC)')
ylim([0 1])
xlim([0.75 xmax])
set(gca, 'XTickLabels', {'pre train','','post train','','post ext'})
title([group, ' Training, 60 mW'])

subplot(2,3,4)
[signrankp.int60.prob, ttestp.int60.prob] = makePlots_RBxTrain(daystats, daystats.rb.prob, 60, startMouse,...
    endMouse, group,0);
ylabel('RB Prob')
ylim([0 1])
xlim([0.75 xmax])
set(gca, 'XTickLabels', {'pre train','','post train','','post ext'})

subplot(2,3,2)
[signrankp.int30.amp, ttestp.int30.amp] = makePlots_RBxTrain(daystats, daystats.rb.amp, 30, startMouse,...
    endMouse, group,0);
ylabel('RB amp (FEC)')
ylim([0 1])
xlim([0.75 xmax])
set(gca, 'XTickLabels', {'pre train','','post train','','post ext'})
title([group, ' Training, 30 mW'])

subplot(2,3,5)
[signrankp.int30.prob, ttestp.int30.prob] = makePlots_RBxTrain(daystats, daystats.rb.prob, 30, startMouse,...
    endMouse, group,0);
ylabel('RB Prob')
ylim([0 1])
xlim([0.75 xmax])
set(gca, 'XTickLabels', {'pre train','','post train','','post ext'})

subplot(2,3,3)
[signrankp.int15.amp, ttestp.int15.amp] = makePlots_RBxTrain(daystats, daystats.rb.amp, 15, startMouse,...
    endMouse, group,0);
ylabel('RB amp (FEC)')
ylim([0 1])
xlim([0.75 xmax])
set(gca, 'XTickLabels', {'pre train','','post train','','post ext'})
title([group, ' Training, 15 mW'])

subplot(2,3,6)
[signrankp.int15.prob, ttestp.int15.prob] = makePlots_RBxTrain(daystats, daystats.rb.prob, 15, startMouse,...
    endMouse, group,0);
ylabel('RB Prob')
ylim([0 1])
xlim([0.75 xmax])
set(gca, 'XTickLabels', {'pre train','','post train','','post ext'})

end