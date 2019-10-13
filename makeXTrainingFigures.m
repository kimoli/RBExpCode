close all
clear all

% go to data directory
machine = 'OREK';
if strcmpi(machine, 'OREK')
    basedir = 'E:\pcp2ChR2 data\rebound';
elseif strcmpi(machine, 'COMPUPITAR')
    basedir = 'D:\pcp2ChR2 data\rebound';
else
    disp('Please specify computer so we know what directory to use')
    basedir = '';
end
cd(basedir)
load('acrossTrainingData.mat')


% set up time bins & windows
timeVector = 1:size(xtraindata.eyelidpos,2);
timeVector = timeVector * 0.00488372;
timeVector = timeVector - 0.2;

crblwinstart = 1;
temp = timeVector - 0.21; % find first bin timestamped at or after 210 ms (time of puff onset)
temp2 = find(temp>=0);
crwinend = temp2(1);
temp = timeVector - 0.14; % find first bin timestamped at or after 140 ms (end of beta startle window [?])
temp2 = find(temp>=0);
crwinstart = temp2(1);
crblwinend = temp2(1)-1;
temp = timeVector - 0.850; % find first bin timestamped at or after 850 ms (time of laser offset)
temp2 = find(temp>=0);
rbwinstart = temp2(1);
rbwinend = length(timeVector);
rbblwinend = temp2(1)-1;
temp = timeVector - 0.8; % find first bin timestamped at or after 800 ms (50 ms before laser offset)
temp2 = find(temp>=0);
rbblwinstart = temp2(1);
clear temp temp2

% find animals
mice = unique(xtraindata.mouse);

% cycle through animals and get all relevant data
for m = 1:length(mice)
    mouseidx = find(xtraindata.mouse == mice(m,1));
    cr.day = [];
    cr.adjamp = [];
    cr.usdur =[];
    rb.day = [];
    rb.adjamp = [];
%     figure
    for i = 1:length(mouseidx)
        if max(xtraindata.eyelidpos(mouseidx(i),crblwinstart:crblwinend)) < 0.3 % arbitrary threshold for determining if trial baseline is stabls
            if xtraindata.csdur(mouseidx(i),1) > 0 % if trial has cs, consider CR
                crbl = mean(xtraindata.eyelidpos(mouseidx(i),crblwinstart:crblwinend));
                cr.adjamp(end+1,1) = max(xtraindata.eyelidpos(mouseidx(i),crwinstart:crwinend)) - crbl;
                cr.day(end+1,1) = xtraindata.day(mouseidx(i),1);
                cr.usdur(end+1,1) = xtraindata.usdur(mouseidx(i),1);
%                 plot(timeVector, xtraindata.eyelidpos(mouseidx(i),:))
%                 hold on
%                 plot([timeVector(1), timeVector(end)], [crbl crbl])
%                 plot([timeVector(1), timeVector(end)], [cr.adjamp(end)+crbl, cr.adjamp(end)+crbl])
%                 plot([0 0], [0 1])
%                 plot([0.21 0.21], [0 1])
%                 pause
%                 hold off
            elseif xtraindata.lasdur(mouseidx(i),1) > 0 % if trial has laser, consider RB
                rbbl = mean(xtraindata.eyelidpos(mouseidx(i),rbblwinstart:rbblwinend));
                rb.adjamp(end+1,1) = max(xtraindata.eyelidpos(mouseidx(i),rbwinstart:rbwinend)) - rbbl;
                rb.day(end+1,1) = xtraindata.day(mouseidx(i),1);
%                 plot(timeVector, xtraindata.eyelidpos(mouseidx(i),:))
%                 hold on
%                 plot([timeVector(1), timeVector(end)], [rbbl rbbl])
%                 plot([timeVector(1), timeVector(end)], [rb.adjamp(end)+rbbl, rb.adjamp(end)+rbbl])
%                 plot([0.85 0.85], [0 1])
%                 pause
%                 hold off
            end
        end
    end
    figure
    subplot(2,1,1)
    scatter(1:length(cr.adjamp), cr.adjamp)
    title(num2str(mice(m,1)))
    ylabel('cr amp')
    ylim([0 1])
    subplot(2,1,2)
    scatter(1:length(rb.adjamp), rb.adjamp)
    ylabel('rb amp')
    
    % get daily values
    days = [unique(cr.day); unique(rb.day)];
    days = unique(days);
    days = sort(days);
    daily.date = nan(length(days),1);
    daily.cradjamp = nan(length(days),1);
    daily.crprob = nan(length(days),1);
    daily.rbadjamp = nan(length(days),1);
    daily.rbprob = nan(length(days),1);
    daily.sesstype = nan(length(days),1); % 0 is laser test, 1 is training, 2 is extinction, 3 is savings
    extd = 0;
    for d = 1:length(days)
        daily.date(d,1) = days(d,1);
        
        cridx = cr.day == days(d,1);
        if sum(cridx)>0 % if there were cs trials that day
            daily.cradjamp(d,1) = mean(cr.adjamp(cridx,1));
            daily.crprob(d,1) = sum(cr.adjamp(cridx,1)>=0.1)./sum(cridx);
        end
        
        rbidx = rb.day == days(d,1);
        if sum(rbidx)>0
            daily.rbadjamp(d,1) = mean(rb.adjamp(rbidx,1));
            daily.rbprob(d,1) = sum(rb.adjamp(rbidx,1)>=0.1)./sum(rbidx);
        end
        
        if sum(rbidx)>0 && sum(cridx)==0 % laser
            daily.sesstype(d,1) = 0;
        elseif sum(rbidx)>0 && sum(cridx)>0 % training, extinction, or savings
            temp = find(cridx);
            if cr.usdur(temp(end),1)==0
                extd = 1;
                daily.sesstype(d,1) = 2;
            else
                if extd == 0
                    daily.sesstype(d,1) = 1;
                elseif extd == 1
                    daily.sesstype(d,1) = 3;
                end
            end
        end  
    end
    figure
    subplot(2,1,1)
    scatter(1:length(days), daily.cradjamp)
    hold on
    scatter(1:length(days), daily.rbadjamp)
    legend('cr', 'rb')
    ylabel('Amplitude (FEC)')
    title(mice(m,1))
    ylim([0 1])
    subplot(2,1,2)
    title(mice(m,1))
    scatter(1:length(days), daily.crprob)
    hold on
    scatter(1:length(days), daily.rbprob)
    legend('cr', 'rb')
    ylabel('Probability')
    ylim([0 1])
    
    
    trainidx = daily.sesstype==1;
    extidx = daily.sesstype==2;
    figure
    scatter(daily.cradjamp(trainidx), daily.rbadjamp(trainidx))
    hold on
    scatter(daily.cradjamp(extidx), daily.rbadjamp(extidx))
    title(mice(m,1))
    xlabel('cradjamp')
    ylabel('rbadjamp')
    xlim([0 1])
    ylim([0 1])
    legend('training', 'extinction')
    text(0.8, 0.1, ['r2=',num2str(corr(daily.cradjamp(trainidx),daily.rbadjamp(trainidx)))])
    text(0.8, 0.2, ['r2=',num2str(corr(daily.cradjamp(trainidx|extidx),daily.rbadjamp(trainidx|extidx)))])

    
    figure
    scatter(daily.crprob(trainidx), daily.rbprob(trainidx))
    hold on
    scatter(daily.crprob(extidx), daily.rbprob(extidx))
    title(mice(m,1))
    xlabel('crprob')
    ylabel('rbprob')
    xlim([0 1])
    ylim([0 1])
    legend('training', 'extinction')
    text(0.8, 0.1, ['r2=',num2str(corr(daily.crprob(trainidx),daily.rbprob(trainidx)))])
    text(0.8, 0.2, ['r2=',num2str(corr(daily.crprob(trainidx|extidx),daily.rbprob(trainidx|extidx)))])
    
    
    if m == 2
        figure
        for t = 4634:-1:1400
            plot(timeVector, xtraindata.eyelidpos(mouseidx(t),:))
            title([num2str(xtraindata.csdur(mouseidx(t),1)), ' ', num2str(xtraindata.lasdur(mouseidx(t),1))])
            hold on
            text(0.7, 0.9, num2str(t))
            ylim([0 1])
            pause
            hold off
        end
        
        % early laser no rebound
        figure
        plot(timeVector, xtraindata.eyelidpos(mouseidx(83),:))
        hold on
        plot([0 0], [0 1], 'LineStyle', '--', 'Color', [0 0 1])
        plot([0.85 0.85], [0 1], 'LineStyle', '--', 'Color', [0 0 1])
        ylim([0 1])
        
        % early cs no cr
        figure
        plot(timeVector, xtraindata.eyelidpos(mouseidx(73),:))
        hold on
        plot([0 0], [0 1], 'LineStyle', '--', 'Color', [1 0 0])
        plot([0.21 0.21], [0 1], 'LineStyle', '--', 'Color', [1 0 0])
        ylim([0 1])
        
        % middle CR with CR
        figure
        plot(timeVector, xtraindata.eyelidpos(mouseidx(1035),:))
        hold on
        plot([0 0], [0 1], 'LineStyle', '--', 'Color', [1 0 0])
        plot([0.21 0.21], [0 1], 'LineStyle', '--', 'Color', [1 0 0])
        ylim([0 1])
        
        % middle laser with RB
        figure
        plot(timeVector, xtraindata.eyelidpos(mouseidx(1037),:))
        hold on
        plot([0 0], [0 1], 'LineStyle', '--', 'Color', [0 0 1])
        plot([0.85 0.85], [0 1], 'LineStyle', '--', 'Color', [0 0 1])
        ylim([0 1])
        
        % late CS with CR
        figure
        plot(timeVector, xtraindata.eyelidpos(mouseidx(1200),:))
        hold on
        plot([0 0], [0 1], 'LineStyle', '--', 'Color', [1 0 0])
        plot([0.21 0.21], [0 1], 'LineStyle', '--', 'Color', [1 0 0])
        ylim([0 1])
        
        % late laser with RB
        figure
        plot(timeVector, xtraindata.eyelidpos(mouseidx(1381),:))
        hold on
        plot([0 0], [0 1], 'LineStyle', '--', 'Color', [0 0 1])
        plot([0.85 0.85], [0 1], 'LineStyle', '--', 'Color', [0 0 1])
        ylim([0 1])
        
        % ext CS without CR
        figure
        plot(timeVector, xtraindata.eyelidpos(mouseidx(4443),:))
        hold on
        plot([0 0], [0 1], 'LineStyle', '--', 'Color', [1 0 0])
        plot([0.21 0.21], [0 1], 'LineStyle', '--', 'Color', [1 0 0])
        ylim([0 1])
        
        % ext laser with RB
        figure
        plot(timeVector, xtraindata.eyelidpos(mouseidx(4493),:))
        hold on
        plot([0 0], [0 1], 'LineStyle', '--', 'Color', [0 0 1])
        plot([0.85 0.85], [0 1], 'LineStyle', '--', 'Color', [0 0 1])
        ylim([0 1])
        
        
    end
end