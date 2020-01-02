basedir = 'E:\pcp2ChR2 data\rebound';
mice = {'OK211';'OK213';'OK214';'OK215';'OK216';'OK217';'OK218'};

daily.mouse = {};
daily.crprob = [];
daily.cradjamp = [];
daily.cradjampHit = [];
daily.sesstype = []; % 0 = training, 1 = extinction, 2 = laser, 3 = laser extinction trial
for m = 1:length(mice)
    mousedir = [basedir,'\',mice{m,1}];
    cd(mousedir)
    
    days = dir('19*');
    for d = 1:length(days)
        
        daydir = [mousedir,'\',days(d,1).name];
        cd(daydir)
        
        % check if there is a trialdata file present
        loaded = false;
        if exist('newTrialdata.mat','file')==2
            load('newTrialdata.mat')
            loaded = true;
        elseif exist('trialdata.mat','file')==2
            load('trialdata.mat')
            loaded = true;
        end
        
        if loaded == true
            % check if this day had CS trials or if it was a laser test day
            numCSTrials = sum(trials.c_csdur>0);
            
            if numCSTrials > 0
                CSUSidx = find(trials.c_csdur>0 & trials.c_usdur>0);
                CSOnlyidx = find(trials.c_csdur>0 & trials.c_usdur==0);
                laserIdx = find(trials.laser.dur>0);               
                if ~isempty(CSUSidx)
                    checkidx = CSUSidx;
                    if length(laserIdx)>length(CSUSidx)                        
                        daily.sesstype(end+1,1) = 3;
                    else
                        daily.sesstype(end+1,1) = 0;
                    end
                else
                    checkidx = CSOnlyidx;
                    daily.sesstype(end+1,1) = 1;
                end
                
                baseline = nan(length(checkidx),1);
                cradjamp = nan(length(checkidx),1);
                stable = nan(length(checkidx),1);
                %eyeadj = nan(length(checkidx),size(trials.eyelidpos,2));
                for t = 1:length(checkidx)
                    baseline(t,1) = mean(trials.eyelidpos(checkidx(t,1),1:40));
                    stable(t,1) = max(trials.eyelidpos(checkidx(t,1),1:40))<0.3;
                    cradjamp(t,1) = max(trials.eyelidpos(checkidx(t,1),76:85))-baseline(t,1);
                    %eyeadj(t,:)=trials.eyelidpos(checkidx(t,1),:)-baseline(t,1);
                end
                
                meancradjamp = median(cradjamp(stable==1,1));
                meancradjampHit = median(cradjamp(stable==1 & cradjamp>=0.1));
                crprob = sum(cradjamp(stable==1,1)>=0.1)./sum(stable==1);
                
%                 figure
%                 plot(eyeadj(stable==1,:)')
%                 hold on
%                 plot(mean(eyeadj(stable==1,:)),'Color',[0 0 0],'LineWidth',3)
%                 title(num2str(meancradjamp))
%                 pause
                
                daily.mouse{end+1,1} = mice{m,1};
                daily.cradjamp(end+1,1) = meancradjamp;
                daily.cradjampHit(end+1,1) = meancradjampHit;
                daily.crprob(end+1,1) = crprob;
                
                clear trials meancradjamp crprob baseline cradjamp stable CSUSidx...
                    checkidx CSOnlyidx
            else
                daily.mouse{end+1,1} = mice{m,1};
                daily.sesstype(end+1,1) = 2;
                daily.cradjamp(end+1,1) = nan;
                daily.crprob(end+1,1) = nan;
                daily.cradjampHit(end+1,1) = nan;
            end
            
            clear numCSTrials
            
        end
        
        clear loaded daydir
    end
    
    clear mousedir
end


% for m = 1:length(mice)
%     thisMouse = mice(m,1);
%     idx = [];
%     for i = 1:length(daily.mouse)
%         if strcmpi(thisMouse, daily.mouse{i,1})
%             idx(end+1)=i;
%         end
%     end
%     viewData = [daily.sesstype(idx), daily.crprob(idx), daily.cradjamp(idx), daily.cradjampHit(idx)];
%     pause
%     clear idx
% end

% let performance threshold be >80% CRs and >50% median amp of hit trials
threshProb = 0.8;
threshAmp = 0.5;
numdays.mouse = mice;
numdays.acq = nan(7,1);
numdays.reacq = nan(7,1);
for m = 1:length(mice)
    thisMouse = mice(m,1);
    idx = [];
    for i = 1:length(daily.mouse)
        if strcmpi(thisMouse, daily.mouse{i,1})
            idx(end+1)=i;
        end
    end
    
    tempArr = [daily.sesstype(idx), daily.crprob(idx), daily.cradjamp(idx), ...
        daily.cradjampHit(idx)];
    
    % get days to acquisition
    counter = 1;
    iteridx = 2;
    while tempArr(iteridx,1)==0
        if tempArr(iteridx,2)>=threshProb && tempArr(iteridx,4)>=threshAmp
            break
        else
            counter = counter+1;
        end
        iteridx = iteridx + 1;
        if iteridx > length(tempArr) || ~tempArr(iteridx,1)==0
            counter = nan;
            break
        end
    end
    numdays.acq(m,1) = counter;
    

    % get days to reacquisition
    trainidx = find(tempArr(:,1)==2);
    iteridx = trainidx(3)+1;
    counter = 1;    
    while tempArr(iteridx,1)==0
        if tempArr(iteridx,2)>=threshProb && tempArr(iteridx,4)>=threshAmp
            break
        else
            counter = counter+1;
        end
        iteridx = iteridx + 1;
        if iteridx > length(tempArr) || ~tempArr(iteridx,1)==0
            counter = nan;
            break
        end
    end
    numdays.reacq(m,1) = counter;
    
end

savingsRatio = numdays.reacq./numdays.acq;