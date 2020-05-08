function [rbstats, crstats, daystats, daycrstats] = getrbprops(inputData, phasename, rbstats,...
    daystats, crstats, daycrstats, m, mice, timeVector, checkfigs)

checkidx = inputData.mouse == mice(m,1) & inputData.lasdur>0 &...
    inputData.csdur == 0 & inputData.usdur==0;

% are there multiple laser intensities for this day?
lasints = unique(inputData.lasamp(checkidx)); % no calibration trial will be included in the data sent to this function, so # unique vals will be the # of laser intensities in the session
% mice(m,1)
% length(lasints)
% sum(checkidx)
% pause
%% CR stuff
checkidx = inputData.mouse == mice(m,1) & inputData.lasdur==0 &...
        inputData.csdur > 0;
cstrials = find(checkidx);
clear checkidx
temp = find(timeVector > 0);
cswinstart = temp(1);
clear temp
temp = find((timeVector - 0.05)>0);
crwinstart = temp(1);
clear temp
temp = find((timeVector - 0.215)>0);
crwinend = temp(1);
clear temp

baselines = nan(length(cstrials),1);
stable = nan(length(cstrials),1);
cramps = nan(length(cstrials),1);
for c = 1:length(cstrials)
    baselines(c,1) = mean(inputData.eyelidpos(cstrials(c),1:cswinstart-1));
    if max(inputData.eyelidpos(cstrials(c),1:cswinstart-1)) < 0.3
        stable(c,1) = 1;
    else
        stable(c,1) = 0;
    end
    cramps(c,1) = max(inputData.eyelidpos(cstrials(c),crwinstart:crwinend))-baselines(c,1);
%     figure
%     plot(inputData.eyelidpos(cstrials(c),:))
%     hold on
%     plot([0 200], [baselines(c,1) baselines(c,1)])
%     plot([0 200], [cramps(c,1)+baselines(c,1) cramps(c,1)+baselines(c,1)])
%     plot([crwinstart crwinstart], [0 1])
%     plot([crwinend crwinend], [0 1])
%     pause
%     close all
end

alltraces = inputData.eyelidpos(cstrials,:);
stabletraces = alltraces(stable==1,:);
meantrace = nanmean(stabletraces);
if size(meantrace,2)<340
    addcols = nan(1,340-size(meantrace,2));
    meantrace = [meantrace,addcols];
end

%% update output variables with CR information

addmice = ones(length(cramps),1)*mice(m,1);
crstats.mouse = [crstats.mouse; addmice];
addphase = cell(length(cramps),1);
[addphase{1:end}]=deal(phasename);
crstats.phase = [crstats.phase;addphase];
crstats.amp = [crstats.amp;cramps];

daycrstats.mouse = [daycrstats.mouse; mice(m,1)];
daycrstats.phase = [daycrstats.phase; phasename];
daycrstats.amp = [daycrstats.amp; mean(cramps(stable==1,1))];
daycrstats.prob = [daycrstats.prob;sum(cramps(stable==1,1)>=0.1)./sum(stable)];
daycrstats.meantr = [daycrstats.meantr; meantrace];


clear addmice addphase cramps meantrace alltraces stabletrials

for I = 1:length(lasints)
    
    %% rebound data
    checkidx = inputData.mouse == mice(m,1) & inputData.lasdur>0 &...
        inputData.csdur == 0 & inputData.usdur==0 & inputData.lasamp==lasints(I);
    
    % get the indices for checking rb properties
    lasOffTime = (inputData.lasdur(checkidx) + inputData.lasdel(checkidx))./1000;
    rbwinstart = nan(length(lasOffTime),1);
    rbwinend = nan(length(lasOffTime),1); % rb appears to occur in the 600 ms after the laser turns off\
    blwinstart = nan(length(lasOffTime),1);
    % the window end is sometimes beyond the duration of the video, but I
    % think it should be ok as long as I remember to exclude the nans from
    % the analysis
    for t = 1:length(lasOffTime)
        tempvals = timeVector - lasOffTime(t);
        tempidx = find(tempvals>0);
        rbwinstart(t,1) = tempidx(1);
        clear tempvals tempidx
        
        tempvals = timeVector - (lasOffTime(t) + 0.6);
        tempidx = find(tempvals>0);
        rbwinend(t,1) = tempidx(1);
        clear tempvals tempidx
        
        tempvals = timeVector - (lasOffTime(t) - 0.05);
        tempidx = find(tempvals>0);
        blwinstart(t,1) = tempidx(1);
        clear tempvals tempidx
    end
    
    % collect information about the following rebound properties:
    %       amplitude, latency, probability
    theseData = inputData.eyelidpos(checkidx,:);
    bls = nan(length(lasOffTime),1);
    rbamps = nan(length(lasOffTime),1);
    rblats = nan(length(lasOffTime),1);
    for t = 1:length(lasOffTime)
        
        bls(t,1) = mean(theseData(t,blwinstart(t,1):rbwinstart(t,1)-1));
        % amplitude: max FEC within the window
        rbamps(t,1) = max(theseData(t,rbwinstart(t,1):rbwinend(t,1)))-bls(t,1);
        
        % latency: time for FEC to exceed the baseline (the 50 ms before
        % window onset) by 0.05
        temp = (theseData(t,:) - bls(t,1))-0.05;
        if sum(temp>0)>0
            temp2 = find(temp>0);
            latbin = temp2(1);
            if latbin < rbwinstart(t,1)
                latbin = nan;
            end
            while latbin < rbwinstart(t,1) && length(temp2)>1
                temp2 = temp2(2:end);
                latbin = temp2(1);
                if length(temp2) == 1
                    latbin = nan;
                    break
                end
            end
            if ~isnan(latbin)
                rblats(t,1) = timeVector(latbin)-lasOffTime(t,1);
            end
            clear temp temp2 latbin
        else
            clear temp
        end
        
        if checkfigs==1
            %check whether pulling out reasonable values
            figure
            plot(timeVector, theseData(t,:))
            hold on
            plot([timeVector(1) timeVector(end)], [rbamps(t,1)+bls(t,1), rbamps(t,1)+bls(t,1)], ...
                'LineStyle', '--')
            plot([timeVector(1) timeVector(end)], [bls(t,1)+0.05, bls(t,1)+0.05], ...
                'LineStyle', '--')
            plot([timeVector(1) timeVector(end)], [bls(t,1), bls(t,1)], ...
                'LineStyle', '--')
            if ~isnan(rblats(t,1))
                plot([rblats(t,1)+timeVector(rbwinstart(t,1))...
                    rblats(t,1)+timeVector(rbwinstart(t,1))], [0 1], ...
                    'LineStyle', '--')
            end
            pause
            close all
        end
    end
    
    % probability: # rbs >= 0.1 FEC over baseline / total trials with rb
    rbprob = sum(rbamps>=0.1)./length(rbamps);
    
    %% update the output variables

    addmice = ones(length(rbamps),1)*mice(m,1);
    addphase = cell(length(rbamps),1);
    [addphase{1:end}]=deal(phasename);
    addlasamps = ones(length(rbamps),1)*lasints(I);
    addlaspows = nan(length(rbamps),1);
    laspow = nan;
    %mice(m,1)
    if length(lasints)>=3
        switch I
            case 1
                addlaspows = ones(length(rbamps),1)*15;
                laspow = 15;
            case 2
                addlaspows = ones(length(rbamps),1)*30;
                laspow = 30;
            case 3
                addlaspows = ones(length(rbamps),1)*60;
                laspow = 60;
            case 4
                addlaspows = ones(length(rbamps),1)*70;
                laspow = 70;
        end
    end
    
    rbstats.mouse = [rbstats.mouse; addmice];
    rbstats.phase = [rbstats.phase; addphase];
    rbstats.amp = [rbstats.amp; rbamps];
    rbstats.lat = [rbstats.lat; rblats];
    rbstats.lasamp = [rbstats.lasamp; addlasamps];
    rbstats.laspow = [rbstats.laspow; addlaspows];
    rbstats.winstart = [rbstats.winstart; rbwinstart];
    
    clear addmice addphase addlasamps addlaspows
    
    daystats.mouse = [daystats.mouse; mice(m,1)];
    daystats.phase = [daystats.phase; phasename];
    daystats.rb.amp = [daystats.rb.amp; mean(rbamps)];
    daystats.rb.hitamp = [daystats.rb.hitamp; mean(rbamps(rbamps>0.1))];
    daystats.rb.lat = [daystats.rb.lat; nanmean(rblats)];
    daystats.rb.prob = [daystats.rb.prob; rbprob];
    daystats.lasamp = [daystats.lasamp; lasints(I)];
    daystats.laspow = [daystats.laspow; laspow];
    if size(theseData,1)==1
        %theseData= theseData-bls;
        daystats.meanRBTr = [daystats.meanRBTrHit; theseData];
        if rbamps>0.1
            daystats.meanRBTrHit = [daystats.meanRBTrHit; theseData];
        else
            daystats.meanRBTrHit = [daystats.meanRBTrHit; nan(1,size(theseData,2))];
        end
    else
%         for b = 1:length(bls)
%             theseData(b,:) = theseData(b,:)-bls(b,1);
%         end
        daystats.meanRBTr = [daystats.meanRBTr; mean(theseData)];
        if sum(rbamps>=0.1)==1
            daystats.meanRBTrHit = [daystats.meanRBTrHit; theseData(rbamps>=0.1,:)];
        elseif sum(rbamps>=0.1)>1
            daystats.meanRBTrHit = [daystats.meanRBTrHit; mean(theseData(rbamps>=0.1,:))];
        else
            daystats.meanRBTrHit = [daystats.meanRBTrHit; nan(1,size(theseData,2))];
        end
    end
    daystats.winstart = [daystats.winstart; rbwinstart(1)];
    
    clear rbamps rblats rbprob theseData bls laspow
end

end