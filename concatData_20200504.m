function [data] = concatData_20200504(basedir, mice, dates)

data.eyelidpos = [];
data.mouse = [];
data.day = [];
data.csdur = [];
data.lasdur = [];
data.usdur = [];
data.lasdel = [];
data.lasamp = [];
for i = 1:length(dates)
    mouse = mice{i,1};
    day = dates(i,1);
    if ~isnan(day)
        temp = [basedir, '\', mouse, '\', num2str(day)];
        cd(temp)
        clear temp
        try
            load('newTrialdata.mat')
        catch ME
            disp('newTrialdata missing')
            pause
        end
        
        % make eyelid position arrays concatenate-able (trials might have a
        % different duration for different animals)
        prevcols = size(data.eyelidpos,2);
        curcols = size(trials.eyelidpos,2);
        % I think the if statements here are technically unnecessary but they
        % make the code a little more legible
        if prevcols > curcols
            while size(trials.eyelidpos,2) < prevcols
                trials.eyelidpos = [trials.eyelidpos, ...
                    nan(size(trials.eyelidpos,1),1)];
            end
        elseif curcols > prevcols
            while size(data.eyelidpos,2) < curcols
                data.eyelidpos = [data.eyelidpos, ...
                    nan(size(data.eyelidpos,1),1)];
            end
        end
        
        % concatenate current eyelidposition & other information arrays to the
        % end of the previously collected data
        data.eyelidpos = [data.eyelidpos; trials.eyelidpos];
        data.csdur = [data.csdur; trials.c_csdur];
        data.lasdur = [data.lasdur; trials.laser.dur];
        data.usdur = [data.usdur; trials.c_usdur];
        data.lasdel = [data.lasdel; trials.laser.delay];
        data.lasamp = [data.lasamp; trials.laser.amp];
        
        while size(data.day,1) < size(data.csdur,1)
            data.mouse(end+1,1) = str2double(mouse(3:end));
            data.day(end+1,1) = day;
        end
    end
end

end