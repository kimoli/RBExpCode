function [data] = concatData(basedir, xlsdata, col, idx)

data.eyelidpos = [];
data.mouse = {};
data.day = [];
data.csdur = [];
data.lasdur = [];
data.usdur = [];
data.lasdel = [];
data.lasamp = [];
for i = 1:length(idx)
    mouse = xlsdata{idx(i,1),1};
    day = cell2mat(xlsdata(idx(i,1),col));
    temp = [basedir, '\', mouse, '\', num2str(day)];
    cd(temp)
    clear temp
    load('newTrialdata.mat')
    
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
        data.mouse{end+1,1} = mouse;
        data.day(end+1,1) = day;
    end
    
end

end