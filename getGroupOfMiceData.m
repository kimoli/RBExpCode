function [pairedphasecrprob, pairedphasecradjamp, unpairedphasecrprob,...
    unpairedphasecradjamp, pairedphaseeyelidpos, unpairedphaseeyelidpos]=...
    getGroupOfMiceData(mice, paired, unpaired)

pairedphasecrprob = nan(20,2);
pairedphasecradjamp = nan(20,2);
unpairedphasecrprob = nan(20,2);
unpairedphasecradjamp = nan(20,2);
pairedphaseeyelidpos = nan(20,440);
unpairedphaseeyelidpos = nan(20,440);

for s = 1:20
    tempdata = nan(length(mice),1);
    for m = 1:length(mice)
        if ~isempty(paired.crprob(paired.session==s & ...
                strcmpi(paired.mouse,mice{m,1}),1))
            tempdata(m,1) = paired.crprob(paired.session==s & ...
                strcmpi(paired.mouse,mice{m,1}),1);
        end
    end
    pairedphasecrprob(s,1) = nanmedian(tempdata);
    pairedphasecrprob(s,2) = mad(tempdata,1);
    
    tempdata = nan(length(mice),1);
    for m = 1:length(mice)
        if ~isempty(paired.crprob(paired.session==s & ...
                strcmpi(paired.mouse,mice{m,1}),1))
            tempdata(m,1) = paired.cradjamp(paired.session==s & ...
                strcmpi(paired.mouse,mice{m,1}),1);
        end
    end
    pairedphasecradjamp(s,1) = nanmedian(tempdata);
    pairedphasecradjamp(s,2) = mad(tempdata,1);
    
    tempdata = nan(length(mice),440);
    for m = 1:length(mice)
        if ~isempty(paired.crprob(paired.session==s & ...
                strcmpi(paired.mouse,mice{m,1}),1))
            tempdata(m,:) = paired.eyelidposadj(paired.session==s & ...
                strcmpi(paired.mouse,mice{m,1}),:);
        end
    end
    pairedphaseeyelidpos(s,:) = nanmean(tempdata);
    
    
    tempdata = nan(length(mice),1);
    for m = 1:length(mice)
        if ~isempty(unpaired.crprob(unpaired.session==s & ...
                strcmpi(unpaired.mouse,mice{m,1}),1))
            tempdata(m,1) = unpaired.crprob(unpaired.session==s & ...
                strcmpi(unpaired.mouse,mice{m,1}),1);
        end
    end
    unpairedphasecrprob(s,1) = nanmedian(tempdata);
    unpairedphasecrprob(s,2) = mad(tempdata,1);
    
    tempdata = nan(length(mice),1);
    for m = 1:length(mice)
        if ~isempty(unpaired.crprob(unpaired.session==s & ...
                strcmpi(unpaired.mouse,mice{m,1}),1))
            tempdata(m,1) = unpaired.cradjamp(unpaired.session==s & ...
                strcmpi(unpaired.mouse,mice{m,1}),1);
        end
    end
    unpairedphasecradjamp(s,1) = nanmedian(tempdata);
    unpairedphasecradjamp(s,2) = mad(tempdata,1);
    
    tempdata = nan(length(mice),440);
    for m = 1:length(mice)
        if ~isempty(unpaired.crprob(unpaired.session==s & ...
                strcmpi(unpaired.mouse,mice{m,1}),1))
            tempdata(m,:) = unpaired.eyelidposadj(unpaired.session==s & ...
                strcmpi(unpaired.mouse,mice{m,1}),:);
        end
    end
    unpairedphaseeyelidpos(s,:) = nanmean(tempdata);
end
end