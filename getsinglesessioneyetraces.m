timeVector = 1:size(trials.eyelidpos,2);
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

eyelidtraces = [];
for i = 1:length(trials.c_csdur)
    if max(trials.eyelidpos(i,1:40))<=0.3 & trials.c_csdur>0 & trials.c_usdur>0
        eyelidtraces = [eyelidtraces; trials.eyelidpos];
    end
end

figure
plot(timeVector, nanmean(eyelidtraces), 'Color', [0.5 0.5 0.5])