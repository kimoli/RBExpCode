function [datevector] = getDateVector(startDate, stopDate)
% assumes that you only ever let one month roll over

startstring = num2str(startDate);
startmonth = startstring(3:4);
stopstring = num2str(stopDate);
stopmonth = stopstring(3:4);
if startmonth==stopmonth
    % don't need to roll over months to make date vector
    datevector = [startDate:stopDate]';
else
    monthnum = str2double(startstring);
    if monthnum==9 || monthnum==4 || monthnum==6 || monthnum==11
        % month has 30 days
        lastDayInMonth = [startstring(1:4),'30'];
    elseif monthnum==2
        % month days depends on if year is a leap year (if it was 2019 or
        % 2020 when I did the experiment, for simplification)
        if strcmpi(startstring(1:2),'19')
            lastDayInMonth = [startstring(1:4),'27'];
        elseif strcmpi(startstring(1:2),'20')
            lastDayInMonth = [startstring(1:4),'28'];
        else
            disp('script not designed to accommodate this year')
        end
    else
        % month has 31 days
        lastDayInMonth = [startstring(1:4),'31'];
    end
    firstDayInMonth = [stopstring(1:4),'01'];
    datevector = [startDate:str2double(lastDayInMonth)]';
    temp = [str2double(firstDayInMonth):stopDate]';
    datevector = [datevector;temp];
end
end