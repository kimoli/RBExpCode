function writeDataForRBatchPairwise(basestring, basefilename, header, headeridx,...
    compdata)
% headeridx: which column does the header name go into? this will also be
% the column that compdata(:,1) goes into
    for i = 2:size(compdata,2)
        if i<size(compdata,2)-1
            headerstring = [basestring, num2str(i-1)];
        elseif i==size(compdata,2)-1
            if contains(basestring, 'Reacq')
                headerstring = [basestring, num2str(i-1)];
            else
                headerstring = [basestring, 'nMinus1'];
            end
        else
            headerstring = [basestring, 'Last'];
        end
        if headeridx==1
            headers = {header, headerstring};
            values = [compdata(:,1), compdata(:,i)];
        else
            headers = {headerstring, header};
            values = [compdata(:,i), compdata(:,1)];
        end
        tempcsv = [headers;num2cell(values)];
        tempcsv = cell2table(tempcsv(2:end,:),'VariableNames',tempcsv(1,:));
        filename = [basefilename, num2str(i-1),'_forPairwise.csv'];
        writetable(tempcsv,filename)
    end   
end