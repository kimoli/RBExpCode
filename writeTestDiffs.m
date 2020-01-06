function writeTestDiffs(data)
temp = data(:,1)-data(:,2);
csvwrite('temp.csv',temp')
end