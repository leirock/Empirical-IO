clear
clc
filename = 'E:\Dropbox\Documents\Academic\06_IO\04_CU_Emprical\hw2_BLP\Data\autoblp.csv';
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
autoblp = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;
save autoblp