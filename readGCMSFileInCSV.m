function [data, massZ] = readGCMSFileInCSV(filename)
%
% function to read the GC file that's written into csv format
%
% -- filename: string for the filename of the GC file you want to read
%
% The data output are: 
% -- data : nTimeIndex x nMassSpec
% -- massZ : nMassSpec x 1
%
% XRW, CSIRO
% updated 14 Feb 2014

% open file
fid = fopen(filename);
    
% this is the csv file of converting QTOF data into csv, using matlab
% code. thus the save csv is readable as straight read in matlab
%   The first row records the masses, except the first cell, which is
% NaN. The first column records the timestamps of each measurement.
% Subsequent cells records the measurement of each (time, mass) tuple.
%   The other files type save the data transposed from this data, so
% need to transpose this file.

dataAll = load(filename);
massZ = dataAll(1, 2:end);
data = dataAll(2:end, :)';

fclose(fid);