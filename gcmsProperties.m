%gcmsProperties.m

% Switch which data set we are looking at:
global projectInvestigate;
projectInvestigate = 'Malaria-HealthyControl';

global pathToRawData    % directory where raw GC-MS data is kept, 
global pathToRawDataSub % sub-directory for raw GC-MS data, 
global pathToPeakData   % directory where peak detection results goes to
global pathToPeakDataSub    % sub-directory for peak detection result
% peak detection parameters
global normaliseData    % should the data be first normalised before peak detection
                    % 0 - raw data only
                    % 1 - normalise by sample, 
                    %     this one is no longer implemented
                    % 2 - normalise by blank sample of the same day
                    % 3 - normalise by all data analysed at the same time,
                    %     this requires a file "normFactorAllData.txt" in
                    %     the directory of the raw data
global numbBlank      % number of blank samples
global sigThresholdNormFactor     % factor to redefine the signal threshold 
                    % using the overall noise level of all data.
                    % [] - no redefining, 
                    % rule of thumb: 100 - easy data (e.g. plant); 
                    % between 2 and 10 - complex data (e.g. human breath)


if (strcmp(projectInvestigate, 'Malaria-HealthyControl')) 
	%==========    
	% file info
	pathToRawData = '/data/Malaria/HealthyControl/csvFiles/';
    pathToRawDataSub = '';
    pathToPeakData = '/results/Malaria/HealthyControl/';
    pathToPeakDataSub = '1.rawData/sigThreshold_10/';
    %==========
    normaliseData = 0;
    numbBlank = 0;
    sigThresholdNormFactor = 10;
end
