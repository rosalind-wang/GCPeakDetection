% mainMS_peakDetection.m
% 
% Main script to process GC-MS data for peak detection. 
%
% Released under GNU General Public Licence v3. 
% 
% X. Rosalind Wang
% CSIRO
% October 2013

clear all

global pathToRawData
global pathToRawDataSub
global pathToPeakData  
global pathToPeakDataSub  
global normaliseData 
global numbBlank     
global sigThresholdNormFactor

% get all the global variables from the properties file
gcmsProperties

dirDataSub = [];
dirData = [pathToRawData, pathToRawDataSub, dirDataSub, '/'];
dirSave = [pathToPeakData, pathToPeakDataSub];

% get all the files in the raw data directory
files = dir(sprintf('%s/*.csv', dirData));
Nfiles = length(files);

% set up all the parameters for inputting into the peak detection algorithm
splitSignalSize = 100;
degPoly = 3;
thresholdSD = [0 10]; %5
thresholdFD = [0 5]; %5
sigThreshold = 2e4;
thresholdSignal = [1 sigThreshold]; %2
splitCoeultion = false;
showPlot = true;

% flag to save results 
% -- turn this off if we want to look at data and results as we go.
flagSave = false;
% other show plot flags
showPlotMZ = true;         % ind ms data and baseline
showPlotMain = true;        % result of each peak detection


% write all command window output to a text file
if flagSave
    diary([dirSave, 'allScreenOutputs.txt']);
end

% do we want to look at each mass spec channel's data and it's baseline?
if showPlotMZ
    hFigMZ = figure(4);
    props = {'LineWidth', 1};
    propsFont = {'FontSize', 12, 'FontWeight', 'bold'};
end

% do we want to look at the result of each peak detection?
if showPlotMain
    hFigMain = figure(5);   % simply because the function is using 1 and 2.
    props = {'LineWidth', 1};
    propsMarker = {'MarkerSize', 8};
    propsFont = {'FontSize', 12, 'FontWeight', 'bold'};
end

% results column
colPeakStart = 1;
colPeakEnd = 2;
colPeakMax = 3;
colPeakHeight = 4;
colPeakArea = 5;

% if we don't want to see the peak results, then output the time of start
% and end of this script
if ~showPlotMain && ~showPlot
    timeCurrent = clock;
    fprintf('Current time: %d-%d %d:%d:%.1f \n\n', ...
        timeCurrent(2), timeCurrent(3), ...
        timeCurrent(4), timeCurrent(5), timeCurrent(6));
end


% if we want to normalise the data by the blank sample
if normaliseData == 2
  dataNew = [];
  fprintf('\n\n');
  for nb = 1 : numbBlank
    fname = files(nb).name;
    [dataRaw, massZ] = readGCMSFileInCSV([dirData, fname]);
    fprintf('Loaded blank sample file: %s\n', fname);
    mzDataExpt = dataRaw(2:end, :);
    dataNew = [dataNew; mzDataExpt(mzDataExpt > 0)];
  end
  meanBlankSample = mean(dataNew);
  fprintf('    mean value of all non-zero measurement in blank sample is: %.4f\n', meanBlankSample);
end

if normaliseData == 3
    meanAllSamples = abs(load([dirData, 'normFactorAllData.txt']));    % TODO: check this, not sure about the abs()
    fprintf('Using previously calculated mean values for all data as normalisation factor. \n');
    fprintf('     normalisation factor is: %.4f\n', meanAllSamples);
end



%%

for fn = 1 : Nfiles

    % open each file 
    fname = files(fn).name;
    fprintf('\n\nPeak detection for data sample %s\n', fname);
    % put time in too: 
    timeCurrent = clock;
    fprintf('Current time: %d-%d %d:%d:%.1f \n\n', ...
        timeCurrent(2), timeCurrent(3), ...
        timeCurrent(4), timeCurrent(5), timeCurrent(6));
    
    % read the GCMS file
    [dataRaw, massZ] = readGCMSFileInCSV([dirData, fname]);
    
    Nmz = length(massZ);
    
    % initialise the results
    resMSPeak = cell(Nmz, 1);
    peakData = cell(Nmz, 1);
    
    % separate the raw data into time stamps and intensity
    mzDataExpt = dataRaw(2:end, :);
    timesteps = dataRaw(1,:);
    
    % if we want to normalise the data
%     if normaliseData
%         meanAllData = mean(mzDataExpt(:));
%         sigThresholdNorm = sigThreshold / meanAllData;
%         thresholdSignal = [1 sigThresholdNorm];   % TODO: error input?
%     end
    if normaliseData == 2
        % normalise the data
        mzDataExpt = mzDataExpt / meanBlankSample;
    end
    
    if normaliseData == 3
        mzDataExpt = mzDataExpt / meanAllSamples;    
    end
            
    % redefine the signal threshold using the overall noise level of all data
    if ~isempty(sigThresholdNormFactor)
        neighbour1 = mzDataExpt(:, 1:end-2);
        neighbour2 = mzDataExpt(:, 3:end);
        meanNeighbours = (neighbour1 + neighbour2) / 2;
        diffWithNeighbour = mzDataExpt(:, 2:end-1) - meanNeighbours;
        noiseValMZ = median(abs(diffWithNeighbour), 2);
        % TODO for the following, I think this signal threshold should be applied to all data
        sigThresholdNorm = mean(noiseValMZ) * sigThresholdNormFactor;     
        thresholdSignal = [1 sigThresholdNorm];
    end
    
    
    NPeaksOverAllMz = 0;
    
    for mzi = 1 : Nmz   
        signalMZ = mzDataExpt(mzi, :);
        if max(signalMZ) <= 10
            lambdaALS = 2e3;
        else
            lambdaALS = 5e3;
        end
        pALS = 1e-4;
        
%         % if we want to normalise the data
%         if normaliseData
%             signalMZ = signalMZ / meanAllData;
%         end
        
        % find baseline
        baselineMZ = baseCorrALS(signalMZ', lambdaALS, pALS);
        if showPlotMZ
            figure(hFigMZ);
            clf;
            hold on;
            plot(timesteps, signalMZ, 'b', props{:});
            plot(timesteps, baselineMZ, 'r', props{:});
            strTitle = sprintf('M/Z = %d', massZ(mzi));
            title(strTitle, propsFont{:});
            fprintf('press any key to continue......\n');
            pause;
        end
        
        % it seems it's a good idea to always take the baseline drift off
        % the reason for checking the drift was mainly to account for one
        % channel which had a large bleed for 1 minute. 
        signalMZbc = signalMZ - baselineMZ';
        
        % check to see if the signal is greater than noise or not
        if normaliseData == 2
            sigAboveNoiseLogical = (signalMZbc > sigThresholdNorm);
        else
            sigAboveNoiseLogical = (signalMZbc > sigThreshold);
        end
        totalSigAboveNoise = sum(sigAboveNoiseLogical);
        if totalSigAboveNoise > 5
            fprintf('------------------\n');
            fprintf('Now doing peak detection for M/Z = %d\n', massZ(mzi));
            resMSPeak{mzi} = GCPeakDetection(signalMZbc, timesteps, ...
                splitSignalSize, degPoly, thresholdSD, thresholdFD, ...
                thresholdSignal, splitCoeultion, showPlot, max(signalMZbc));
            % fprintf('\nFinished, press any key to continue...\n\n');
            % pause;
            
            %  total peaks found
            NPeaksTotal = size(resMSPeak{mzi}, 1);
            
            % if the total number of peaks isn't 0, then save this result
            if NPeaksTotal > 0
                % need to map the detected peaks into a vector of the same size as the
                % original data, that's all zeros except at peak maxima with values of
                % the peak area
                peakData{mzi} = zeros(length(timesteps), 1);
                for npeak = 1 : NPeaksTotal
                    % find the index in the timestep with the peak maximum time
                    indPeakMax = find(timesteps == resMSPeak{mzi}(npeak, colPeakMax));
                    peakData{mzi}(indPeakMax) = resMSPeak{mzi}(npeak, colPeakArea);
                end
                
                % save the peaks for this mz into the appropriate folder
                if flagSave
                    fnameSaveMZ = sprintf('%smz%d/%s', dirSave, massZ(mzi), fname);
                    csvwrite(fnameSaveMZ, peakData{mzi});
                end
            end
            
            % TODO: if showPlotMain -> plot all the resulting peaks. 
            if showPlotMain
                % plot all the peaks on the full signal
                figure(hFigMain);
                clf;
                hold on;
                plot(timesteps, signalMZ, 'b', props{:});
                grid on;
                propsMarker = {'MarkerSize', 8};
                plot(resMSPeak{mzi}(:, colPeakStart), zeros(NPeaksTotal, 1), ...
                    'm^', propsMarker{:}, 'MarkerFaceColor', 'm');
                plot(resMSPeak{mzi}(:, colPeakEnd), zeros(NPeaksTotal, 1), ...
                    'rv', propsMarker{:}, 'MarkerFaceColor', 'r');
                plot(resMSPeak{mzi}(:, colPeakMax), zeros(NPeaksTotal, 1), ...
                    'go', propsMarker{:}, 'MarkerFaceColor', 'g');
                strTitle = sprintf('M/Z = %d', massZ(mzi));
                title(strTitle, propsFont{:});
                fprintf('press any key to continue......\n');
                pause;
            end
            
            % keep a running total of peaks found for this sample
            NPeaksOverAllMz = NPeaksOverAllMz + NPeaksTotal;
            
        end
        
    end

    % save the full results
    if flagSave
        fnameSaveFull = [dirSave, '0fullResults/', fname(1:end-4), '.mat'];
        save(fnameSaveFull, 'resMSPeak', 'peakData', 'dataRaw', 'massZ');
    end
    
    % put some lines and space in the command window output
    fprintf('\n\n\n');
    fprintf('Finished peak detection of all mass/z for %s\n', fname);
    fprintf('Total number of peaks found over all mass was %d\n', NPeaksOverAllMz);
    fprintf('\n=========================================================');
    fprintf('\n\n\n');

end

% if we don't want to see the peak results, then output the time of start
% and end of this script
if ~showPlotMain && ~showPlot
    fprintf('\n\nAll files processed. \n\n');
    timeCurrent = clock;
    fprintf('Current time: %d-%d %d:%d:%.1f \n\n', ...
        timeCurrent(2), timeCurrent(3), ...
        timeCurrent(4), timeCurrent(5), timeCurrent(6));
end


% close diary
diary off

