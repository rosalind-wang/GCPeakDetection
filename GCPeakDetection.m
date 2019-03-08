function resPeaks = GCPeakDetection(gcData, timesteps, splitSignalSize, ...
    degPoly, thresholdSD, thresholdFD, thresholdSignal, splitCoeultion, ...
    showPlot, maxYLimPlot)
% 
% function to find 1D Chromatogram peaks using Savitzky-Golay smoothing.
%
% INPUTS:
% -- gcData : a vector containing the values of GC measurements
% -- timesteps : a vector of equal length to gcData, containing the elution
%      times of the recording. If this input is not given or is empty, then 
%      the default is time index. 
% -- splitSignalSize : 0 if not splitting the input signal into smaller
%      chunks for peak detection, default is 100 (input is empty).
% -- degPoly : the degree of polynomial for the SG smoothing, default is 3.
% -- thresholdSD : threshold for the second derivative, any SD below this
%      threshold is considered a candidate for peak. The input should be in
%      a vector of two values: [v1, v2], where v1 indicates
%      whether it's a absolute value threshold (1) or a relative value to
%      the noise (0), v2 is a value for the actual threshold. Default is [0
%      5], i.e. 5 times the noise. 
% -- thresholdFD : threshold for the first derivative, where we consider 
%      the start and end of any peak. The input is in the same format as
%      thresholdSD. Default is [0 2]. Note, [1 0.05] seems to be
%      reasonable for normalised data. 
% -- thresholdSignal : threshold above which the data is considered a
%      signal and not some impurity or noise. The input is in the same
%      format as thresholdSD, except if v1 = 0, then threshold is 
%      (baseline + v2 * noise). Default is [0 2].
% -- splitCoeultion : true if want to split the strong co-elution peaks,
%      use with caution. default if false
% -- showPlot : true if want to see the peak detection in action. default
%      is false. 
% -- maxYLimPlot : value for the maximum y-axis for the plot. default is
%      min(1000, max(gcData))
%
% OUTPUTS:
% -- resPeaks : a mx5 matrix of m peaks found in the input data, each row
%      is in the order of [peak start; peak end; peak maximum; peak height;
%      peak area], the first 3 values are either time index or elution time
%      depending on whether the second input is given or not. 
%
% Written according to the paper: Vivo-Truyols et al. "Automatic program
% for peak detection and deconvolution of multi-overlapped chromatographic
% signals, Part I: Peak detection", Journal of Chromatography A, vol. 1096
% (2005) pp 133-145
%
% Released under GNU General Public Licence v3. 
% 
% X. Rosalind Wang
% CSIRO
% September 2013

%--
% Initialisation, checking input data, and setting default values
%--
if size(gcData, 1) > 1 && size(gcData, 2) > 1
    error('data input should be a 1D vector.');
end
N = length(gcData);

if nargin < 2 || isempty(timesteps)
    timesteps = 1 : N;
end
if size(timesteps, 1) > 1 && size(timesteps, 2) > 1
    error('time step input should be a 1D vector.');
elseif length(timesteps) ~= N
    error('input data and time steps do not match in size.');
end

if nargin < 3 || isempty(splitSignalSize)
    splitSignalSize = 100;    
end
if splitSignalSize > N
    splitSignalSize = 0;    % use the whole signal
end

if nargin < 4 || isempty(degPoly)
    degPoly = 3;
end

if nargin < 5 || isempty(thresholdSD)
    % default is 5 * noise of SD
    thresholdSD = [0 5];
elseif length(thresholdSD) == 2
    % check to make sure the input threshold values are valid
    if ~(thresholdSD(1) ~= 0 || thresholdSD(1) ~= 1)
        error('Threshold for second derivative input 1 should be 0 or 1');
    end
else
    error('Threshold for second derivative not valid format');
end

if nargin < 6 || isempty(thresholdFD)
    % default is 2 * noise of FD
    thresholdFD = [0 2];  
elseif length(thresholdFD) == 2
    % check to make sure the input threshold values are valid
    if ~(thresholdFD(1) ~= 0 || thresholdFD(1) ~= 1)
        error('Threshold for first derivative first input should be 0 or 1');
    end
else
    error('Threshold for first derivative not valid format');
end

if nargin < 7 || isempty(thresholdSignal)
    % default is baseline + 2 * noise
    thresholdSignal = [0 2];    
elseif length(thresholdSignal) == 2
    % check to make sure the input threshold values are valid
    if ~(thresholdSignal(1) ~= 0 || thresholdSignal(1) ~= 1)
        error('Threshold for signal input 1 should be 0 or 1');
    end
else
    error('Threshold for signal not valid format');
end

if nargin < 8 || isempty(splitCoeultion)
    splitCoeultion = false;
end

if nargin < 9 
    showPlot = false;
end

if showPlot && nargin < 10
    if max(gcData) > 1000
        yULim = 1000;
    else
        yULim = max(gcData);
    end
elseif ~showPlot && nargin < 10
    yULim = max(gcData);
else 
    yULim = maxYLimPlot;
end

% initialise the output matrix
resPeaks = zeros(0, 5);     
colPeakStart = 1;
colPeakEnd = 2;
colPeakMax = 3;
colPeakHeight = 4;
colPeakArea = 5;

if showPlot
    props = {'LineWidth', 1};
    hFig = figure(1);
    clf;
    plot(timesteps, gcData, props{:});
    set(gca, 'YLim', [0 yULim]);
    grid on;
    hold on;
end

%==================
% STEP 1 : Dividing the input signal into managable chunks
%  -- run an initial SG smoothing to find signal's noise
%  -- divide the signal into roughly size of splitSignalSize, but making
%  sure that no actual signal is divided up. 

% initial smooth, using window size 5. 
smdInitial = sgolayfilt(gcData, degPoly, 5, 0);
% calculate the noise
% turns out using the noise criterion defined in the calcNoise function
% below doesn't work well here, so just going to use the median value
noiseSmdInitial = calcNoise(smdInitial);
baselineSmdInitial = median(smdInitial);
if ~thresholdSignal(1)
    thrSig = baselineSmdInitial + thresholdSignal(2) * noiseSmdInitial;
else
    thrSig = thresholdSignal(2);
end
% find all indices where the smoothed signal is less than noise x 2
% indNoiseSmdInitial = find(smdInitial < thrSig);
% plot the threshold for noise 
if showPlot
    figure(hFig);
    plot(timesteps, smdInitial, 'r:', props{:});
     plot([timesteps(1), timesteps(end)], [thrSig, thrSig], ...
        'b--', 'LineWidth', 1);
end

% find all the partitions, save the start of the partition indices
if splitSignalSize == 0
    % TODO: make sure this bit makes sense. 
    partStIndex = 1;
    partEndIndex = N;
else
    [partStIndex, partEndIndex] = partitionSignal(smdInitial, ...
        thrSig, N, splitSignalSize);
end
% total number of partitions
NPartition = length(partStIndex);
% plot these partitions out
if showPlot
    figure(hFig);
    for i = 1 : length(partStIndex)
        plot([timesteps(partStIndex(i)) timesteps(partStIndex(i))], [0 yULim], 'k--');
    end
end


% if the total number of partition found is too small, i.e. the signal
% threshold was set too low, then we need to reset that. 
% nb. this isn't in the referneced paper either. 
if splitSignalSize && NPartition <= 3
    % if the threshold was set to a certain value, then use noise as the
    % threshold, if we used noise as a threshold, then double the threshold
    % value
    if thresholdSignal(1) == 1
        thresholdSignal = [0 2];
    elseif thresholdSignal(1) == 0
        thresholdSignal(2) = thresholdSignal(2) * 2;
    else
        % we shouldn't get in here, but let's be safe
        error('Threshold for signal input 1 should be 0 or 1');
    end
    thrSig = baselineSmdInitial + thresholdSignal(2) * noiseSmdInitial;
    
    % recalculate the partitions
    [partStIndex, partEndIndex] = partitionSignal(smdInitial, ...
        thrSig, N, splitSignalSize);
    NPartition = length(partStIndex);
    
    % replot if the flag was on
    if showPlot
        figure(hFig);
        for i = 1 : length(partStIndex)
            plot([timesteps(partStIndex(i)) timesteps(partStIndex(i))], ...
                [0 yULim], 'k--');
        end
    end    
    
end

if showPlot
    fprintf('Chromatogram split into smaller sections. \n');
    fprintf('Press any key to continue... \n\n');
    pause;
end


%==================
% Find peaks for each partition of the data
NPeaksTotal = 0;
for partNumb = 1 : NPartition
    
    % get the data associated with this section
    dataPartOrig = gcData( partStIndex(partNumb) : partEndIndex(partNumb) );
    smdInitPart = smdInitial( partStIndex(partNumb) : partEndIndex(partNumb) );
    NDataPart = length(dataPartOrig);
    timePart = timesteps( partStIndex(partNumb) : partEndIndex(partNumb) );
    
    
    % If any signal in this partition is larger than the threshold for
    % signal we defined earlier, then we'll check for peaks
    if any(smdInitPart - thrSig > 0)
    
        %==================
        % STEP 2: find the best window size for the smoothing
        %  -- search through all odd size between 5 and 41
        %  -- use the Durbin-Watson (DW) test
        %  -- best ws is one which give DW closest to 2
        
        maxWinSize = 41;
        if partNumb == NPartition && length(dataPartOrig) < maxWinSize+4
            maxWinSize = length(dataPartOrig) - 5;
        end
        windowSize = 5 : 2 : maxWinSize;
        valDW = zeros(length(windowSize), 1);
        for wsi = 1 : length(windowSize)
            dataPartSmd = sgolayfilt(dataPartOrig, degPoly, windowSize(wsi), 0);
            diffSmdOrig = dataPartOrig - dataPartSmd;
            valDW(wsi) = ...
                sum( ( diffSmdOrig(2:end) - diffSmdOrig(1:end-1) ) .^2 ) ./ ...
                sum( diffSmdOrig .^2 ) * ( NDataPart / (NDataPart-1) );
        end
        % find the difference between the DW values and 2
        diffDWn2 = abs(valDW - 2);
        % find the actual window size
        [temp, indWS] = min(diffDWn2); %#ok<ASGLU>
        WS = windowSize(indWS);
%         fprintf('window size: %d\n', WS);
        
        % the smoothed data using the discovered window size
        % ZD - zero-th derivative
        smdZD = sgolayfilt(dataPartOrig, degPoly, WS, 0);
        
        if showPlot
            propsPart = {'LineWidth', 2};
            hFigPart = figure(2);
            clf;
            title(sprintf('Partition Number %d', partNumb), 'FontSize', 14);
            hold on
            plot(timePart, dataPartOrig, 'b-', propsPart{:});
            plot(timePart, smdZD, 'r:', propsPart{:});
            if max(dataPartOrig) > yULim
                set(gca, 'YLim', [0 yULim]);
            end
            grid on
        end
        
    
        %==================
        % STEP 3 : Peak Detection
        
        % We'll only do the peak detection if the signal in the partition is
        % more than the threshold set
    %---- comment out the previous "if any( > 0)" line and uncomment
    %the following if you want to see what the signal in the partition
    %looks like
    %if any(smdZD - thrSig > 0)
        
        %==================
        % STEP 3a. calculate derivatives
        
        % first derivative
        smdFD = sgolayfilt(dataPartOrig, degPoly, WS, 1);
        % second derivative
        smdSD = sgolayfilt(dataPartOrig, degPoly, WS, 2);
        % third derivative -- only to be used for moderate co-elutions
        smdTD = sgolayfilt(dataPartOrig, degPoly, WS, 3);
        % plot first and second derivative
        if showPlot
            figure(hFigPart);
            plot(timePart, smdFD, 'g--', propsPart{:});
            plot(timePart, smdSD, 'm-.', propsPart{:});
            legend('original data', 'SG smoothed', 'SG 1st derivative', ...
                'SG 2nd derivative');
        end
        
        %==================
        % STEP 3b. use noise of the SD to threshold the zones of -ve SD
        
        % calculate noise
        noiseSD = calcNoise(smdSD);
        % threshold is either some multiple of the noise, or a set value
        if ~thresholdSD(1)
            thrSD = thresholdSD(2) * noiseSD;
        else
            thrSD = thresholdSD(2);
        end        
        % plot the threshold for SD and the signal
        if showPlot
            figure(hFigPart);
            plot([timePart(1), timePart(end)], [-thrSD, -thrSD], ...
                'b--', 'LineWidth', 1);
            plot([timePart(1), timePart(end)], [thrSig, thrSig], ...
                '--', 'LineWidth', 1, 'Color', [0.33 0 0]);   % dark red
        end
        
        %==================
        % STEP 3c. find the negative zones of second derivative
        
        % find SD less than -thrSD
        indSD = find(smdSD < -thrSD);
        % check if there are any -ve zones
        if ~isempty(indSD)
            % find, within these indices, the zones of -ve SD, i.e. those with 
            % continues indices
            indIndSD_NegZoneEnd = [find(diff(indSD) > 1), length(indSD)];
            indIndSD_NegZoneSta = [1, indIndSD_NegZoneEnd(1:end-1)+1];
            % get the actual indices of these zones
            indSD_NegZoneEnd = indSD(indIndSD_NegZoneEnd);
            indSD_NegZoneSta = indSD(indIndSD_NegZoneSta);
            % number of negative zones found
            nSDNegZones = length(indSD_NegZoneEnd);
        else
            nSDNegZones = 0;
        end
        
        %==================
        % STEP 3d. use the negative zones of SD to find the peaks
        
        % don't think this is necessary. 
%         % recalculate the signal threshold for this partition
%         noiseSmd = calcNoise(smdZD);
%         baselineSmd = median(smdZD);
%         if ~thresholdSignal(1)
%             thrSig2 = baselineSmd + thresholdSignal(2) * noiseSmd;
%         else
%             thrSig2 = thresholdSignal(2);
%         end
        
        % need a threshold for the first derivative for start and end of
        % peak
        noiseFD = calcNoise(smdFD);
        if ~thresholdFD(1)
            thrFD = thresholdFD(2) * noiseFD;
        else
            thrFD = thresholdFD(2);
        end
        % check to make sure this threshold isn't too big 
        % don't think this is wise -- only applicable to cases where the
        % data is normalised. 
%         if thrFD > 0.05
%             thrFD = 0.05;
%         end
        % plot FD threshold
        if showPlot
            figure(hFigPart);
            plot([timePart(1), timePart(end)], [thrFD, thrFD], ...
                '--', 'LineWidth', 1, 'Color', [0 0.33 0]);   % dark green
            plot([timePart(1), timePart(end)], [-thrFD, -thrFD], ...
                '--', 'LineWidth', 1, 'Color', [0 0.33 0]);   % dark green
        end
        
        % need change of sign in FD, for checking end of the peak
        % -- this isn't in the ref'd paper.
        changeInSignFD = abs(diff(smdFD > 0));
        
        % initialise peak results
        numbPeaksPart = 0;      % number of peaks found in this partition
        resPeaksPart = zeros(0, 5);
        
        for nZone = 1 : nSDNegZones
%             % first check that the ZD value with the minimum SD is actually
%             % above the threshold set for a valid signal 
%             [temp, indMinSD] = min( ...
%                 smdSD(indSD_NegZoneSta(nZone):indSD_NegZoneEnd(nZone)) ); %#ok<ASGLU>
%             indPeakMax = indSD_NegZoneSta(nZone)+indMinSD-1; 

            % get the ZD values within the negative zone
            smdZDinNegZone = smdZD(indSD_NegZoneSta(nZone):indSD_NegZoneEnd(nZone));
            
            % the peak is the maximum ZD value
            [temp, indMaxZD] = max(smdZDinNegZone); %#ok<ASGLU>
            indPeakMax = indSD_NegZoneSta(nZone)+indMaxZD-1; 
            
            % also check the first derivative in this zone is above the
            % threshold for FD  -- this bit is not in the ref. paper
            smdFDinNegZone = smdFD(indSD_NegZoneSta(nZone):indSD_NegZoneEnd(nZone));

            % check to make sure that the previous peak doesn't end after
            % the current peak maximum. if so, then we need to do something
            % about it in the peak detection loop. 
            % -- this isn't in the ref'd paper. 
            prevPeakEndB4CurrPeak = false;
            if numbPeaksPart == 0
                prevPeakEndB4CurrPeak = true;
            elseif timePart(indPeakMax) > resPeaksPart(numbPeaksPart, colPeakEnd)
                prevPeakEndB4CurrPeak = true;
            end
            
            % conditions for considering this zone: 
            %  -- the smoothed signal in this zone above thrZD, and 
            %  -- any abolute value of FD in this zone is above thrFD.
            % else go to the next zone
            %  -- if the length of the zone is greater than one (not in
            % ref'd paper)
            if (any(smdZDinNegZone > thrSig)) && (any(abs(smdFDinNegZone) > thrFD)) ...
                    && (length(smdZDinNegZone) > 1)
                
                % increment the count of peak numbers in this partition
                numbPeaksPart = numbPeaksPart + 1;
                
                % we can write down the time index for the peak maximum and
                % its value 
                resPeaksPart(numbPeaksPart, colPeakMax) = timePart(indPeakMax);
                resPeaksPart(numbPeaksPart, colPeakHeight) = smdZD(indPeakMax);
                
                % if previous peak's end is after the current peak maximum,
                % then we need to change the peak ending time for the
                % previous peak
                if ~prevPeakEndB4CurrPeak
                    % get the index for the end of the previous negative SD
                    % zone. 
                    indEndPrevSD = indSD_NegZoneEnd(nZone-1);
                    % starting from this index, we search in the third
                    % derivative for the change in sign between two time
                    % indices. 
                    indTemp = indEndPrevSD;
                    while true
                        if indTemp > indSD_NegZoneSta(nZone)
                            indEndPrevSD = indTemp - 1;
                            break;
                        elseif diff(smdTD(indTemp:indTemp+1) > 0) == 0
                            indTemp = indTemp+1;
                        else
                            indEndPrevSD = indTemp+1;
                            break;
                        end
                    end
                    resPeaksPart(numbPeaksPart-1, colPeakEnd) = ...
                        timePart(indEndPrevSD);
                end
                
                % search in the first derivative for start of peak
                indFD_startOfPeak = indSD_NegZoneSta(nZone);
                while true
                    % if the current value is at the very first index of
                    % the partition, then we've found peak
                    if indFD_startOfPeak == 1
                        break;
                    end
                    % if the current value and the previous value have
                    % different signs, then we've found start of peak 
                    if changeInSignFD(indFD_startOfPeak-1) 
                        break;
                    end
                    % use absolute value of the FD, for those cases where
                    % the FD at the beginning of the zone isn't positive
                    if abs(smdFD(indFD_startOfPeak)) > thrFD
                        % go to the previous FD value
                        indFD_startOfPeak = indFD_startOfPeak - 1;
                    else
                        % found the start of the peak
                        break;
                    end
                    % also check that the previous value is: 
                    % -- the very first value of the partition, or 
                    % -- the same as the last peak's end.
                    % if so, start found.
                    if indFD_startOfPeak == 1
                        break;
                    elseif numbPeaksPart>1 && ...
                            timePart(indFD_startOfPeak) == resPeaksPart(numbPeaksPart-1, colPeakEnd)
                        break;
                    end
                end
                % write down the start of the peak 
                resPeaksPart(numbPeaksPart, colPeakStart) = ...
                    timePart(indFD_startOfPeak);
                
                % search in the first derivative for the end of the peak
                indFD_endOfPeak = indSD_NegZoneEnd(nZone);
                while true
                    % also check that the current value is not the very 
                    % last value of the partition, if so, end found.
                    if indFD_endOfPeak == NDataPart
                        break;
                    end
                    % if the current value and the next value have
                    % different signs, then we've found peak end
                    if changeInSignFD(indFD_endOfPeak) 
                        break;
                    end
                    % absolute value of the FD above the threshold, for
                    % those cases that the FD at the end of the zone isn't
                    % negative.
                    if abs(smdFD(indFD_endOfPeak)) > thrFD
                        % go to the next FD value
                        indFD_endOfPeak = indFD_endOfPeak + 1;
                    else
                        % found the end of the peak
                        break;
                    end
                end
                % write down the end of the peak 
                resPeaksPart(numbPeaksPart, colPeakEnd) = ...
                    timePart(indFD_endOfPeak);
                
                
                %==================
                % STEP 3e. 
                % check if the SD in the negative zone has more than 1 dip,
                % if so then we have a case of strong co-elusion, and
                % require a split of this peak into the seperate peaks
                
                % first calculate the number of change of sign in the
                % third derivative
                % i. get the third derivative values in the negative zone
                smdTDinNegZone = smdTD(indSD_NegZoneSta(nZone):indSD_NegZoneEnd(nZone));
                % ii. get logical index of TD values greater than 0 (+ve
                % and -ve values of TD), any difference between consequent
                % indices that's not 0 shows there's a change in sign
                indChangeInSignTD = find(diff(smdTDinNegZone > 0));
                % iii. count the number of change in sign
                nChangeInSignTD = length(indChangeInSignTD);
                % iv. calculate the number of minimum in the SD
                if mod(nChangeInSignTD, 2)
                    % if n is odd
                    nDipsInSD = (nChangeInSignTD + 1) / 2;
                else
                    % if n is even
                    % nb. don't know how this could happen but I'll leave
                    % this in here
                    nDipsInSD = nChangeInSignTD / 2;
                end
                
                % If there are more than one dips in the SD -ve zone, then
                % we'll need to split the current peak
                if nDipsInSD > 1 && splitCoeultion
                    % set aside another matrix for these results
                    multiPeakZoneRes = zeros(nDipsInSD, 5);
                    % the start and end of the peak region now goes into
                    % the start of peak 1 and end of last peak
                    multiPeakZoneRes(1, colPeakStart) = ...
                        resPeaksPart(numbPeaksPart, colPeakStart);
                    multiPeakZoneRes(nDipsInSD, colPeakEnd) = ...
                        resPeaksPart(numbPeaksPart, colPeakEnd);
                    % keep track the indices for the start and end
                    mpi1 = zeros(nDipsInSD, 1);
                    mpi2 = zeros(nDipsInSD, 1);
                    mpi1(1) = indFD_startOfPeak;
                    mpi2(end) = indFD_endOfPeak;
                                        
                    % v. find all the minimum SD in the -ve zone, these
                    % will be the peak maxima, this is when the TD change
                    % from -ve to +ve
                    indMinsSD = find(diff(smdTDinNegZone > 0) == 1);
                    indPeakMaxInMultiZone = indSD_NegZoneSta(nZone)+indMinsSD-1; 
                    multiPeakZoneRes(:, colPeakMax) = timePart(indPeakMaxInMultiZone)';
                    multiPeakZoneRes(:, colPeakHeight) = smdZD(indPeakMaxInMultiZone)';
                    
                    % vi. when the TD change from +ve to -ve, that's when 
                    % the split should be
                    indMaxsSD = find(diff(smdTDinNegZone > 0) == -1);
                    indPeakSplit = indSD_NegZoneSta(nZone)+indMaxsSD-1;
                    if ~mod(nChangeInSignTD, 2)
                        % if n is even, then we need to make sure the split
                        % happens after the peak
                        indPeakSplit_indActual = ...
                            find(indPeakSplit > indPeakMaxInMultiZone(1));
                        indPeakSplit = indPeakSplit(indPeakSplit_indActual);
                    end
                    multiPeakZoneRes(1:nDipsInSD-1, colPeakEnd) = ...
                        timePart(indPeakSplit)';
                    multiPeakZoneRes(2:nDipsInSD, colPeakStart) = ...
                        timePart(indPeakSplit)';
                    mpi1(2:nDipsInSD) = indPeakSplit;
                    mpi2(1:nDipsInSD-1) = indPeakSplit;
                    
                    %==================
                    % STEP 3f. calculate the area under the curve
                    for nDip = 1 : nDipsInSD
                        multiPeakZoneRes(nDip, colPeakArea) = ...
                            calcAUC( smdZD( mpi1(nDip):mpi2(nDip) ), ...
                            timePart( mpi1(nDip):mpi2(nDip) ) );
                    end
                    % write the whole temporary matrix into the matrix for
                    % this partition's result
                    resPeaksPart(numbPeaksPart:numbPeaksPart+nDipsInSD-1, :) = ...
                        multiPeakZoneRes;
                    
                    % update the number of peaks found in this partition
                    numbPeaksPart = numbPeaksPart + nDipsInSD - 1;

                else
                    % if there's only one dip
                    
                    %==================
                    % STEP 3f. calculate the area under the curve
                    % first check that the peak is not just one point, if
                    % it is, then delete the current row, and decrement the
                    % count
                    if length(smdZD(indFD_startOfPeak:indFD_endOfPeak)) > 1
                        areaPeak = calcAUC(smdZD(indFD_startOfPeak:indFD_endOfPeak), ...
                            timePart(indFD_startOfPeak:indFD_endOfPeak));
                        % write down the area under the peak
                        resPeaksPart(numbPeaksPart, colPeakArea) = areaPeak;
                    else
                        resPeaksPart(numbPeaksPart, :) = [];
                        numbPeaksPart = numbPeaksPart - 1;
                    end
                end
                    

            end % end if : the -ve zone is valid for peak
            
        end % end for : all the -ve zones in this partition
        
        % plot all the peaks in this partition on the figure
        if showPlot
            figure(hFigPart);
            propsMarker = {'MarkerSize', 12};
            plot(resPeaksPart(:, colPeakStart), zeros(numbPeaksPart, 1), ...
                'm^', propsMarker{:}, 'MarkerFaceColor', 'm');
            plot(resPeaksPart(:, colPeakEnd), zeros(numbPeaksPart, 1), ...
                'rv', propsMarker{:}, 'MarkerFaceColor', 'r');
            plot(resPeaksPart(:, colPeakMax), zeros(numbPeaksPart, 1), ...
                'go', propsMarker{:}, 'MarkerFaceColor', 'g');
        end
        
        % Output number of peaks found in the partition
        fprintf('Partition number %d, number of peaks found: %d.   ', ...
            partNumb, numbPeaksPart);
        if showPlot
            fprintf('Press any key to continue. \n');
            pause;
        else
            fprintf('\n');
        end
        
        % update peaks for the full chromatogram
        NPeaksTotal = NPeaksTotal + numbPeaksPart;
        resPeaks = [resPeaks; resPeaksPart];
        
    end     % end of peak detection for each partition
end     % analysed all data presented. 


if showPlot
    % plot all the peaks on the full signal
    figure(hFig);
    clf;
    hold on;
    plot(timesteps, gcData, props{:});
    set(gca, 'YLim', [0 yULim]);
    grid on;
    propsMarker = {'MarkerSize', 8};
    plot(resPeaks(:, colPeakStart), zeros(NPeaksTotal, 1), ...
        'm^', propsMarker{:}, 'MarkerFaceColor', 'm');
    plot(resPeaks(:, colPeakEnd), zeros(NPeaksTotal, 1), ...
        'rv', propsMarker{:}, 'MarkerFaceColor', 'r');
    plot(resPeaks(:, colPeakMax), zeros(NPeaksTotal, 1), ...
        'go', propsMarker{:}, 'MarkerFaceColor', 'g');
end

fprintf('------------------\n');
fprintf('Total number of peaks found: %d. \n', NPeaksTotal);


%%-------------------------------------
% subfunction to calculate the noise of a signal. 
%   The noise is the median of all abolute differences between the signal
% at i and the mean of its immediate neighbours i-1 and i+1
function noiseVal = calcNoise(signal)

neighbour1 = signal(1:end-2);
neighbour2 = signal(3:end);
meanNeighbours = (neighbour1 + neighbour2) / 2;
diffWithNeighbour = signal(2:end-1) - meanNeighbours;
noiseVal = median(abs(diffWithNeighbour));


%%-------------------------------------
% subfunction to calculate the area under the curve of a signal. 
%   The area is calculated using the Trapezoidal rule since some of the
% peaks are very sharp to use first principle. 
%   The inputs are the part of the signal with the peak, and the elution
% time or the timesteps corresponding to the signal
function resArea = calcAUC(signal, etime)

val1 = signal(1:end-1);
val2 = signal(2:end);
incTime = diff(etime);
resArea = sum((val1 + val2) .* incTime ./ 2);


%%-------------------------------------
% subfunction to partition the whole signal
function [partStIndex, partEndIndex] = partitionSignal(smdInitial, ...
    thrSig, N, splitSignalSize)

nPart = ceil(N/splitSignalSize);
partStIndex = zeros(nPart, 1);
k = 1;
partStIndex(k) = 1;
flgPartition = 1;
while flgPartition
    k = k + 1;
    % get to the next start point
    stNext = partStIndex(k-1) + splitSignalSize;
    % check if this point and it's surrounding 10 points are all below
    % the threshold, if not, then go to the next point along the line
    % and search until a section is found.
    surrPoints = 10;
    while true
        %secVal = smdInitial(stNext-surrPoints : stNext+surrPoints);
        if smdInitial(stNext-surrPoints : stNext+surrPoints) < thrSig
            break;
        else
            stNext = stNext + 1;
        end
        % make sure that the index stNext is not near the end of the
        % data
        if stNext+surrPoints >= N
            stNext = 0;
            break;
        end
    end
    partStIndex(k) = stNext;
    % check if the number of points reminding is smaller than
    % splitSignalSize or not, exit while loop if true
    if N - stNext - surrPoints <= splitSignalSize || stNext == 0
        flgPartition = 0;
        break;
    end
end
% make sure there's no zeros in the indices, delete these if so
partStIndex(find(partStIndex == 0)) = []; %#ok<FNDSB>

% the indices for the end of the partitions
partEndIndex = [partStIndex(2:end)-1; N];
