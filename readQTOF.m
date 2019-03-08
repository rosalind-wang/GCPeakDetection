function readQTOF()
% code to read Q-TOF xml data file. 

specBegin = 35;
specEnd = 350;
spectra = specBegin : specEnd;
nSpectra = length(spectra);

mainDir = '/data/Malaria/';
fileName = 'S025 D5 AM R1.mzdata.xml';
fileWriteName = 'test.mzdata.csv';

% files = dir(sprintf('%s/*.xml', mainDir));
% fileName = files(nfile).name;
% fileWriteName = [fileName(1:end-4), '.csv'];

fid = fopen([mainDir, fileName], 'rt');
fidWrite = fopen([mainDir, fileWriteName], 'wt');
fprintf(fidWrite, 'nan,');
for i = 1 : nSpectra
    fprintf(fidWrite, '%d,', spectra(i));
end
fprintf(fidWrite, '\n');
fclose(fidWrite);

fidWrite = fopen([mainDir, fileWriteName], 'at');
count = 0;
while ~feof(fid)
    inputLine = fgetl(fid);
    
    % look for "spectrumList count", this will tell us how many time steps
    % we've got for this data file
    % not necessary
%     if strfind(inputLine, 'spectrumList count')
%         indQuote = strfind(inputLine, '"');
%         nTimesteps = str2num(inputLine( indQuote(1)+1 : indQuote(2)-1 ));
%     end
    
    % if we find the pattern "spectrum id" in the line, then we're starting
    % the data for a spectrum line, otherwise, move on
    if strfind(inputLine, 'spectrum id')
        count = count + 1;
        if mod(count, 100) == 0
            fprintf('Reading spectrum line %d\n', count);
        end
        
        [timestamp, mzValues, intValues] = readSpectrumLine(fid);
        
        % create a vector to put the respective intensity values in
        intensities = zeros(nSpectra, 1);
        mzValuesRound = round(mzValues);
        for j = 1 : length(mzValuesRound)
            intensities(mzValuesRound(j)-specBegin+1) = ...
                intensities(mzValuesRound(j)-specBegin+1) + intValues(j);
        end
        
        % write the values into the file
        fprintf(fidWrite, '%.3f,', timestamp);
        for i = 1 : nSpectra
            fprintf(fidWrite, '%d,', intensities(i));
        end
        fprintf(fidWrite, '\n');
    end
        
end

fclose(fid);
fclose(fidWrite);


%%%%%%%
% sub-function for read the data for each spectrum line

function [timestamp, mzValues, intValues] = readSpectrumLine(fid)

flagReadEnd = 0;

while ~flagReadEnd
    inputLine = fgetl(fid);
    
    % find the entry for the time in minutes
    if strfind(inputLine, 'TimeInMinutes')
        indValueTag = strfind(inputLine, 'value=');
        tempLine = inputLine(indValueTag:end);
        indQuote = strfind(tempLine, '"');
        timestamp = str2num(tempLine( indQuote(1)+1 : indQuote(2)-1 ));
    end
    
    % find the entry for the mz values
    if strfind(inputLine, '<mzArrayBinary>')
        inputLine = fgetl(fid);
        indRAngle = strfind(inputLine, '>');
        indLAngle = strfind(inputLine, '<');
        dataLine = inputLine(indRAngle(1)+1 : indLAngle(2)-1);
        
        mzValues = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(uint8(dataLine)), 'double');
    end
    
    % find the entry for the measurements
    if strfind(inputLine, '<intenArrayBinary>')
        inputLine = fgetl(fid);
        indRAngle = strfind(inputLine, '>');
        indLAngle = strfind(inputLine, '<');
        dataLine = inputLine(indRAngle(1)+1 : indLAngle(2)-1);
        
        intValues = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(uint8(dataLine)), 'single');
        
        % set flag for read this particular time step
        flagReadEnd = 1;
    end
    
end