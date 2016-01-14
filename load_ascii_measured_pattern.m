%% User input variables
% Enter all needed variables here to convert to gain
% Note that frequency is read from the files

dist = 1.5; % Distance from feed to antenna in meters
Gt = 6; % Gain of transmitter in dB at one freq

%% Read in data from file

[fileName, pathName] = uigetfile('*.txt', 'Pick your ASCII file');
fullFileName = fullfile(pathName, fileName);
fid = fopen(fullFileName,'r');

if fid == -1
    error('Problem opening the file');
end

configCount = 0;
valueCount = 0;
line = fgetl(fid);
while ischar(line)
    pieces = regexp(line,'\t','split');
    if isnan(str2double(pieces{1}))
        switch pieces{1}
            case 'File Name:'
                configCount = configCount + 1;
                valueCount = 0;
            case 'Feed/Prob.'
                % Get the angle of rotation for feed
                feedRot(configCount) = str2double(pieces{2});
            case 'Frequency'
                % Get the frequency of experiment
                subPieces = regexp(pieces{2}, ' ', 'split');
                freq(configCount) = str2double(subPieces{1});
            case 'AUT_Roll  ' % Two spaces after AUT_Roll
                % Get angle of rotation for the AUT
                autRoll(configCount) = str2double(pieces{2});
            case 'Ch.'
                % Not sure
            case 'beam'
                % Not sure
            case 'Switch'
                % S21? Is this needed?
            case 'Azimuth [dB]'
                % Just labels, might be needed if units change
            otherwise
                % Blank line most likely
        end
    else
        % It was converted properly => Beginning of data
        valueCount = valueCount + 1;
        azimuth(valueCount, configCount) = str2double(pieces{1});
        amplitude(valueCount, configCount) = str2double(pieces{2});
        phase(valueCount, configCount) = str2double(pieces{3});
    end
    line = fgetl(fid);
end

fclose(fid);
%% Update user on what was read
fprintf('Read file: %s\n', fullFileName);
fprintf('Read in %d configurations\n', configCount);
fprintf('Each configuration had %d values\n', valueCount); % Same count?
fprintf('AUT_Roll:\t\t');
fprintf('%.2f, ', autRoll(:));
fprintf('\n');
fprintf('Feed Rotation:\t');
fprintf('%.2f, ', feedRot(:));
fprintf('\n');

%% Convert amplitudes to gains using Friss Transmission Equation
if length(unique(freq)) ~= 1
    error(['Different frequencies in same file?',...
           'This is not supported yet.',...
           'Break up the files by frequency and this MIGHT work.']);
end

% =-GT +AMP  +(20*LOG(12.56*Radius / WAVELENGTH)) 

wavelength = 0.3 / freq(1);
gain = -Gt + amplitude + (20*log10(12.56*dist./wavelength));

%% Save the file, with headers
[p,n,e,v] = fileparts(fullFileName);
newName = sprintf('%s_GainOutput.csv', n);
saveFileName = fullfile(p,newName);

%dB(RealizedGainPhi) [] - Freq='0.725GHz' Phi='0deg'
header = {'Theta [deg]',...
          sprintf('dB(RealizedGainPhi) [] - Freq=''%.3fGHz'' Phi=''0deg''', freq(1)),...
          sprintf('dB(RealizedGainPhi) [] - Freq=''%.3fGHz'' Phi=''0deg''', freq(1)),...
          sprintf('dB(RealizedGainPhi) [] - Freq=''%.3fGHz'' Phi=''0deg''', freq(1)),...
          sprintf('dB(RealizedGainPhi) [] - Freq=''%.3fGHz'' Phi=''0deg''', freq(1))};
xlswrite(saveFileName, header) % by defualt starts from A1
xlswrite(saveFileName, [azimuth(:,1) gain],'sheet1','A2') % array under the header.

fprintf('Saved data to %s\n', saveFileName);
fprintf('COMPLETED\n');
fprintf('**************************\n\n');