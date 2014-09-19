%% Introduction
% For this project we need to plot the following MIMO parameters:
%   - Capacity loss (C_loss)
%   - Total Active Refelection Coefficient (TARC)
%   - Correlation factor/coefficient (rho)
%   - Mean Effective Gain (MEG)

clear all;

%% Import S-parameters and Power Gain Patterns
% Import the file
isolation = false;
if isolation
    fileNameReal = 'F:\Grad_School\EE740\01_MIMO_Configuration_Mesh-Sparams_real.csv';
    fileNameImag = 'F:\Grad_School\EE740\01_MIMO_Configuration_Mesh-Sparams_imag.csv';
    fileNameGain1p8 = 'F:\Grad_School\EE740\01_MIMO_Configuration_Mesh-gain1.8GHz.csv';
    fileNameGain5p5 = 'F:\Grad_School\EE740\01_MIMO_Configuration_Mesh-gain5.5GHz.csv';
else
    fileNameReal = 'F:\Grad_School\EE740\02_MIMO_NoIsolation-Sparams_real.csv';
    fileNameImag = 'F:\Grad_School\EE740\02_MIMO_NoIsolation-Sparams_imag.csv';
    fileNameGain1p8 = 'F:\Grad_School\EE740\02_MIMO_NoIsolation-gain1.8GHz.csv';
    fileNameGain5p5 = 'F:\Grad_School\EE740\02_MIMO_NoIsolation-gain5.5GHz.csv';
end

importedReal = importdata(fileNameReal);
importedImag = importdata(fileNameImag);
importedGain1p8 = importdata(fileNameGain1p8);
importedGain5p5 = importdata(fileNameGain5p5);

% Create new variables in the base workspace from those fields.
vars = fieldnames(importedReal);
for i = 1:length(vars)
    assignin('base', vars{i}, importedReal.(vars{i}));
end
vars = fieldnames(importedImag);
for i = 1:length(vars)
    assignin('base', vars{i}, importedImag.(vars{i}));
end
vars = fieldnames(importedImag);
for i = 1:length(vars)
    assignin('base', vars{i}, importedGain1p8.(vars{i}));
end
vars = fieldnames(importedImag);
for i = 1:length(vars)
    assignin('base', vars{i}, importedGain5p5.(vars{i}));
end

freq = importedReal.data(:,1); % Doesn't matter which it takes it from
S11 = complex(importedReal.data(:,2),importedImag.data(:,2));
S12 = complex(importedReal.data(:,3),importedImag.data(:,3));
S21 = complex(importedReal.data(:,4),importedImag.data(:,4));
S22 = complex(importedReal.data(:,5),importedImag.data(:,5));

bands = [1.7, 1.9, 4.5, 5.9]; % From isolated version
gainFreqs = [1.8, 5.5];
phi = importedGain1p8.data(:,1);
gainPhi(1,:) = importedGain1p8.data(:,2);
gainTheta(1,:) = importedGain1p8.data(:,3);

gainPhi(2,:) = importedGain5p5.data(:,2);
gainTheta(2,:) = importedGain5p5.data(:,3);

% Plots
figs.sParams = figure;
plot(freq, 20*log10(abs(S11)));
hold on;
plot(freq, 20*log10(abs(S22)),'r');
plot(freq, 20*log10(abs(S12)),'g');
plot(freq, 20*log10(abs(S21)),'black');
ylims = get(gca, 'Ylim');
plot([bands(1), bands(1)],[-100, 100], '--black');
plot([bands(2), bands(2)],[-100, 100], '--black');
plot([bands(3), bands(3)],[-100, 100], '--black');
plot([bands(4), bands(4)],[-100, 100], '--black');
plot([freq(1) freq(end)], [-10, -10], '--green');
set(gca,'Ylim',ylims);
title('S-parameters (dB)');
legend('S11','S12','S21','S22','Location','best');
xlabel('Frequency (GHz)');
ylabel('S-params (dB)');

figs.gainLow = figure;
plot(phi, 10*log10(gainPhi(1,:)));
hold on;
plot(phi, 10*log10(gainTheta(1,:)),'r');
plot([phi(1) phi(end)], [0 0], '--g');
title('Low band - 1.8 GHz');
legend('gainPhi(dB)', 'gainTheta(dB)','Location','best');
xlabel('Phi (degrees)');
ylabel('Realized Gain (dB)');

figs.gainHigh = figure;
plot(phi, 10*log10(gainPhi(2,:)));
hold on;
plot(phi, 10*log10(gainTheta(2,:)),'r');
plot([phi(1) phi(end)], [0 0], '--g');
title('High band - 5.5 GHz');
legend('gainPhi(dB)', 'gainTheta(dB)','Location','best');
xlabel('Phi (degrees)');
ylabel('Realized Gain (dB)');

%% Capacity Loss
% Calculation
psiR = zeros(2,2);
CL = zeros(1,length(freq));
for i = 1:length(freq)
    psiR(1,1) = 1-(abs(S11(i)).^2 + abs(S12(i)).^2);
    psiR(1,2) = -(conj(S11(i)).*S12(i) + conj(S21(i)).*S22(i));
    psiR(2,1) = -(conj(S22(i)).*S21(i) + conj(S12(i)).*S11(i));
    psiR(2,2) = 1-(abs(S22(i)).^2 + abs(S21(i)).^2);
    CL(i) = -log2(det(psiR));
end

% Plots
 % Note: Warning about imaginary part can be ignored - it looks to be on 
 % order of 5e-21.
figs.cl = figure;
plot(freq,CL);
hold on;
ylims = get(gca, 'Ylim');
plot([bands(1), bands(1)],[-100, 100], '--black');
plot([bands(2), bands(2)],[-100, 100], '--black');
plot([bands(3), bands(3)],[-100, 100], '--black');
plot([bands(4), bands(4)],[-100, 100], '--black');
set(gca,'Ylim',ylims);
title('Capacity Loss');
xlabel('Frequency (GHz)');
ylabel('Capacity loss');

%% TARC
% Calculation
thetaVals = [0, 30, 60, 90, 180];
tarc = zeros(length(thetaVals),length(freq));
for i = 1:length(thetaVals)
    for j = 1:length(freq)
        theta = thetaVals(i);
        tarc(i,j) = sqrt(abs(S11(j)+S12(j).*exp(complex(0,1)*theta))^2 +...
            abs(S21(j)+S22(j).*exp(complex(0,1)*theta))^2)/sqrt(2);
    end
end

% Plots
figs.tarc = figure;
plot(freq, 20*log10(tarc));
hold on;
ylims = get(gca, 'Ylim');
plot([bands(1), bands(1)],[-100, 100], '--black');
plot([bands(2), bands(2)],[-100, 100], '--black');
plot([bands(3), bands(3)],[-100, 100], '--black');
plot([bands(4), bands(4)],[-100, 100], '--black');
plot([freq(1) freq(end)], [-6, -6], '--green');
set(gca,'Ylim',ylims);
legendStr = [];
for i = 1:length(thetaVals)
    legendStr{i} = sprintf('theta = %d',thetaVals(i));
end
title('Total Active Reflection Coefficient');
legend(legendStr,'Location','best');
xlabel('Frequency (GHz)');
ylabel('TARC (dB)');

%% Correlation factor
rho = abs(conj(S11).*S12+conj(S21).*S22).^2./(1-abs(S11).^2-abs(S21).^2)./(1-abs(S22).^2-abs(S12).^2);

% Plots
figs.rho = figure;
plot(freq, rho);
hold on;
ylims = get(gca, 'Ylim');
plot([bands(1), bands(1)],[-100, 100], '--black');
plot([bands(2), bands(2)],[-100, 100], '--black');
plot([bands(3), bands(3)],[-100, 100], '--black');
plot([bands(4), bands(4)],[-100, 100], '--black');
plot([freq(1) freq(end)], [0.3 0.3], '--g');
set(gca,'Ylim',ylims);
title('Correlation Factor');
xlabel('Frequency (GHz)');
ylabel('Correlation Factor');

%% Mean Effective Gain (MEG)
xpdVals = [1, 0.5];
deltaPhi = 2/180*pi; % radians
meg = zeros(length(xpdVals),length(gainFreqs));
for i = 1:length(xpdVals)
    for j = 1:2 % number of frequencies for gains
        xpd = xpdVals(i);
        meg(i,j) = (1/2/pi)*sum(xpd/(1+xpd)*gainTheta(j,:)+1/(1+xpd)*gainPhi(j,:))*deltaPhi;
    end
end

fprintf('Mean Effective Gain Table (dB)\r\n');
fprintf('Freq(GHz)\t\tXPD = 1\t\tXPD = 0.5\r\n');
fprintf('%.1f\t\t\t%.3f\t\t%.3f\r\n',gainFreqs(1),20*log10(meg(1,1)),20*log10(meg(1,2)));
fprintf('%.1f\t\t\t%.3f\t\t%.3f\r\n',gainFreqs(2),20*log10(meg(2,1)),20*log10(meg(2,2)));

%% Clean up
names = fieldnames(figs);
for i = 1:length(names)
    close(figs.(names{i}));
end