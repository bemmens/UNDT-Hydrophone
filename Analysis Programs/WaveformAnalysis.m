%% Plot waveforms
clear all
%% Load Data
% Consider clearing workspace to relieve RAM
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'MatchTest16';
path = strcat(folder_path,file_name,'.mat');
load(path)

%% Slect Data
size(scanData)
x_index = 5;
y_index = 5;
data = squeeze(scanData(x_index,y_index,:,1));

%% Plot
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us
figure(1)
plot(t,data)
xlim([0,0.2*1000])
hold on
x = raster.xs(x_index);
y = raster.ys(y_index);
title(strcat('Waveform at [x,y]=[',string(x),',',string(y),'] mm'))
xlabel('Time [us]');
ylabel('Amplitude [V]');
tx_rcvd = 32; % us
%xline(tx_rcvd)
hold off

%% Calculate travel times

c = 1481; % m/s speed of sound in water
d = c*tx_rcvd/1e3; % mm
x = c*t/1e3;
disp(strcat('Distance to sensor:', string(d),'mm.'))

%% Spectrogram
figure(5)
spectrogram(data,3e2,[],1e3,scpSettings.SampleFrequency,'yaxis');
ylim([0,5])
colorbar('off')
title('Bandpassed Waveform Spectrogram')

%% Band pass filter
bpass = bandpass(data,[0.5e6,1.5e6],scpSettings.SampleFrequency); % 1 +/- 0.1MHz
figure(4)
%plot(t,bpass)
plot(x,bpass)
xlim([0,300])
%xline([32])
%xlabel('Time [us]')
xlabel('Distance [mm]')
ylabel('Voltage [V]')
title('Single Cycle Pulse with Bandpass Filter: 1+/-0.5MHz')

