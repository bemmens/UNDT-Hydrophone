%% Plot waveforms

%% Load Data
% Consider clearing workspace to relieve RAM
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'PinCyclePin4';
path = strcat(folder_path,file_name,'.mat');
load(path)

%% Slect Data
size(scanData)
x_index = 13;
y_index = 13;
data = squeeze(scanData(x_index,y_index,:,:));

%% Plot
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us
figure(1)
plot(t,squeeze(scanData(x_index,y_index,:,1)))
xlim([25,60])
hold on
x = raster.xs(x_index);
y = raster.ys(y_index);
title(strcat('Waveform at [x,y]=[',string(x),',',string(y),'] mm'))
xlabel('Time [us]');
ylabel('Amplitude [V]');
tx_rcvd = 32; % us
xline(tx_rcvd)
hold off

%% Calculate travel times

c = 1481; % m/s speed of sound in water
d = c*tx_rcvd/1e3; % mm
disp(strcat('Distance to sensor:', string(d),'mm.'))

