%% Plot waveforms

%% Load Data
% Consider clearing workspace to relieve RAM
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'Echocheck2';
path = strcat(folder_path,file_name,'.mat');
load(path)

%% Slect Data
size(scanData)
x_index = 3;
y_index = 3;
data = squeeze(scanData(x_index,y_index,:,:));

%% Calculate travel times

c = 1481; % m/s speed of sound in water
x = c*t/1e3;

%% 3Plot
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us
figure(1)
%plot(t,squeeze(scanData(x_index,y_index,:,1)))
plot(x,squeeze(scanData(x_index,y_index,:,1)))
%xlim([0,0.2*1000])
xlim([0,350])
ylim([-0.3,0.3])
xline([20,60,100,140,180,220]-4)
hold on
title(strcat('Pulse With F28 - xline 40mm spacing'))
%xlabel('Time [us]');
xlabel('Distance [mm]');
ylabel('Amplitude [V]');
%xline(tx_rcvd)
hold off

%% Load Data
% Consider clearing workspace to relieve RAM
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'PZTDiscBeamform5';
path = strcat(folder_path,file_name,'.mat');
load(path)

%% Slect Data
size(scanData)
x_index = 50;
y_index = 50;
data = squeeze(scanData(x_index,y_index,:,:));

%% Plot
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us
figure(2)
%plot(t,squeeze(scanData(x_index,y_index,:,1)))
plot(x,squeeze(scanData(x_index,y_index,:,1)))
%xlim([0,0.2*1000])
xlim([0,350])
ylim([-0.3,0.3])
title(strcat('Pusle Without F28 - xline 40mm spacing'))
%xlabel('Time [us]');
xlabel('Distance [mm]');
ylabel('Amplitude [V]');
%xline(tx_rcvd)
xline([20,60,100,140,180,220]+1)
hold off



