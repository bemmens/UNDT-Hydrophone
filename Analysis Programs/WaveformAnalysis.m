%% Plot waveforms

%% Load Data
% Consider clearing workspace to relieve RAM
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'test';
path = strcat(folder_path,file_name,'.mat');
load(path)

%% Plot
plot(squeeze(scanData(1,1,1:1e6,:)))