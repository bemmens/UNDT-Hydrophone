%% Plot 2dv waveforms
clear all
%% Load Data
% Consider clearing workspace to relieve RAM
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'MatchTest16';
path = strcat(folder_path,file_name,'.mat');
load(path)

in = scanData(:,:,20000,1);
imagesc(raster.xs,flip(raster.ys),rot90(in))
