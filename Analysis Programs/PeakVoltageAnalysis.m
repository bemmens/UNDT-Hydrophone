%% Peak Voltage Analysis
% For reference: scanData = [x,y,samples,chanel];
clear all
%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'PinCyclePin4';
path = strcat(folder_path,file_name,'.mat');
load(path)

%% Extract Data
Vpk = squeeze(max(scanData(:,:,:,1),[], 3)); % peak voltage

%% Plot
%figure(1)
%(xs,ys,Vpk)
%colorbar
%xlabel('x (mm)')
%ylabel('y (mm)')
%title('Volts')

%% To MPa
mVperMPa = 170.12; % CHECK
MPa = Vpk*1e3/mVperMPa; 

figure(2)
imagesc(raster.xs,flip(raster.ys),rot90(MPa))
axis image;
set(gca,'YDir','normal') % rot90, flip to get stage coords to match image
colorbar
xlabel('x (mm)')
ylabel('y (mm)')
title('MPa')

%% Histogram
flatMpa = reshape(MPa,[],1);
figure(3)
hist(flatMpa)
xlabel('MPa')
ylabel('Count')
title('Peak Pressure Distribution')