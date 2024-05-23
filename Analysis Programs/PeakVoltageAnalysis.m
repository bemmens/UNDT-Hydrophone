%% Peak Voltage Analysis
% For reference: scanData = [x,y,samples,chanel];

%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'test';
path = strcat(folder_path,file_name,'.mat');
load(path)

%% Extract Data
Vpk = squeeze(max(scanData(:,:,:,1),[], 3)); % peak voltage
mVperMPa = 170; % CHECK
MPa = Vpk*1e3/mVperMPa; 

%% Plot
figure(1)
imagesc(xs,ys,Vpk)
colorbar
xlabel('x (mm)')
ylabel('y (mm)')
title('Volts')

figure(2)
imagesc(xs,ys,MPa)
colorbar
xlabel('x (mm)')
ylabel('y (mm)')
title('MPa')