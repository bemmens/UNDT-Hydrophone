%% Peak Voltage Analysis
% For reference: scanData = [x,y,samples,chanel];

%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'test';
path = strcat(folder_path,file_name,'.mat');
load(path)

%% Extract Data
Vpk = squeeze(max(scanData(:,:,:,1),[], 3));

%% Plot
imagesc(xs,ys,Vpk)
colorbar
xlabel('x (mm)')
ylabel('y (mm)')