%% Peak Voltage Analysis
% For reference: scanData = [x,y,z,samples,channel];
clear all

%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = '3DTest';
path = strcat(folder_path,file_name,'.mat');
load(path)
disp('Data Timestamp:')
disp(scpSettings.timestamp)
disp('Data Size:')
disp(size(scanData))

raster.relxs = raster.xs - raster.home(1);
raster.relys = raster.ys - raster.home(2);
raster.relzs = raster.zs - raster.home(3);

%% Extract Peak Voltages
pkrange = [10,15]; % us - time range to look for peak 
pkrangeidx = pkrange*scpSettings.SampleFrequency/1e6; % corresponding array index
Vpk = squeeze(max(scanData(:,:,:,pkrangeidx(1):pkrangeidx(2),1),[],4));

%% Check Waveform
x_index = 1;
y_index = 2;
z_index = 3;
pks = find(scanData(x_index,y_index,z_index,:,1) == Vpk(x_index,y_index,z_index));
wvfmData = squeeze(scanData(x_index,y_index,z_index,:,1));
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us
figure(1)
plot(t,wvfmData)
xlim([0,200])
hold on
x = raster.relxs(x_index);
y = raster.relys(y_index);
z = raster.relys(z_index);
title(strcat('Waveform at [x,y]=[',string(x),',',string(y),',',string(z),'] mm'))
xlabel('Time [us]');
ylabel('Amplitude [V]');
hold off
xline(pkrange)
xline(t(pks),'--r')


%% Plot 
% Double check axes orientation
zidx = 3;
figure(1)
imagesc(raster.relxs,raster.relys,squeeze(Vpk(:,:,zidx)))
colorbar
xlabel('x (mm)')
ylabel('y (mm)')
title('Volts')

xidx = 3;
figure(2)
imagesc(raster.relxs,raster.relys,squeeze(Vpk(xidx,:,:)))
colorbar
xlabel('y (mm)')
ylabel('z (mm)')
title('Volts')

yidx = 3;
figure(3)
imagesc(raster.relxs,raster.relys,squeeze(Vpk(:,yidx,:)))
colorbar
xlabel('x (mm)')
ylabel('z (mm)')
title('Volts')

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