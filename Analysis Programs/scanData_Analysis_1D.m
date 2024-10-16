%% Peak Voltage Analysis
% For reference: scanData = [x,y,z,samples,channel];
clear all

%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'DIYMk1Test5';
path = strcat(folder_path,file_name,'.mat');
load(path)
disp('Data Timestamp:')
disp(scpSettings.timestamp)
% check size of data array - note whether one or two channel
disp('Data Size:')
disp(size(scanData))

%% Check Waveform and Extract Peak Voltages
x_index = 1;
y_index = 1;
z_index = 300;

pkrange = [1,50]; % us - time range to look for peak 
pkrangeidx = pkrange*scpSettings.SampleFrequency/1e6; % corresponding array index

% remove bias
scanData_noBias = scanData(:,:,:,:,1) - mean(scanData(:,:,:,:,1),4);
Vpk = squeeze(max(scanData_noBias(:,:,:,pkrangeidx(1):pkrangeidx(2),1),[],4)); % max voltage at [x,y,z]

pks = find(scanData_noBias(x_index,y_index,z_index,:,1) == Vpk(z_index));

wvfmData = squeeze(scanData_noBias(x_index,y_index,z_index,:,1));

t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

figure(1)
plot(t,wvfmData)
xlim([0,200])
hold on
x = raster.xs(x_index);
y = raster.ys(y_index);
z = raster.zs(z_index);
title(strcat('Waveform at [x,y,z]=[',string(x),',',string(y),',',string(z),'] mm'))
xlabel('Time [us]');
ylabel('Amplitude [V]');
hold off
xline(pkrange)
xline(t(pks),'--r')



%% 

plot(raster.zs,Vpk)
xlabel('Z Position [mm]')
ylabel('Peak Voltage [V]')
x = raster.xs(x_index);
y = raster.ys(y_index);
z = raster.zs(z_index);
title(strcat('Scan Home: [x,y,z]=[',string(raster.home(1)),',',string(raster.home(2)),',',string(raster.home(3)),'] mm'))

%% To MPa
mVperMPa = 170.12; % CHECK
MPa = Vpk*1e3/mVperMPa; 
