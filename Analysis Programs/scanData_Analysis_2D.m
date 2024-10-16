%% Peak Voltage Analysis
% For reference: scanData = [x,y,z,samples,channel];
clear all

%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'DIYMk1Test26';
path = strcat(folder_path,file_name,'.mat');
load(path)
disp('Data Timestamp:')
disp(scpSettings.timestamp)
% check size of data array - note whether one or two channel
disp('Data Size:')
disp(size(scanData))

scanData = squeeze(scanData(:,1,:,:,1));
disp(size(scanData))

%% Check Waveform and Extract Peak Voltages
x_index = 1;
y_index = 1;

pkrange = [1,100]; % us - time range to look for peak 
pkrangeidx = pkrange*scpSettings.SampleFrequency/1e6; % corresponding array index

% remove bias
scanData_noBias = scanData - mean(scanData,3);
disp(size(scanData_noBias))


Vpk = squeeze(max(scanData_noBias(:,:,pkrangeidx(1):pkrangeidx(2)),[],3)); % max voltage at [x,y,z]
disp(size(Vpk))

pks = find(scanData_noBias(x_index,y_index,:) == Vpk(x_index,y_index));

wvfmData = squeeze(scanData_noBias(x_index,y_index,:));

t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

figure(1)
plot(t,wvfmData)
%xlim([0,200])
%hold on
x = raster.xs(x_index);
y = raster.zs(y_index);
title(strcat('Waveform at [x,z]=[',string(x),',',string(y),'] mm'))
xlabel('Time [us]');
ylabel('Amplitude [V]');
hold off
xline(pkrange)
xline(t(pks),'--r')

%% To MPa
mVperMPa = 153.23; % CHECK
MPa = Vpk*1e3/mVperMPa; 

figure(2)
imagesc(raster.xs,flip(raster.zs),rot90(MPa))
axis image;
set(gca,'YDir','normal') % rot90, flip to get stage coords to match image
a=colorbar;
a.Label.String = 'MPa';
xlabel('x (mm)')
ylabel('z (mm)')
title(strcat('y= ',string(raster.ys),'mm'))
%clim([0.1,0.8])

