%% Peak Voltage Analysis
% For reference: scanData = [x,y,z,samples,channel,nrepeats];
clear all

analysisVersion = 1
%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'repeats_test';
path = strcat(folder_path,file_name,'.mat');
load(path)
disp('Data Timestamp:')
disp(scpSettings.timestamp)
if analysisVersion == 0
    disp('Skipping compatability check. Data dimensions may be incompatible with requested plots.')
elseif analysisVersion ~= scpSettings.scanVersion
    warning("Scan and analysis programs have mismatched versions.")
end
% check size of data array - note whether one or two channel
disp('Data Size:')
disp(size(scanData))

%%
% Useful Constants
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

%% Post-Process Data

pkrange = [1,60]; % us - time range to look for peak 
pkrangeidx = pkrange*scpSettings.SampleFrequency/1e6; % corresponding array index

% remove trigger (2nd) channel
scanData_noTrigger = squeeze(scanData(:,:,:,:,1,:));

% remove bias
scanData_noBias = scanData(:,:,:,:,1,:) - mean(scanData(:,:,:,:,1,:),4);

% take mean & std
[scanData_std,scanData_mean] = std(scanData_noBias,0,6);
%% PROGRESS MARK

Vpk = squeeze(max(scanData_noBias(:,:,:,pkrangeidx(1):pkrangeidx(2),1),[],4)); % max voltage at [x,y,z]

%% Check Waveform 
x_index = 1;
y_index = 1;
z_index = 1;

pks = find(scanData_noBias(x_index,y_index,z_index,:,1) == Vpk(x_index,y_index,z_index));
wvfmData = squeeze(scanData_noBias(x_index,y_index,z_index,:,1));

figure(1)
plot(t,wvfmData)
%xlim([0,200])
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

%% To MPa
mVperMPa = 153.23; % CHECK
MPa = Vpk*1e3/mVperMPa; 
figure(10)
imagesc(raster.xs, flip(raster.ys), rot90(MPa(:,:,z_index)))
axis image;
set(gca, 'YDir', 'normal')
a = colorbar;
a.Label.String = 'MPa';
xlabel('x (mm)')
ylabel('y (mm)')
title(strcat('z= ', string(raster.zs(z_index)), 'mm'))

% Add Slider
nZs = length(raster.zs);
slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', nZs, 'Value', z_index, 'Position', [20 20 200 20]);
addlistener(slider, 'Value', 'PostSet', @(~,~) updatePlot(MPa,round(slider.Value), raster));

%%
figure(4)
isosurface(smooth3(MPa,"box",3))
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z(mm)')
title('Pressure Isosurface')

%% functions

% Update Plot Function
function updatePlot(MPa,z_index, raster)
    figure(10)
    imagesc(raster.xs, flip(raster.ys), rot90(MPa(:,:,z_index)))
    axis image;
    set(gca, 'YDir', 'normal')
    a = colorbar;
    a.Label.String = 'MPa';
    xlabel('x (mm)')
    ylabel('y (mm)')
    title(strcat('z= ', string(raster.zs(z_index)), 'mm'))
end