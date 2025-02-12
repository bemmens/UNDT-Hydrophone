%% Peak Voltage Analysis
% For reference: scanData = [x,y,z,samples,channel];
clear all

%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'DIYMk1Test24';
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
z_index = 1;

pkrange = [1,60]; % us - time range to look for peak 
pkrangeidx = pkrange*scpSettings.SampleFrequency/1e6; % corresponding array index

% remove bias
scanData_noBias = scanData(:,:,:,:,1) - mean(scanData(:,:,:,:,1),4);
Vpk = squeeze(max(scanData_noBias(:,:,:,pkrangeidx(1):pkrangeidx(2),1),[],4)); % max voltage at [x,y,z]

pks = find(scanData_noBias(x_index,y_index,z_index,:,1) == Vpk(x_index,y_index,z_index));

wvfmData = squeeze(scanData_noBias(x_index,y_index,z_index,:,1));

t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

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
isosurface(MPa,0);
hold on
isosurface(MPa,0.5)
isosurface(MPa,0.75)
hold off
legend('0MPa','0.5MPa','0.75MPa')

%% Isosurface
isosurface(smooth3(MPa,"box",3))

%% Scatter
scatter3()

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
