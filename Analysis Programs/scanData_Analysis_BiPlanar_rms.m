%% Peak Voltage Analysis
% For reference: scanData.AB = [A,B,samples,channel,nrepeats];
clear all

analysisVersion = 3;

%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'NearSurface_DIYMk1_5';
path = strcat(folder_path,file_name,'.mat');
load(path)
disp('Data Timestamp:')
disp(scpSettings.timestamp)
if analysisVersion == 0
    disp('Skipping compatability check. Data dimensions may be incompatible with requested plots.')
elseif analysisVersion ~= scpSettings.scanVersion
    warning("Scan and analysis programs have mismatched versions.")
end


%%
% Useful Constants
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

%% Post-Process Data

pkrange = [1,100]; % us - time range to look for peak 
pkrangeidx = pkrange*scpSettings.SampleFrequency/1e6; % corresponding array index

% remove trigger (2nd) channel
scanData_noTrigger.XY = squeeze(scanData.XY(:,:,:,1,:));
scanData_noTrigger.YZ = squeeze(scanData.YZ(:,:,:,1,:));

% remove bias
scanData_noBias.XY = scanData_noTrigger.XY - mean(scanData_noTrigger.XY,3);
scanData_noBias.YZ = scanData_noTrigger.YZ - mean(scanData_noTrigger.YZ,3);


%% Bandpass filter 

x_index = 4;
y_index = 4;
z_index = 10;

Fs = scpSettings.SampleFrequency; % Sampling Frequency
F0 = 1*1e6; % Centre
width = 0.25*1e6;
Fpass1 = F0-width; % First Passband Frequency
Fpass2 = F0+width; % Second Passband Frequency

% Apply the bandpass filter
figure(1)
bandpass(squeeze(scanData_noBias.XY(x_index,y_index,:))', [Fpass1 Fpass2], Fs)
[~,bpfilter] = bandpass(squeeze(scanData_noBias.XY(x_index,y_index,:))', [Fpass1 Fpass2], Fs);

scanData_bpf.XY = filter( bpfilter.Coefficients, 1, scanData_noBias.XY, [], 3);
scanData_bpf.YZ = filter( bpfilter.Coefficients, 1, scanData_noBias.YZ, [], 3);

figure(100)
plot(t,squeeze(scanData_bpf.XY(x_index,y_index,:,1)))
hold on
%plot(t,squeeze(scanData_noBias.XY(10,10,:)))
hold off
%xlim([100,150])
%% With nRpeats
Vrms.XY = mean(squeeze(rms(scanData_noBias.XY,3)),3); % rms voltage at [x,y,z0]
Vrms.YZ = mean(squeeze(rms(scanData_noBias.YZ,3)),3); % rms voltage at [x,y,z0]

% Vrms.XY = squeeze(rms(scanData_bpf.XY,3)); % rms voltage at [x,y,z0]
% Vrms.YZ = squeeze(rms(scanData_bpf.YZ,3)); % rms voltage at [x,y,z0]
% Vrms.XZ = squeeze(rms(scanData_bpf.XZ,3)); % rms voltage at [x,y,z0]

%% To MPa
mVperMPa = 170.49; % CHECK
MPa.XY = Vrms.XY*1e3/mVperMPa; 
MPa.YZ = Vrms.YZ*1e3/mVperMPa; 

%% Check Waveform 

wvfmData_raw1 = squeeze(scanData_noBias.XY(x_index,y_index,:,1))*1e3/mVperMPa;

%% Plots
figure(1)
plot(t,wvfmData_raw1)
hold on
%plot(t,wvfmData_raw2)
%plot(t,wvfmData_raw3)
%xlim([0,200])
x = raster.xs(x_index);
y = raster.ys(y_index);
z = raster.zs(z_index);
title(strcat('Raw Data at [x,y,z] = [',string(x),', ',string(y),', ',string(z),'] mm'))
xlabel('Time [us]');
ylabel('Amplitude [MPa]');
hold off
%legend('Raw Waveform','Mean waveform','pkrange min','pkrange max','Vrms')
% xlim([0,10])

%% Coords relative to plot

relX = raster.xs - raster.home(1);
relY = raster.ys - raster.home(2);
relZ = flip(raster.zs - raster.home(3));

%% Plot 3D orthogonal views - CoPilot
figure(2)
hold on

% XY plane
[X1, Y1] = meshgrid(raster.xs, raster.ys);
Z1 = ones(size(X1)) * raster.home(3);
surf(X1, Y1, Z1, MPa.XY', 'EdgeColor', 'none')

% YZ plane
[Y2, Z2] = meshgrid(raster.ys, raster.zs);
X2 = ones(size(Y2)) * raster.home(1);
surf(X2, Y2, Z2, MPa.YZ', 'EdgeColor', 'none')

cb = colorbar;
cb.Label.String = 'Pressure (MPa)';
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
%xlim([15,30])
%ylim([15,33])
title('Bi-Planar Scan of Acoustic Field')
subtitle('Stage Coordinates')
view(3)
axis vis3d
hold off
rotate3d

%% Plot 3D orthogonal views - Relative to scan centre
figure(3)
hold on

% XY plane
[X1, Y1] = meshgrid(relX, relY);
Z1 = zeros(size(X1));
surf(X1, Y1, Z1, MPa.XY', 'EdgeColor', 'none')

% YZ plane
[Y2, Z2] = meshgrid(relY, relZ);
X2 = zeros(size(Y2));
surf(X2, Y2, Z2, MPa.YZ', 'EdgeColor', 'none')

cb = colorbar;
cb.Label.String = 'Pressure (MPa)';
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
%xlim([15,30])
%ylim([15,33])
title('Bi-Planar Scan of Acoustic Field')
subtitle('Scan Coordinates')
view(3)
axis vis3d
hold off
rotate3d






