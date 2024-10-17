%% Peak Voltage Analysis
% For reference: scanData.AB = [A,B,samples,channel,nrepeats];
clear all

analysisVersion = 2;

%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'tri_planar_test';
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

pkrange = [1,60]; % us - time range to look for peak 
pkrangeidx = pkrange*scpSettings.SampleFrequency/1e6; % corresponding array index

% remove trigger (2nd) channel
scanData_noTrigger.XY = squeeze(scanData.XY(:,:,:,1,:));
scanData_noTrigger.YZ = squeeze(scanData.YZ(:,:,:,1,:));
scanData_noTrigger.XZ = squeeze(scanData.XZ(:,:,:,1,:));

% remove bias
scanData_noBias.XY = scanData_noTrigger.XY - mean(scanData_noTrigger.XY,3);
scanData_noBias.YZ = scanData_noTrigger.YZ - mean(scanData_noTrigger.YZ,3);
scanData_noBias.XZ = scanData_noTrigger.XZ - mean(scanData_noTrigger.XZ,3);

% take mean & std of repeats
[scanData_std.XY,scanData_mean.XY] = std(scanData_noBias.XY,0,4);
[scanData_std.YZ,scanData_mean.YZ] = std(scanData_noBias.YZ,0,4);
[scanData_std.XZ,scanData_mean.XZ] = std(scanData_noBias.XZ,0,4);

Vpk.XY = squeeze(max(scanData_mean.XY(:,:,pkrangeidx(1):pkrangeidx(2)),[],3)); % max voltage at [x,y,z0]
Vpk.YZ = squeeze(max(scanData_mean.YZ(:,:,pkrangeidx(1):pkrangeidx(2)),[],3)); % max voltage at [x,y,z0]
Vpk.XZ = squeeze(max(scanData_mean.XZ(:,:,pkrangeidx(1):pkrangeidx(2)),[],3)); % max voltage at [x,y,z0]

%% Check Waveform 
x_index = 1;
y_index = 1;

pks = find(scanData_mean.XY(x_index,y_index,:,1) == Vpk.XY(x_index,y_index));
wvfmData = squeeze(scanData_noBias.XY(x_index,y_index,:,1));

figure(1)
plot(t,wvfmData)
%xlim([0,200])
hold on
x = raster.xs(x_index);
y = raster.ys(y_index);
z = raster.home(3);
title(strcat('Raw Data at [x,y,z] = [',string(x),', ',string(y),', ',string(z),'] mm'))
xlabel('Time [us]');
ylabel('Amplitude [V]');
hold off
xline(pkrange)
xline(t(pks),'--r')
legend('Waveform','pkrange min','pkrange max','Vpk')

%% To MPa
mVperMPa = 153.23; % CHECK
MPa.XY = Vpk.XY*1e3/mVperMPa; 
MPa.YZ = Vpk.YZ*1e3/mVperMPa; 
MPa.XZ = Vpk.XZ*1e3/mVperMPa; 

%% Plot 3D orthogonal views - CoPilot
figure(2)
hold on

% XY plane
[X, Y] = meshgrid(raster.xs, raster.ys);
Z = ones(size(X)) * raster.home(3);
surf(X, Y, Z, MPa.XY, 'EdgeColor', 'none')

% YZ plane
[Y, Z] = meshgrid(raster.ys, raster.zs);
X = ones(size(Y)) * raster.home(1);
surf(X, Y, Z, MPa.YZ, 'EdgeColor', 'none')

% XZ plane
[X, Z] = meshgrid(raster.xs, raster.zs);
Y = ones(size(X)) * raster.home(2);
surf(X, Y, Z, MPa.XZ, 'EdgeColor', 'none')

cb = colorbar;
cb.Label.String = 'Pressure (MPa)';
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title('Tri-Planar Scan of Acoustic Field')
view(3)
axis vis3d
hold off


