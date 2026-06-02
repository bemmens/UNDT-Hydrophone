%% Peak Voltage Analysis
% For reference: scanData.AB = [A,B,samples,channel,nrepeats];
clear all

analysisVersion = 3;

%% Load Data using matfile for memory efficiency
folder_path = '/Users/gv19838/Library/CloudStorage/OneDrive-UniversityofBristol/PhD/Hydrophone/UNDT-Hydrophone/DataOut/';
file_name = 'DIY_Acrylic_5';
path = strcat(folder_path,file_name,'.mat');

% Open matfile without loading entire file into memory
mf = matfile(path, 'Writable', false);

% Load metadata (small variables - loaded entirely)
scpSettings = mf.scpSettings;
raster = mf.raster;

% Load XY and YZ sequentially to stay under 16GB memory limit
% First load XY (approximately 15GB for single channel)
disp('Loading scanData.XY (channel 1 only)...')
scanData_full = mf.scanData;
scanData.XY = squeeze(scanData_full.XY(:,:,:,1,:));
disp(['Loaded XY: ', num2str(numel(scanData.XY)), ' elements'])

% Clear and reload for YZ to avoid holding both in memory simultaneously
clear scanData_full
disp('Clearing memory...')
pause(0.1)

% Now load YZ
disp('Loading scanData.YZ (channel 1 only)...')
scanData_full = mf.scanData;
scanData.YZ = squeeze(scanData_full.YZ(:,:,:,1,:));
disp(['Loaded YZ: ', num2str(numel(scanData.YZ)), ' elements'])
clear scanData_full mf

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

% Pre-allocate RMS output arrays
Vrms.XY = zeros(length(raster.xs), length(raster.ys));
Vrms.YZ = zeros(length(raster.ys), length(raster.zs));

% Track maximum value and location across all data
maxValueXY = 0;
maxIndexXY = [];
max_x_idx_XY = 0;
max_y_idx_XY = 0;
max_time_idx_XY = 0;
max_repeat_XY = 0;

maxValueYZ = 0;
maxIndexYZ = [];
max_y_idx_YZ = 0;
max_z_idx_YZ = 0;
max_time_idx_YZ = 0;
max_repeat_YZ = 0;

% Open matfile for chunked loading
mf = matfile(path, 'Writable', false);

% Process XY data: load and process one X-slice at a time
disp('Processing XY data in chunks...')
for x_idx = 1:length(raster.xs)
    % Load single X-slice: squeeze removes singleton dimensions
    xy_slice = squeeze(mf.scanData.XY(x_idx, :, :, 1, :));  % [y, samples, repeats]
    
    % Remove bias (subtract mean across samples for this slice)
    xy_slice_noBias = xy_slice - mean(xy_slice, 2);
    
    % Compute RMS across samples for this X-slice
    Vrms.XY(x_idx, :) = mean(squeeze(rms(xy_slice_noBias, 2)), 2);  % [y]
    
    % Find local maximum in this slice
    [local_max, local_max_idx] = max(abs(xy_slice_noBias(:)));
    if local_max > maxValueXY
        maxValueXY = local_max;
        [max_y_idx_XY, max_time_idx_XY, max_repeat_XY] = ind2sub(size(xy_slice_noBias), local_max_idx);
        max_x_idx_XY = x_idx;
    end
    
    if mod(x_idx, 10) == 0
        disp(['  Processed X slice ', num2str(x_idx), ' of ', num2str(length(raster.xs))])
    end
    clear xy_slice xy_slice_noBias
end

% Process YZ data: load and process one Y-slice at a time
disp('Processing YZ data in chunks...')
for y_idx = 1:length(raster.ys)
    % Load single Y-slice
    yz_slice = squeeze(mf.scanData.YZ(y_idx, :, :, 1, :));  % [z, samples, repeats]
    
    % Remove bias
    yz_slice_noBias = yz_slice - mean(yz_slice, 2);
    
    % Compute RMS across samples for this Y-slice
    Vrms.YZ(y_idx, :) = mean(squeeze(rms(yz_slice_noBias, 2)), 2);  % [z]
    
    % Find local maximum in this slice
    [local_max, local_max_idx] = max(abs(yz_slice_noBias(:)));
    if local_max > maxValueYZ
        maxValueYZ = local_max;
        [max_z_idx_YZ, max_time_idx_YZ, max_repeat_YZ] = ind2sub(size(yz_slice_noBias), local_max_idx);
        max_y_idx_YZ = y_idx;
    end
    
    if mod(y_idx, 10) == 0
        disp(['  Processed Y slice ', num2str(y_idx), ' of ', num2str(length(raster.ys))])
    end
    clear yz_slice yz_slice_noBias
end

disp('Chunk processing complete.')
disp(size(Vrms.XY))

% %% Bandpass filter 
% 
x_index = 13;
y_index = 13;
z_index = 5;


% Fs = scpSettings.SampleFrequency; % Sampling Frequency
% F0 = 1*1e6; % Centre
% width = 0.25*1e6;
% Fpass1 = 5e5; % First Passband Frequency
% Fpass2 = F0+width; % Second Passband Frequency
% 
% % Apply the bandpass filter
% figure(1)
% bandpass(squeeze(scanData_noBias.XY(x_index,y_index,:))', [Fpass1 Fpass2], Fs)
% [~,bpfilter] = bandpass(squeeze(scanData_noBias.XY(x_index,y_index,:))', [Fpass1 Fpass2], Fs);
% 
% scanData_bpf.XY = filter( bpfilter.Coefficients, 1, scanData_noBias.XY, [], 3);
% scanData_bpf.YZ = filter( bpfilter.Coefficients, 1, scanData_noBias.YZ, [], 3);
% 
% figure(100)
% plot(t,squeeze(scanData_bpf.XY(x_index,y_index,:,1)))
% hold on
% %plot(t,squeeze(scanData_noBias.XY(10,10,:)))
% hold off
% %xlim([100,150])
%% To MPa (Vrms already computed in chunk processing loop above)
mVperMPa = 170.49; % CHECK
MPa.XY = Vrms.XY*1e3/mVperMPa; 
MPa.YZ = Vrms.YZ*1e3/mVperMPa; 

%% Check Waveform 
% Reload only the slice containing the maximum value
if maxValueXY >= maxValueYZ
    disp('Maximum voltage found in XY scan data');
    disp(['Location: x=', num2str(raster.xs(max_x_idx_XY)), ' mm, y=', num2str(raster.ys(max_y_idx_XY)), ' mm, z=', num2str(raster.zPlane), ' mm']);
    disp(['Maximum absolute value: ', num2str(maxValueXY), ' V at time index ', num2str(max_time_idx_XY)]);
    
    % Reload only the X-slice containing the max
    xy_max_slice = squeeze(mf.scanData.XY(max_x_idx_XY, :, :, 1, :));  % [y, samples, repeats]
    xy_max_slice_noBias = xy_max_slice - mean(xy_max_slice, 2);
    wvfmData_max = squeeze(xy_max_slice_noBias(max_y_idx_XY, :, max_repeat_XY))*1e3/mVperMPa;
    max_location = ['x=', num2str(raster.xs(max_x_idx_XY)), ', y=', num2str(raster.ys(max_y_idx_XY)), ', z=', num2str(raster.zPlane)];
else
    disp('Maximum voltage found in YZ scan data');
    disp(['Location: x=', num2str(raster.home(1)), ' mm, y=', num2str(raster.ys(max_y_idx_YZ)), ' mm, z=', num2str(raster.zs(max_z_idx_YZ)), ' mm']);
    disp(['Maximum absolute value: ', num2str(maxValueYZ), ' V at time index ', num2str(max_time_idx_YZ)]);
    
    % Reload only the Y-slice containing the max
    yz_max_slice = squeeze(mf.scanData.YZ(max_y_idx_YZ, :, :, 1, :));  % [z, samples, repeats]
    yz_max_slice_noBias = yz_max_slice - mean(yz_max_slice, 2);
    wvfmData_max = squeeze(yz_max_slice_noBias(max_z_idx_YZ, :, max_repeat_YZ))*1e3/mVperMPa;
    max_location = ['x=', num2str(raster.home(1)), ', y=', num2str(raster.ys(max_y_idx_YZ)), ', z=', num2str(raster.zs(max_z_idx_YZ))];
end

clear mf xy_max_slice xy_max_slice_noBias yz_max_slice yz_max_slice_noBias

%% Plots
figure(1)
plot(t, wvfmData_max)
hold on
title(['Maximum Amplitude Waveform at [', max_location, '] mm'])
xlabel('Time [us]');
ylabel('Amplitude [MPa]');
hold off

%% Coords relative to plot

relX = raster.xs - raster.home(1);
relY = raster.ys - raster.home(2);
relZ = flip(raster.zs - raster.home(3));

%% Plot 3D orthogonal views - CoPilot
figure(2)
hold on

% XY plane
[X1, Y1] = meshgrid(raster.xs, raster.ys);
Z1 = ones(size(X1)) * raster.zPlane;
surf(X1, Y1, Z1, MPa.XY', 'EdgeColor', 'none')

% YZ plane
[Y2, Z2] = meshgrid(raster.ys, raster.zs);
X2 = ones(size(Y2)) * raster.home(1);
surf(X2, Y2, Z2, MPa.YZ', 'EdgeColor', 'none')

set(gca, 'ZDir', 'reverse')

cb = colorbar;
cb.Label.String = 'RMS Pressure (MPa)';
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

% %% Plot 3D orthogonal views - Relative to scan centre
% figure(3)
% hold on
% 
% % XY plane
% [X1, Y1] = meshgrid(relX, relY);
% Z1 = zeros(size(X1));
% surf(X1, Y1, Z1, MPa.XY', 'EdgeColor', 'none')
% 
% % YZ plane
% [Y2, Z2] = meshgrid(relY, relZ);
% X2 = zeros(size(Y2));
% surf(X2, Y2, Z2, MPa.YZ', 'EdgeColor', 'none')
% 
% cb = colorbar;
% cb.Label.String = 'Pressure (MPa)';
% xlabel('x (mm)')
% ylabel('y (mm)')
% zlabel('z (mm)')
% %xlim([15,30])
% %ylim([15,33])
% title('Bi-Planar Scan of Acoustic Field')
% subtitle('Scan Coordinates')
% view(3)
% axis vis3d
% hold off
% rotate3d

%% 
%% Plot 2D XY Plane
figure(3)
subplot(1, 2, 1); % Create a subplot for XY Plane
imagesc(raster.xs, raster.ys, MPa.XY')
cb = colorbar;
cb.Label.String = 'RMS Pressure (MPa)';
xlabel('X (mm)');
ylabel('Y (mm)');
title(strcat('Pressure in XY Plane at z =',{' '}, string(round(raster.zPlane, 1)), ' mm'));
axis xy; % Correct the axis direction
axis image; % Set aspect ratio suitable for images

subplot(1, 2, 2); % Create a subplot for YZ Plane
imagesc(raster.ys, raster.zs, MPa.YZ') % Flip the z-axis
cb = colorbar;
cb.Label.String = 'RMS Pressure (MPa)';
xlabel('Y (mm)');
ylabel('Z (mm)');
title(strcat('Pressure in YZ Plane at x =',{' '}, string(round(raster.xs(x_index), 1)), ' mm'));
axis xy; % Correct the axis direction
axis image; % Set aspect ratio suitable for images
set(gca, 'YDir', 'reverse'); % Reverse the z-axis