%% Peak Voltage Analysis
% For reference: scanData.AB = [A,B,samples,channel,nrepeats];
clear all

analysisVersion = 3;

%% Load Data
folder_path = '/Users/gv19838/Library/CloudStorage/OneDrive-UniversityofBristol/PhD/Hydrophone/UNDT-Hydrophone/DataOut/';
file_name = 'DIYMk1_Char2';
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

disp(size(scanData_noBias.XY))

% %% Bandpass filter 
% 
x_index = 20;
y_index = 40;
z_index = 5;


%% With nRpeats
Vrms.XY = mean(squeeze(rms(scanData_noBias.XY,3)),3); % rms voltage at [x,y,z0]
Vrms.YZ = mean(squeeze(rms(scanData_noBias.YZ,3)),3); % rms voltage at [x,y,z0]

%% To MPa
mVperMPa = 170.49; % CHECK
MPa.XY = Vrms.XY*1e3/mVperMPa; 
MPa.YZ = Vrms.YZ*1e3/mVperMPa; 

%% Check Waveform 

wvfmData_raw1 = squeeze(scanData_noBias.XY(x_index,y_index,:,1))*1e3/mVperMPa;

% Find maximum voltage across both XY and YZ datasets
[maxValueXY, maxIndexXY] = max(abs(scanData_noBias.XY(:)));
[maxValueYZ, maxIndexYZ] = max(abs(scanData_noBias.YZ(:)));

% Determine which dataset has the higher maximum
if maxValueXY >= maxValueYZ
    % Get the indices for the maximum value in XY data
    [max_x_idx, max_y_idx, max_time_idx, max_repeat] = ind2sub(size(scanData_noBias.XY), maxIndexXY);
    
    % Display information about the maximum value
    disp('Maximum voltage found in XY scan data');
    disp(['Location: x=', num2str(raster.xs(max_x_idx)), ' mm, y=', num2str(raster.ys(max_y_idx)), ' mm, z=', num2str(raster.zPlane), ' mm']);
    disp(['Maximum absolute value: ', num2str(maxValueXY), ' V at time index ', num2str(max_time_idx)]);
    
    % Extract the complete waveform at the location of maximum voltage
    wvfmData_max = squeeze(scanData_noBias.XY(max_x_idx, max_y_idx, :, max_repeat))*1e3/mVperMPa;
else
    % Get the indices for the maximum value in YZ data
    [max_y_idx, max_z_idx, max_time_idx, max_repeat] = ind2sub(size(scanData_noBias.YZ), maxIndexYZ);
    
    % Display information about the maximum value
    disp('Maximum voltage found in YZ scan data');
    disp(['Location: x=', num2str(raster.home(1)), ' mm, y=', num2str(raster.ys(max_y_idx)), ' mm, z=', num2str(raster.zs(max_z_idx)), ' mm']);
    disp(['Maximum absolute value: ', num2str(maxValueYZ), ' V at time index ', num2str(max_time_idx)]);
    
    % Extract the complete waveform at the location of maximum voltage
    wvfmData_max = squeeze(scanData_noBias.YZ(max_y_idx, max_z_idx, :, max_repeat))*1e3/mVperMPa;
end
%% Plots

figure()
plot(t, wvfmData_max)
hold on
% Get location information
if maxValueXY >= maxValueYZ
    x = raster.xs(max_x_idx);
    y = raster.ys(max_y_idx);
    max_location = ['x=', num2str(raster.xs(max_x_idx)), ', y=', num2str(raster.ys(max_y_idx)), ', z=', num2str(raster.zPlane)];
else
    y = raster.ys(max_y_idx);
    z = raster.zs(max_z_idx);
    max_location = ['x=', num2str(raster.home(1)), ', y=', num2str(raster.ys(max_y_idx)), ', z=', num2str(raster.zs(max_z_idx))];
end
% Plot horizontal line at RMS of the pressure line
yline(rms(wvfmData_max(pkrangeidx(1):pkrangeidx(2))));
% Add horizontal lines for the RMS of a sine wave with the same amplitude as the peak
Amp = max(wvfmData_max);
sine_rms = Amp / sqrt(2); % RMS of a sine wave with amplitude A
yline(sine_rms,'--');
legend('Waveform', 'RMS', 'Location', 'best');
% Calculate RMS and max values for the waveform
wvfm_rms = rms(wvfmData_max(pkrangeidx(1):pkrangeidx(2)));
wvfm_max = max(abs(wvfmData_max));
% Add RMS and max values to the subtitle
subtitle(sprintf('RMS = %.3f MPa, Max = %.3f MPa', wvfm_rms, wvfm_max));
title(['Maximum Amplitude Waveform at [', max_location, '] mm'])
xlabel('Time [us]');
ylabel('Amplitude [MPa]');
legend('Waveform', 'RMS', 'RMS of Sine', 'Location', 'Best');
hold off

%%

% Plot pressure along the YZ line that contains the overall maximum
% Determine y index of overall maximum if not already available

[maxValueYZ, maxIndexYZ] = max(abs(scanData_noBias.YZ(:)));
[max_y_idx, ~, ~, ~] = ind2sub(size(scanData_noBias.YZ), maxIndexYZ);

y_line_idx = max_y_idx;
pressure_line = MPa.YZ(y_line_idx, :);   % MPa.YZ rows = y, cols = z
zvals = raster.zs;

figure()
plot(zvals, pressure_line)
hold on
[~, zmax_idx] = max(abs(pressure_line));
plot(zvals(zmax_idx), pressure_line(zmax_idx), 'r*', 'MarkerSize',10)

hold off
grid on
xlabel('Z (mm)')
ylabel('RMS Pressure (MPa)')
title(sprintf('Pressure vs Z at y = %.2f mm (index %d)', raster.ys(y_line_idx), y_line_idx))
legend('Pressure', 'Max (abs)', 'Location', 'Best')

% Fit a sine to the pressure_line vs zvals using dominant spatial frequency (FFT)
z = zvals(:);
p = pressure_line(:);

% Estimate dominant spatial frequency from FFT
dz = mean(diff(z));
% Estimate dominant spatial frequency and wavelength
P = fft(p - mean(p));
[~, maxIdx] = max(abs(P(2:floor(end/2))));
f_dom = maxIdx / (length(P) * dz); % cycles per mm
L = 1 / f_dom; % wavelength in mm

% Fit sine wave: A*sin(2*pi*z/L + phi) + C
coeff = [sin(2*pi*z / L), cos(2*pi*z / L), ones(size(z))] \ p;
A = norm(coeff(1:2)); % amplitude
phi = atan2(coeff(2), coeff(1)); % phase
C0 = coeff(3); % offset

% Generate a smoother fit using more points
z = linspace(min(z), max(z), 10 * length(z)); % Increase resolution
p_fit = A * sin(2*pi*z / L + phi) + C0;

% Plot fitted curve on current figure
figure()
hold on
plot(z, p_fit)
hold off

% Add subtitle with wavelength
subtitle(sprintf('Fitted wavelength = %.3f mm (Amplitude = %.3f MPa)', L, A))

%% Coords relative to plot

relX = raster.xs - raster.home(1);
relY = raster.ys - raster.home(2);
relZ = flip(raster.zs - raster.home(3));

%% Plot 3D orthogonal views - CoPilot
figure()
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

% %% Plot 3D orthogonal views - Relative to scan centre
% figure()
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
figure()
subplot(1, 2, 1); % Create a subplot for XY Plane
imagesc(raster.xs, raster.ys, MPa.XY')
cb = colorbar;
cb.Label.String = 'Pressure (MPa)';
xlabel('X (mm)');
ylabel('Y (mm)');
title(strcat('Pressure in XY Plane at z =',{' '}, string(round(raster.zPlane, 1)), ' mm'));
axis xy; % Correct the axis direction
axis image; % Set aspect ratio suitable for images

subplot(1, 2, 2); % Create a subplot for YZ Plane
imagesc(raster.ys, raster.zs, MPa.YZ') % Flip the z-axis
cb = colorbar;
cb.Label.String = 'Pressure (MPa)';
xlabel('Y (mm)');
ylabel('Z (mm)');
title(strcat('Pressure in YZ Plane at x =',{' '}, string(x), ' mm'));
axis xy; % Correct the axis direction
axis image; % Set aspect ratio suitable for images
set(gca, 'YDir', 'reverse'); % Reverse the z-axis