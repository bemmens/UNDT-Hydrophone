%% Peak Voltage Analysis
% For reference: scanData.AB = [A,B,samples,channel,nrepeats];
clear all

analysisVersion = 3;

%% Load Data
folder_path = '/Users/gv19838/Library/CloudStorage/OneDrive-UniversityofBristol/PhD/Hydrophone/UNDT-Hydrophone/DataOut/';
file_name = 'DIY_Acrylic_Day4_5';
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

disp(size(scanData_noBias.YZ))

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
%% With nRpeats
Vrms.XY = mean(squeeze(rms(scanData_noBias.XY,3)),3); % rms voltage at [x,y,z0]
Vrms.YZ = mean(squeeze(rms(scanData_noBias.YZ,3)),3); % rms voltage at [x,y,z0]

% Vrms.XY = squeeze(rms(scanData_bpf.XY,3)); % rms voltage at [x,y,z0]
% Vrms.YZ = squeeze(rms(scanData_bpf.YZ,3)); % rms voltage at [x,y,z0]
% Vrms.XZ = squeeze(rms(scanData_bpf.XZ,3)); % rms voltage at [x,y,z0]



%% To MPa
mVperMPa = 151.91; % CHECK
MPa.XY = Vrms.XY*1e3/mVperMPa; 
MPa.YZ = Vrms.YZ*1e3/mVperMPa; 

%% Calculate Phase at 1MHz

% Target frequency
targetFreq = 1.06e6; % 1 MHz

% Get dimensions
[nX, nY, nSamples, nRepeats_XY] = size(scanData_noBias.XY);
[nY_yz, nZ, ~, nRepeats_YZ] = size(scanData_noBias.YZ);

% Frequency vector for FFT
freq = scpSettings.SampleFrequency * (0:floor(nSamples/2)) / nSamples;

% Find index closest to target frequency
[~, freqIdx] = min(abs(freq - targetFreq));
disp(['Calculating phase at nearest FFT bin: ', num2str(freq(freqIdx)/1e6), ' MHz']);

% Initialize phase matrices
Phase.XY = zeros(nX, nY);
Phase.YZ = zeros(nY_yz, nZ);

% Calculate phase for XY plane
for xi = 1:nX
    for yi = 1:nY
        % Average waveform across repeats
        waveform_avg = mean(squeeze(scanData_noBias.XY(xi, yi, :, :)), 2);
        % Compute FFT
        fft_result = fft(waveform_avg);
        % Extract phase at target frequency and convert to degrees
        Phase.XY(xi, yi) = rad2deg(angle(fft_result(freqIdx)));

    end
end

% Calculate phase for YZ plane
for yi = 1:nY_yz
    for zi = 1:nZ
        % Average waveform across repeats
        waveform_avg = mean(squeeze(scanData_noBias.YZ(yi, zi, :, :)), 2);
        % Compute FFT
        fft_result = fft(waveform_avg);
        % Extract phase at target frequency and convert to degrees
        Phase.YZ(yi, zi) = rad2deg(angle(fft_result(freqIdx)));

    end
end

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
figure(1)
plot(t, wvfmData_max)
hold on
% Get location information
if maxValueXY >= maxValueYZ
    max_location = ['x=', num2str(raster.xs(max_x_idx)), ', y=', num2str(raster.ys(max_y_idx)), ', z=', num2str(raster.zPlane)];
else
    max_location = ['x=', num2str(raster.home(1)), ', y=', num2str(raster.ys(max_y_idx)), ', z=', num2str(raster.zs(max_z_idx))];
end
title(['Maximum Amplitude Waveform at [', max_location, '] mm'])
xlabel('Time [us]');
ylabel('Amplitude [MPa]');
hold off

% FFT of the check waveform
N = length(wvfmData_raw1);
Y = fft(wvfmData_raw1);
P2 = abs(Y / N);
P1 = P2(1:floor(N/2)+1);
if N > 1
    P1(2:end-1) = 2 * P1(2:end-1);
end
f = scpSettings.SampleFrequency * (0:floor(N/2)) / N / 1e6; % MHz
figure(5)
plot(f, P1)
xlim([0, max(f)])
xlabel('Frequency (MHz)');
ylabel('Amplitude (MPa)');
title(['FFT of Check Waveform at x=', num2str(raster.xs(x_index)), ' mm, y=', num2str(raster.ys(y_index)), ' mm']);
grid on

% Extract the 1 MHz Fourier component from the check waveform
freqRaw = scpSettings.SampleFrequency * (0:floor(N/2)) / N;
[~, idx1MHz] = min(abs(freqRaw - targetFreq));
component1MHz = Y(idx1MHz);
phase1MHz = angle(component1MHz);
amplitude1MHz = 2 * abs(component1MHz) / N;
reconstruct1MHz = amplitude1MHz * cos(2*pi*targetFreq*(t*1e-6) + phase1MHz);

disp(['1 MHz component amplitude: ', num2str(amplitude1MHz), ' MPa, phase: ', num2str(rad2deg(phase1MHz)), ' degrees']);

figure(6)
wvfmData_raw1_norm = wvfmData_raw1 / max(abs(wvfmData_raw1));
reconstruct1MHz_norm = reconstruct1MHz *0.8/ max(abs(reconstruct1MHz));
plot(t, wvfmData_raw1_norm, 'b')
hold on
plot(t, reconstruct1MHz_norm, 'r')
hold off
xlabel('Time [us]');
ylabel('Normalized amplitude');
legend('Check waveform', '1 MHz component');
title(['Check waveform and extracted 1 MHz component at x=', num2str(raster.xs(x_index)), ' mm, y=', num2str(raster.ys(y_index)), ' mm']);
grid on


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

%% Plot 2D Phase at 1MHz
figure(4)
subplot(1, 2, 1); % Create a subplot for XY Plane
imagesc(raster.xs, raster.ys, Phase.XY')
cb = colorbar;
% cb.Label.String = 'Phase (degrees)';
xlabel('X (mm)');
ylabel('Y (mm)');
title(strcat('Phase at 1MHz in XY Plane at z =',{' '}, string(round(raster.zPlane, 1)), ' mm'));
axis xy; % Correct the axis direction
axis image; % Set aspect ratio suitable for images

subplot(1, 2, 2); % Create a subplot for YZ Plane
imagesc(raster.ys, raster.zs, Phase.YZ')
cb = colorbar;
% cb.Label.String = 'Phase (degrees)';
xlabel('Y (mm)');
ylabel('Z (mm)');
title(strcat('Phase at 1MHz in YZ Plane at x =',{' '}, string(round(raster.xs(x_index), 1)), ' mm'));
axis xy; % Correct the axis direction
axis image; % Set aspect ratio suitable for images
set(gca, 'YDir', 'reverse'); % Reverse the z-axis