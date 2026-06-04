%% Peak Voltage Analysis
% For reference: scanData.AB = [A,B,samples,channel,nrepeats];
clear all

analysisVersion = 4;

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

% remove bias
scanData_noBias.XY = scanData_noTrigger.XY - mean(scanData_noTrigger.XY,3);

disp(size(scanData_noBias.XY))

% %% Bandpass filter 
% 
x_index = 10;
y_index = 10;


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

%% To MPa
mVperMPa = 232.40; % CHECK
MPa.XY = Vrms.XY*1e3/mVperMPa; 

%% Check Waveform 

wvfmData_raw1 = squeeze(scanData_noBias.XY(x_index,y_index,:,1))*1e3/mVperMPa;
wvfmData_mean = squeeze(mean(scanData_noBias.XY(x_index,y_index,:,:),4))*1e3/mVperMPa;
disp(size(wvfmData_mean))

%% Plots
figure(1)
plot(t,wvfmData_raw1)
hold on
%plot(t,wvfmData_raw2)
%plot(t,wvfmData_raw3)
plot(t,wvfmData_mean)

%xlim([0,200])
x = raster.xs(x_index);
y = raster.ys(y_index);
% title(strcat('Raw Data at [x,y,z] = [',string(x),', ',string(y)] mm'))
xlabel('Time [us]');
ylabel('Amplitude [MPa]');
hold off
legend('Raw Waveform','Mean waveform')
% xlim([0,10])

%% Bandpass filter 
Fs = scpSettings.SampleFrequency; % Sampling Frequency
F0 = 4.5*1e6; % Centre
width = 0.3*1e6;
Fpass1 = F0-width; % First Passband Frequency
Fpass2 = F0+width; % Second Passband Frequency
figure(2)
bandpass(wvfmData_raw1, [Fpass1 Fpass2], Fs)

%% Coords relative to plot

relX = raster.xs - raster.home(1);
relY = raster.ys - raster.home(2);


%% Plot 2D XY Plane
figure(3)
imagesc(raster.xs, raster.ys, MPa.XY')
cb = colorbar;
cb.Label.String = 'Pressure (MPa)';
xlabel('X (mm)');
ylabel('Y (mm)');
title(strcat('Pressure in XY Plane at z =',{' '}, string(round(raster.home(3), 1)), ' mm'));
axis xy; % Correct the axis direction
axis image; % Set aspect ratio suitable for images
