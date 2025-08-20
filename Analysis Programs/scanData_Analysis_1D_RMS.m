%% Peak Voltage Analysis
% For reference: scanData = [x,y,z,samples,channel];
clear all

%% Load Data
folder_path = 'C:\Users\Public\Documents\GitHub\UNDT-Hydrophone\DataOut\';
file_name = '1DTest';
path = strcat(folder_path,file_name,'.mat');
load(path)
disp('Data Timestamp:')
disp(scpSettings.timestamp)
% check size of data array - note whether one or two channel
disp('Data Size:')
disp(size(scanData))

t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

%% Check Waveform and Extract RMS Voltage

data = squeeze(scanData(:,:));
disp(size(data))

% Bandpass filter design
Fs = scpSettings.SampleFrequency; % Sampling Frequency
F0 = 2*1e6; % Centre frequency
width = 0.15*1e6;
Fpass1 = F0-width; % First Passband Frequency
Fpass2 = F0+width; % Second Passband Frequency

wvfms_filtered = bandpass(data', [Fpass1 Fpass2], Fs)';

bias = mean(data,2);

data_f = wvfms_filtered;
data_no_bias = data - bias;

plot(data_no_bias(50,:))
hold on
plot(data_f(50,:))
hold off

RMS = squeeze(rms(data_f - bias,2));
disp(size(RMS))

%% 

scatter(raster.zs,RMS)
xlabel('Z Position [mm]')
ylabel('RMS Voltage [V]')
x = raster.xs(1);
y = raster.ys(1);
%title(strcat('Scan Home: [x,y,z]=[',string(raster.home(1)),',',string(raster.home(2)),',',string(raster.home(3)),'] mm'))

%% To MPa
mVperMPa = 170.12; % CHECK
MPa = RMS*1e3/mVperMPa; 
