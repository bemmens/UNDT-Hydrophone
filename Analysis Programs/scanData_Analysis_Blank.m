%% Peak Voltage Analysis
% For reference: scanData = [x,y,z,samples,channel,nrepeats];
clear all

analysisVersion = '3D';

%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\Immasonic\';
file_name = 'Impulsonic_1';
path = strcat(folder_path,file_name,'.mat');
load(path)
disp('Data Timestamp:')
disp(scpSettings.timestamp)
dataSize = size(scanData);
dataSize

versionCheck = strcmp(analysisVersion,scpSettings.scanVersion);

if analysisVersion == 0
    disp('Skipping compatability check. Data dimensions may be incompatible with requested plots.')
elseif versionCheck == 0
    warning("Scan and analysis programs have mismatched versions.")
end

%% Process Data

t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

scanData_flat = squeeze(scanData);
scanData_mean = squeeze(mean(scanData_flat,3));
VPk = squeeze(max(scanData_mean,[],2));
rms = mean(abs(scanData_mean),2);

%%
figure(1)
plot(t,scanData_mean(end,:))

%%

figure(2)
plot(rms);
