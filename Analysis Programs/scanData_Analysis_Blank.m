%% Peak Voltage Analysis
% For reference: scanData = [x,y,z,samples,channel,nrepeats];
clear all

analysisVersion = 'Blank';

%% Load Data
folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = 'AmandaTest';
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


