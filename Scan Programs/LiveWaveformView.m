%% Details
% Description

clear all;
close all;
clc;
fclose all;

%% Connect to HandyScope
disp('Connecting to HandyScope...')
% Open LibTiePie and display library info if not yet opened:
import LibTiePie.Const.*
import LibTiePie.Enum.*

if ~exist('LibTiePie', 'var')
  % Open LibTiePie:
  LibTiePie = LibTiePie.Library;
else
    disp('Library Connection Failed');
end

% Search for devices:
LibTiePie.DeviceList.update()

% Try to open an oscilloscope with block measurement support:
clear scp;
for k = 0 : LibTiePie.DeviceList.Count - 1
    item = LibTiePie.DeviceList.getItemByIndex(k);
    if item.canOpen(DEVICETYPE.OSCILLOSCOPE)
        scp = item.openOscilloscope();
        if ismember(MM.BLOCK, scp.MeasureModes)
            break;
        else
            clear scp;
            disp('Cleared Scope')
        end
    end
end
clear item
disp('Connected')
%% SET OSCILLOSCOPE SETTINGS

if exist('scp', 'var')
    % Set measure mode:
    scp.MeasureMode = 2; % Block Mode

    % Set sample frequency:
    MHz = 50;     % CHECK
    scp.SampleFrequency = MHz*1e6; %  MHz

    % Set record length:
    record_time = 100e-6; % s                % CHECK
    scp.RecordLength = scp.SampleFrequency*record_time; % n Samples: max = 33553920 ~ 3e7 (67107840?)    

    % Set pre sample ratio:
    scp.PreSampleRatio = 0; 

    % For all channels:
    for ch = scp.Channels
        % Enable channel to measure it:
        ch.Enabled = true;
        
        % Set coupling:
        ch.Coupling = 1; % DC Volt

        % Release reference:
        clear ch;
    end
    
    % Trigger settings
    % Set trigger timeout: 
    scp.TriggerTimeOut = 5; % s -> Long delay to indicate trigger not found
    
    % Disable all channel trigger sources:
    for ch = scp.Channels
        ch.Trigger.Enabled = false;
        clear ch;
    end
    
    % Setup channel trigger:
    chN = 2; % trigger channel
    chTr = scp.Channels(chN).Trigger; 
    chTr.Enabled = true;
    chTr.Kind = TK.FALLINGEDGE;
    chTr.Levels(1) = 0.5; % Trigger Level (%: 0.5=>50%)
    chTr.Hystereses(1) = 0.1; % Hystereses (%)
    
    % Release reference:
    clear chTr;
    
    % Set range on each channel (V)
    scp.Channels(1).Range = 0.5 ;     % CHECK
    scp.Channels(2).Range = 5 ;     % CHECK
    
else
    warning('No Scope Detected')
end

scpSettings.RecordLength = scp.RecordLength;
scpSettings.SampleFrequency = scp.SampleFrequency;
scpSettings.scp.PreSampleRatio = scp.PreSampleRatio; 

disp(strcat('Record time per measurement:',string(record_time*1e6),'us.'))

disp('scp.SampleFrequency & redord time [MHz,us]:')
disp(scp.SampleFrequency/1e6)
disp(record_time*1e6)

%%
refreshRate = 1*1e0; % seconds between measurements
maxRunTime = 1*60; % seconds
saveData.data = zeros(maxRunTime/refreshRate,scpSettings.RecordLength); % counter,wvfm
saveData.timestamps = zeros(1,maxRunTime/refreshRate);

t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

scpSettings.timestamp = datetime; % start time of day

%% Bandpass filter 
Fs = scpSettings.SampleFrequency; % Sampling Frequency
F0 = 2*1e6; % Centre
width = 0.5*1e6;
Fpass1 = F0-width; % First Passband Frequency
Fpass2 = F0+width; % Second Passband Frequency

%%
run = 1;
counter = 1;
tic
while run == 1

    % Take measurement
    %disp('Measuring...')
    [scp, measurement] = takeMeasOscilloscope( scp );
    elapsedTime = toc;

    wvfm = measurement(:,1);
    trigger = measurement(:,2);

    %data_f =  bandpass(wvfm, [Fpass1 Fpass2], Fs);
   
    figure(1)
    plot(t,noBias(wvfm))
    %hold on
    %plot(t,trigger)
    %plot(t,data_f)
    %hold off
    xlabel('Time [us]');
    ylabel('Voltage [V]');
    %xlim([0,20])
    ymax = 0.02;
    %ylim([-ymax*1.1,ymax*1.1])
    title(strcat('Elapsed Time:',string(elapsedTime),'s'))

    saveData.data(counter,:) = measurement(:,1);
    saveData.timestamps(counter) = elapsedTime;

    pause(refreshRate)
    counter = counter + 1;
    if counter > maxRunTime/refreshRate
        warning('saveData overflowwing!')
    end
end

%%

File_loc = 'C:\Users\Public\Documents\GitHub\UNDT-Hydrophone\DataOut\'; % CHECK
File_name = 'ImpulsonicsLiveDemo2'; % CHECK
Save_String=strcat(File_loc,File_name,'.mat');
save(Save_String,'saveData','scpSettings',"-v7.3");
disp(strcat('File Saved: Data\',File_name,'.mat'));

%% Calculate rms
load('ImpulsonicsLiveDemo2.mat')
mVpMPa = 187.96; % approx

data = saveData.data(1:end,:);
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us
timestamps = saveData.timestamps;
MPa = data*1e3/mVpMPa;
figure(1)
plot(t,MPa(1,:))
ylabel('MPa')
xlabel('Time [us]')
title('Example Noise')

bias = mean(data,2);
mVrms =  rms(data-bias,2)*1e3;
MPa_rms = mVrms/mVpMPa;
mean_bias = mean(bias);

figure(2)
scatter(squeeze(saveData.timestamps(:,1:end)),MPa_rms)
ylabel('MPa RMS')
xlabel('Time')
title('Noise over Time')

%%
% Bandpass filter design
Fs = scpSettings.SampleFrequency; % Sampling Frequency
F0 = 0.419*1e6; % Centre
width = 0.15*1e6;
Fpass1 = F0-width; % First Passband Frequency
Fpass2 = F0+width; % Second Passband Frequency

wvfms_biases = mean(data,2);

% Apply the bandpass filter
wvfms_filtered = bandpass(data', [Fpass1 Fpass2], Fs)';

filtered_rms = squeeze(rms(wvfms_filtered,2))*1e3/150;

figure(3)
plot(timestamps(:,1:end-1),filtered_rms(1:end-1,:))
hold on
plot(timestamps(:,1:end-1),MPa_rms(1:end-1,:))
hold off
ylim([0,0.14])
ylabel('MPa RMS')
xlabel('Time')
title('Signal RMS - f0=429kHz, b=150kHz')
legend('Bandpass','Raw')
