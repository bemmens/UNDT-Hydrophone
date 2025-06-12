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
    record_time = 100/1e6; % s                % CHECK
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
    scp.TriggerTimeOut = 1; % s -> Long delay to indicate trigger not found
    
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
    scp.Channels(1).Range = 0.03 ;     % CHECK
    scp.Channels(2).Range = 5 ;     % CHECK
    scp.Channels(3).Range = 0.6 ;     % CHECK
%     
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
refreshRate = 1; % seconds (1s reccomended min)

t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

scpSettings.timestamp = datetime; % start time of day

%% Bandpass filter 
Fs = scpSettings.SampleFrequency; % Sampling Frequency
F0 = 420*1e3; % Centre
width = 0.3*F0;
Fpass1 = F0-width; % First Passband Frequency
Fpass2 = F0+width; % Second Passband Frequency

%%
run = 1;

while run == 1

    % Take measurement
    %disp('Measuring...')
    [scp, measurement] = takeMeasOscilloscope( scp );

    wvfm = measurement(:,1);
    trig = measurement(:,2);
    input = measurement(:,3);

    MPa = wvfm*1e3/150;
    MPaRMS = rms(wvfm)*1e3/150;

    wvfm_f = bandpass(wvfm, [Fpass1 Fpass2], Fs);
    MPaRMS_f = rms(wvfm_f)*1e3/150;

    figure(1)
%     bandpass(wvfm, [Fpass1 Fpass2], Fs)
    [p,f] = pspectrum(MPa/1e6,Fs);
    plot(f,log(p))
    xlabel('Frequency [MHz]')
    ylabel('Power [log(p)]')
    title('Power Spectrum')

    figure(2)
    plot(t,wvfm*1e3,':')
    hold on
    plot(t,wvfm_f*1e3)
%     plot(t,input)
    hold off
    xlabel('Time (us)')
    ylabel('Amplitude (mV)')
    title(strcat('Filtered RMS Pressure: ',string(MPaRMS_f),' MPa'))

    pause(refreshRate/2)

end

%%
saveData.wvfm = wvfm;

saveData.bandpass = [[Fpass1 Fpass2], Fs];

File_loc = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\'; % CHECK
File_name = 'test'; % CHECK
Save_String=strcat(File_loc,File_name,'.mat');
save(Save_String,'saveData','scpSettings',"-v7.3");
disp(strcat('File Saved: Data\',File_name,'.mat'));


