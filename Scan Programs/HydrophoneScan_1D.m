%% Details
% This program will take a scan of three orthogonal slices

clear all;
close all;
clc;
fclose all;

%% Check Savefile
File_loc = 'C:\Users\Public\Documents\GitHub\UNDT-Hydrophone\DataOut\'; % CHECK
File_name = '1DTest'; % CHECK
Save_String = strcat(File_loc,File_name,'.mat');

if isfile(Save_String)
    warning('Use a unique savefile name.')
    check = input('Continue anyway? [Enter = Yes / 0 = No]:');
    if check == 0
        error('Canceled')
    end
end

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

%% SET OSCILLOSCOPE SETTINGS

if exist('scp', 'var')
    % Set measure mode:
    scp.MeasureMode = 2; % Block Mode

    % Set sample frequency:
    MHz = 50;     % CHECK
    scp.SampleFrequency = MHz*1e6; %  MHz

    % Set record length:
    record_time = 100/1e6; % seconds                % CHECK
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
    scp.TriggerTimeOut = 0; % s -> Long delay to indicate trigger not found
    
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

% Save aprameters for analysis
scpSettings.RecordLength = scp.RecordLength;
scpSettings.SampleFrequency = scp.SampleFrequency;
scpSettings.nRepeats = 1;           % CHECK
scpSettings.timestamp = datetime;
scpSettings.scanVersion = 3; % CHECK

disp(strcat('Record time per measurement:',string(record_time*1e6),'us.'))

disp('scp.SampleFrequency & redord time [MHz,us]:')
disp(scp.SampleFrequency/1e6)
disp(record_time*1e6)

cont = input('Continue? [Yes = enter, No = 0]:');
if cont == 0
    cont = 1; 
    error('Canceled.')
end

%% Initialise Zaber Satges

import zaber.motion.ascii.Connection;
import zaber.motion.Units;

connection = Connection.openSerialPort('COM4');                         %CHECK
try
    connection.enableAlerts();

    deviceList = connection.detectDevices();
    fprintf('Found %d Zaber devices.\n', deviceList.length);

    for i = 1:length(deviceList)
        device = deviceList(i);
        
        axis = device.getAxis(1);
        if ~axis.isHomed()
            axis.home();
            fprintf('Homing device with address %d.\n', device.getDeviceAddress());
        end
    end

    disp('All axes zeroed.')

    zAxis = deviceList(1).getAxis(1);
    xAxis = deviceList(2).getAxis(1);
    yAxis = deviceList(3).getAxis(1);

%% DEFINE raster
% Use scanVolumeChecker to quickly make sure that the raster parameters are
% correct without having to boot up HandyScope each time.

c_water = 1450; % speed of sound m/s
Hz = 2e6; % CHECK
wavelength = c_water*1e3/Hz; % in mm

raster.start = [6.6819   15.3985   25.0000]; % home position [x,y,z] in mm     % CHECK
raster.end = [6.6819   15.3985   25.0000-15];
raster.resolution = 0.5^3*wavelength; % [dx,dy,dz] mm - must be greater than zero          % CHECK
raster.length = norm(raster.end-raster.start);
NPoints = round(raster.length/raster.resolution);
raster.xs = linspace(raster.start(1),raster.end(1),NPoints);
raster.ys = linspace(raster.start(2),raster.end(2),NPoints);
raster.zs = linspace(raster.start(3),raster.end(3),NPoints);

raster.pause_time = 20/1000; % s - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK

raster.xlims = [min(raster.xs),max(raster.xs)];
raster.ylims = [min(raster.ys),max(raster.ys)];
raster.zlims = [min(raster.zs),max(raster.zs)];

disp('rater.start/end/resolution/Npoints:')
disp(raster.start)
disp(raster.end)
disp(raster.resolution)
disp(NPoints)

cont = input('Continue? [Yes = enter, No = 0]:');
if cont == 0
    cont = 1; 
    error('Canceled.')
end

disp('Moving to raster.start ...')
xAxis.moveAbsolute(raster.start(1), Units.LENGTH_MILLIMETRES)
yAxis.moveAbsolute(raster.start(2), Units.LENGTH_MILLIMETRES)
zAxis.moveAbsolute(raster.start(3), Units.LENGTH_MILLIMETRES)
disp('DONE')

cont = input('Continue? [Yes = enter, No = 0]:');
if cont == 0
    cont = 1; 
    error('Canceled.')
end

%% Create results struct

scanData = zeros(length(raster.xs),scp.RecordLength); % [pos,wvfm]

%% SCAN
disp('Scan Started')
tStart = tic;
pause('on')

prog = 0;
f = waitbar(0,'Scan Starting...');

% Scan
for n = 1: NPoints
    tStartStep = tic;
    
    % Move Sensor
    xAxis.moveAbsolute(raster.xs(n), Units.LENGTH_MILLIMETRES)
    yAxis.moveAbsolute(raster.ys(n), Units.LENGTH_MILLIMETRES)
    zAxis.moveAbsolute(raster.zs(n), Units.LENGTH_MILLIMETRES)

    pause(raster.pause_time) % can tweak this to speed up or slow down scan: risk of shaky sensor

    [scp, measurement] = takeMeasOscilloscope( scp );
  
    % Store the measurement in the data array
    scanData(n,:) = measurement(:,1);

    % Progress tracking
    prog = prog + 1;
    dtStep = toc(tStartStep);
    progFrac = prog/NPoints; 
    NPointsRemaining = NPoints - prog;
    estTimeRemaining = round(NPointsRemaining*dtStep/60); % minutes
    f = waitbar((progFrac),f,strcat("Scan Running... Estimated Time Remaining: ", string(estTimeRemaining),'mins'));
end

connection.close();
catch exception
    connection.close();
    rethrow(exception);
end

raster.scanDuration = toc(tStart);

close(f)
disp('Scan Complete');

%% Saving results

disp('Saving...');
save(Save_String,'scanData','raster','scpSettings',"-v7.3");
disp(strcat('File Saved: DataOut\',File_name,'.mat'));