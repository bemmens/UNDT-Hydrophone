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

%% SET OSCILLOSCOPE SETTINGS

if exist('scp', 'var')
    % Set measure mode:
    scp.MeasureMode = 2; % Block Mode

    % Set sample frequency:
    MHz = 50;     % CHECK
    scp.SampleFrequency = MHz*1e6; %  MHz

    % Set record length:
    record_time = 0.5/1e3; % ms                % CHECK
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
    scp.TriggerTimeOut = 5000 * 1e-3; % ms -> Long delay to indicate trigger not found
    
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
    scp.Channels(1).Range = 2 ;     % CHECK
    scp.Channels(2).Range = 5 ;     % CHECK
    
    else
    warning('No Scope Detected')
end

disp(strcat('Record time per measurement:',string(record_time*1e6),'us.'))

%% Initialise Zaber Satges

import zaber.motion.ascii.Connection;
import zaber.motion.Units;

connection = Connection.openSerialPort('COM5');                         %CHECK
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

raster.home = [10,10,10]; % home position [x,y,x] in mm     % CHECK

raster.size = [2,2,4]; % [X,Y,Z] in mm                      % CHECK
raster.step = 1; % mm                                       % CHECK
raster.pause_time = 50/1000; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK

raster.xs = (raster.home(1) - 0.5*(raster.size(1))) : raster.step : (raster.home(1) + 0.5*(raster.size(1))) ;
raster.ys = (raster.home(2) - 0.5*(raster.size(2))) : raster.step : (raster.home(2) + 0.5*(raster.size(2))) ;
raster.zs = (raster.home(3) - 0.5*(raster.size(3))) : raster.step : (raster.home(3) + 0.5*(raster.size(3))) ;

raster.xlims = [min(raster.xs),max(raster.xs)];
raster.ylims = [min(raster.ys),max(raster.ys)];
raster.zlims = [min(raster.zs),max(raster.zs)];

disp('Moving to raster.home ...')
xAxis.moveAbsolute(raster.home(1), Units.LENGTH_MILLIMETRES)
yAxis.moveAbsolute(raster.home(2), Units.LENGTH_MILLIMETRES)
zAxis.moveAbsolute(raster.home(3), Units.LENGTH_MILLIMETRES)
disp('DONE')

traceScanVolume(xAxis,yAxis,zAxis,raster)

NPoints = length(raster.xs)*length(raster.ys)*length(raster.zs);
step_time = 0/1e3; % s

Scan_time = NPoints*(raster.pause_time + step_time);    % Very approximate (actually wrong atm)
display(strcat('Rasters Defined, V.Approx Scan time =',num2str(Scan_time/60,3),'min'));

if min(raster.home - raster.size/2) < 0
    error('ERROR: raster.size too big')
elseif min(raster.home - raster.size/2) == 0
    warning('RASTER LIMIT = AXIS LIMIT')
end

%% Make scan snake
% Define the array to store the coordinates
snakeCoords = zeros(NPoints,3);
ys = raster.ys;
xs = raster.xs;
zs = flip(raster.zs); % invert z-axis: -ve is with gravity

% Loop through the z-axis
for k = 1:length(raster.zs)
    % Check if the y-axis movement should be reversed
    if ys(end) == max(ys)
        ys = flip(raster.ys);
    else
        ys = raster.ys;
    end
    
    % Loop through the y-axis
    for j = 1:length(ys)
        % Check if the x-axis movement should be reversed
        if xs(end) == max(xs)
            xs = flip(raster.xs);
        else
            xs = raster.xs;
        end

        % Loop through the x-axis
        for i = 1:length(xs)
            % Calculate the index for the current coordinate
            index = (k-1)*length(ys)*length(xs) + (j-1)*length(xs) + i;
            
            % Add the current coordinate to the array
            snakeCoords(index,:) = [xs(i), ys(j), zs(k)];
        end
    end
end

%% Create results struct

scanData = zeros(length(raster.xs),length(raster.ys),length(raster.zs),scp.RecordLength,2); % [x,y,z,wvfm,chanel]

%% SCAN
disp('Scan Started')
tic;
pause('on')

prog = 0;
f = waitbar(0,'Scan Running...');

oldCoords = raster.home;

for n = 1: NPoints
% Only attempt move if position has changed
    if snakeCoords(n,1) ~= oldCoords(1)
        %disp('comX')
        xAxis.moveAbsolute(snakeCoords(n,1), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords(n,2) ~= oldCoords(2)
        %disp('comY')
        yAxis.moveAbsolute(snakeCoords(n,2), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords(n,3) ~= oldCoords(3)
        zAxis.moveAbsolute(snakeCoords(n,3), Units.LENGTH_MILLIMETRES)
    end

    pause(raster.pause_time) % can tweak this to spped up or slow down scan

    % Calculate the indices for the current coordinate
    i = find(raster.xs == snakeCoords(n,1));
    j = find(raster.ys == snakeCoords(n,2));
    k = find(raster.zs == snakeCoords(n,3));

    % Take measurement
    [scp, measurement] = takeMeasOscilloscope( scp );
  
    % Store the measurement in the data array
    scanData(i,j,k,:,:) = measurement;

    % Admin
    oldCoords = snakeCoords(n,:);
    prog = prog + 1;
    f = waitbar((prog/NPoints),f,'Scan Running...');

end

% Return to home
xAxis.moveAbsolute(raster.home(1), Units.LENGTH_MILLIMETRES)
yAxis.moveAbsolute(raster.home(2), Units.LENGTH_MILLIMETRES)
zAxis.moveAbsolute(raster.home(3), Units.LENGTH_MILLIMETRES)

connection.close();
catch exception
    connection.close();
    rethrow(exception);
end

toc

close(f)
disp('Scan Complete.');

%% Saving results

scpSettings.RecordLength = scp.RecordLength;
scpSettings.SampleFrequency = scp.SampleFrequency;
scpSettings.timestamp = datetime;

disp('Saving...');
File_loc = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\'; % CHECK
File_name = '3DTest'; % CHECK

Save_String=strcat(File_loc,File_name,'.mat');
save(Save_String,'scanData','raster','scpSettings',"-v7.3");
disp(strcat('File Saved: Data\',File_name,'.mat'));

