%% Details
% Added repeat functions, overwrite protection

clear all;
close all;
clc;
fclose all;

%% Load Inputs
Input = load('ScanParametersTest');
scpSettings.scanVersion = 1;

if Input.inputVersion ~= scpSettings.scanVersion % Ensure scan, and input file are compatible
    error('Incompatible version')
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
    MHz = Input.MHz;     % CHECK
    scp.SampleFrequency = MHz*1e6; %  MHz

    % Set record length:
    record_time = Input.record_time; % seconds                % CHECK
    scp.RecordLength = scp.SampleFrequency*record_time; % n Samples: max = 33553920 ~ 3e7 (67107840?)    

    % Set pre sample ratio:
    scp.PreSampleRatio = Input.scp.PreSampleRatio; 

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
    scp.TriggerTimeOut = Input.scp.TriggerTimeOut; % ms -> Long delay to indicate trigger not found
    
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
    scp.Channels(1).Range = Input.scp.Channels(1).Range ;     % CHECK
    scp.Channels(2).Range = Input.scp.Channels(2).Range ;     % CHECK
    
    else
    warning('No Scope Detected')
end

scpSettings.RecordLength = scp.RecordLength;
scpSettings.SampleFrequency = scp.SampleFrequency;
scpSettings.nRepeats = Input.scpSettings.nRepeats;           % CHECK
scpSettings.timestamp = datetime;

disp(strcat('Record time per measurement:',string(record_time*1e6),'us.'))

disp('scp.SampleFrequency & record time [MHz,us]:')
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

connection = Connection.openSerialPort(Input.serialPort);                          %CHECK
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
% correct without having to boot up HandyScope each time.#


raster.home = Input.raster.home; % home position [x,y,x] in mm     % CHECK
raster.size = Input.raster.size; % [X,Y,Z] in mm                      % CHECK

raster.step = Input.raster.step; % [dx,dy,dx] mm - must be greater than zero          % CHECK

raster.pause_time = Input.raster.pause_time; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK

raster.xs = (raster.home(1) - 0.5*(raster.size(1))) : raster.step(1) : (raster.home(1) + 0.5*(raster.size(1))) ;
raster.ys = (raster.home(2) - 0.5*(raster.size(2))) : raster.step(2) : (raster.home(2) + 0.5*(raster.size(2))) ;
raster.zs = (raster.home(3) - 0.5*(raster.size(3))) : raster.step(3) : (raster.home(3) + 0.5*(raster.size(3))) ;

raster.xlims = [min(raster.xs),max(raster.xs)];
raster.ylims = [min(raster.ys),max(raster.ys)];
raster.zlims = [min(raster.zs),max(raster.zs)];

disp('rater.home/size/step:')
disp(raster.home)
disp(raster.size)
disp(raster.step)

cont = input('Continue? [Yes = enter, No = 0]:');
if cont == 0
    cont = 1; 
    error('Canceled.')
end

disp('Moving to raster.home ...')
xAxis.moveAbsolute(raster.home(1), Units.LENGTH_MILLIMETRES)
yAxis.moveAbsolute(raster.home(2), Units.LENGTH_MILLIMETRES)
zAxis.moveAbsolute(raster.home(3), Units.LENGTH_MILLIMETRES)
disp('DONE')

cont = input('Continue? [Yes = enter, No = 0]:');
if cont == 0
    cont = 1; 
    error('Canceled.')
end

if min(raster.home - raster.size/2) < 0
    error('ERROR: raster.size too big')
elseif min(raster.home - raster.size/2) == 0
    warning('RASTER LIMIT = AXIS LIMIT')
end

%% Make scan snake
% Define the array to store the coordinates
NPoints = length(raster.xs)*length(raster.ys)*length(raster.zs);
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

scanData = zeros(length(raster.xs),length(raster.ys),length(raster.zs),scp.RecordLength,2,scpSettings.nRepeats); % [x,y,z,wvfm,chanel,nth repeat]

%% SCAN
disp('Scan Started')
tStart = tic;
pause('on')

prog = 0;
f = waitbar(0,'Scan Starting...');

oldCoords = raster.home;

for n = 1: NPoints
    tStartStep = tic;
% Only engage axis if position has changed
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

    pause(raster.pause_time) % can tweak this to speed up or slow down scan: risk of shaky sensor

    % Calculate the indices for the current coordinate
    i = find(raster.xs == snakeCoords(n,1));
    j = find(raster.ys == snakeCoords(n,2));
    k = find(raster.zs == snakeCoords(n,3));

    for r = 1:scpSettings.nRepeats
        % Take measurement
        [scp, measurement] = takeMeasOscilloscope( scp );
      
        % Store the measurement in the data array
        scanData(i,j,k,:,:,r) = measurement;
    end

    % Admin
    oldCoords = snakeCoords(n,:);
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
File_loc = Input.File_loc;
File_name = Input.File_name;

Save_String=strcat(File_loc,File_name,'.mat');
if isfile(Save_String)
    error('Use a unique savefile name.')
else
    save(Save_String,'scanData','raster','scpSettings',"-v7.3");
    disp(strcat('File Saved: Data\',File_name,'.mat'));
end