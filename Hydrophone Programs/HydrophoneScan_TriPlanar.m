%% Details
% This program will take a scan of three orthogonal slices

clear all;
close all;
clc;
fclose all;

%% Check Savefile
File_loc = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\'; % CHECK
File_name = 'DIYMK1_39'; % CHECK
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
    record_time = 0.2/1e3; % seconds                % CHECK
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
    scp.TriggerTimeOut = 5 * 1e-3; % ms -> Long delay to indicate trigger not found
    
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
scpSettings.nRepeats = 10;           % CHECK
scpSettings.timestamp = datetime;
scpSettings.scanVersion = 2; % CHECK

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
% Use scanVolumeChecker to quickly make sure that the raster parameters are
% correct without having to boot up HandyScope each time.

wavelength = 1.5; % in mm

ymin = 0;
ymax = 50;
xmin = 0;
xmax = 50;
zmin = 0;
zmax = 10;

xhome = mean([xmin,xmax]);
yhome = mean([ymin,ymax]);
zhome = mean([zmin,zmax]);

xsize = xmax - xmin;
ysize = ymax-ymin;
zsize = zmax-zmin;

raster.home = [xhome,yhome,zhome]; % home position [x,y,x] in mm     % CHECK
raster.size = [xsize ysize zsize]; % [X,Y,Z] in mm                      % CHECK

raster.home = [20,22.5,6]; % home position [x,y,x] in mm     % CHECK
raster.size = [15 15 10]; % [X,Y,Z] in mm                      % CHECK
raster.step = [0.5,0.5,0.125]*wavelength; % [dx,dy,dx] mm - must be greater than zero          % CHECK

raster.pause_time = 50/1000; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK

raster.xs = (raster.home(1) - 0.5*(raster.size(1))) : raster.step(1) : (raster.home(1) + 0.5*(raster.size(1))) ;
raster.ys = (raster.home(2) - 0.5*(raster.size(2))) : raster.step(2) : (raster.home(2) + 0.5*(raster.size(2))) ;
raster.zs = (raster.home(3) - 0.5*(raster.size(3))) : raster.step(3) : (raster.home(3) + 0.5*(raster.size(3))) ;

raster.xlims = [min(raster.xs),max(raster.xs)];
raster.ylims = [min(raster.ys),max(raster.ys)];
raster.zlims = [min(raster.zs),max(raster.zs)];

NPoints = length(raster.xs)*length(raster.ys) + length(raster.ys)*length(raster.zs) + length(raster.xs)*length(raster.zs);

disp('rater.home/size/step/Npoints:')
disp(raster.home)
disp(raster.size)
disp(raster.step)
disp(NPoints)

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


%% Define Scan Sequence
snakeCoords.XY = zeros(length(raster.xs)*length(raster.ys),3);
snakeCoords.YZ = zeros(length(raster.ys)*length(raster.zs),3);
snakeCoords.XZ = zeros(length(raster.xs)*length(raster.zs),3);

ys = raster.ys;
xs = raster.xs;
zs = flip(raster.zs); % invert z-axis: start at the bottom

% XY
index = 1;
for j = 1:length(ys)
    if xs(end) == max(xs)
        xs = flip(raster.xs);
    else
        xs = raster.xs;
    end
    for i = 1:length(xs)
        snakeCoords.XY(index,:) = [xs(i), ys(j), raster.home(3)];
        index = index +1;
    end
end

% YZ
index = 1;
for j = 1:length(zs)
    if ys(end) == max(ys)
        ys = flip(raster.ys);
    else
        ys = raster.ys;
    end
    for i = 1:length(ys)
        snakeCoords.YZ(index,:) = [raster.home(1), ys(i), zs(j)];
        index = index +1;
    end
end


% XZ
index = 1;
for j = 1:length(zs)
    if xs(end) == max(xs)
        xs = flip(raster.xs);
    else
        xs = raster.xs;
    end
    for i = 1:length(xs)
        snakeCoords.XZ(index,:) = [xs(i), raster.home(2), zs(j)];
        index = index +1;
    end
end

%% Create results struct

scanData.XY = zeros(length(raster.xs),length(raster.ys),scp.RecordLength,2,scpSettings.nRepeats); % [x,y,wvfm,chanel,nth repeat]
scanData.YZ = zeros(length(raster.ys),length(raster.zs),scp.RecordLength,2,scpSettings.nRepeats); % [y,z,wvfm,chanel,nth repeat]
scanData.XZ = zeros(length(raster.xs),length(raster.zs),scp.RecordLength,2,scpSettings.nRepeats); % [x,z,wvfm,chanel,nth repeat]

%% SCAN
disp('Scan Started')
tStart = tic;
pause('on')

prog = 0;
f = waitbar(0,'Scan Starting...');

oldCoords = raster.home;

NPointsXY = length(snakeCoords.XY(:,1));
NPointsYZ = length(snakeCoords.YZ(:,1));
NPointsXZ = length(snakeCoords.XZ(:,1));

% XY Scan
for n = 1: NPointsXY
    tStartStep = tic;
    
    % Move Sensor
    % Only engage axis if position has changed
    if snakeCoords.XY(n,1) ~= oldCoords(1)
        xAxis.moveAbsolute(snakeCoords.XY(n,1), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords.XY(n,2) ~= oldCoords(2)
        yAxis.moveAbsolute(snakeCoords.XY(n,2), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords.XY(n,3) ~= oldCoords(3)
        zAxis.moveAbsolute(snakeCoords.XY(n,3), Units.LENGTH_MILLIMETRES)
    end

    pause(raster.pause_time) % can tweak this to speed up or slow down scan: risk of shaky sensor

    % Calculate the indices for the current coordinate
    i = find(raster.xs == snakeCoords.XY(n,1));
    j = find(raster.ys == snakeCoords.XY(n,2));

    % Take measurement
    for r = 1:scpSettings.nRepeats
        
        [scp, measurement] = takeMeasOscilloscope( scp );
      
        % Store the measurement in the data array
        scanData.XY(i,j,:,:,r) = measurement;
    end

    % Admin
    oldCoords = snakeCoords.XY(n,:);

    % Progress tracking
    prog = prog + 1;
    dtStep = toc(tStartStep);
    progFrac = prog/NPoints; 
    NPointsRemaining = NPoints - prog;
    estTimeRemaining = round(NPointsRemaining*dtStep/60); % minutes
    f = waitbar((progFrac),f,strcat("Scan Running... Estimated Time Remaining: ", string(estTimeRemaining),'mins'));

end

% YZ Scan
for n = 1: NPointsYZ
    tStartStep = tic;
    
    % Move Sensor
    % Only engage axis if position has changed
    if snakeCoords.YZ(n,1) ~= oldCoords(1)
        xAxis.moveAbsolute(snakeCoords.YZ(n,1), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords.YZ(n,2) ~= oldCoords(2)
        yAxis.moveAbsolute(snakeCoords.YZ(n,2), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords.YZ(n,3) ~= oldCoords(3)
        zAxis.moveAbsolute(snakeCoords.YZ(n,3), Units.LENGTH_MILLIMETRES)
    end

    pause(raster.pause_time) % can tweak this to speed up or slow down scan: risk of shaky sensor

    % Calculate the indices for the current coordinate
    i = find(raster.ys == snakeCoords.YZ(n,2));
    j = find(raster.zs == snakeCoords.YZ(n,3));

    % Take measurement
    for r = 1:scpSettings.nRepeats
        
        [scp, measurement] = takeMeasOscilloscope( scp );
      
        % Store the measurement in the data array
        scanData.YZ(i,j,:,:,r) = measurement;
    end

    % Admin
    oldCoords = snakeCoords.YZ(n,:);

    % Progress tracking
    prog = prog + 1;
    dtStep = toc(tStartStep);
    progFrac = prog/NPoints; 
    NPointsRemaining = NPoints - prog;
    estTimeRemaining = round(NPointsRemaining*dtStep/60); % minutes
    f = waitbar((progFrac),f,strcat("Scan Running... Estimated Time Remaining: ", string(estTimeRemaining),'mins'));

end

% XZ Scan
for n = 1: NPointsXZ
    tStartStep = tic;
    
    % Move Sensor
    % Only engage axis if position has changed
    if snakeCoords.XZ(n,1) ~= oldCoords(1)
        xAxis.moveAbsolute(snakeCoords.XZ(n,1), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords.XZ(n,2) ~= oldCoords(2)
        yAxis.moveAbsolute(snakeCoords.XZ(n,2), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords.XZ(n,3) ~= oldCoords(3)
        zAxis.moveAbsolute(snakeCoords.XZ(n,3), Units.LENGTH_MILLIMETRES)
    end

    pause(raster.pause_time) % can tweak this to speed up or slow down scan: risk of shaky sensor

    % Calculate the indices for the current coordinate
    i = find(raster.xs == snakeCoords.XZ(n,1));
    j = find(raster.zs == snakeCoords.XZ(n,3));

    % Take measurement
    for r = 1:scpSettings.nRepeats
        
        [scp, measurement] = takeMeasOscilloscope( scp );
      
        % Store the measurement in the data array
        scanData.XZ(i,j,:,:,r) = measurement;
    end

    % Admin
    oldCoords = snakeCoords.XZ(n,:);

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
save(Save_String,'scanData','raster','scpSettings','snakeCoords',"-v7.3");
disp(strcat('File Saved: Data\',File_name,'.mat'));
