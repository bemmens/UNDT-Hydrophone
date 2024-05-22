%% Details
% Simplest code to run hydrophone scan.
% Produces a single array with the raw hydrophone data at each point in the raster.
% Warning can create LARGE files

clear all;
close all;
clc;
fclose all;


%% CONNECT TO STAGE (may have to change serial port name)
disp('Connecting to stage...')
addpath('Stage Control')
    if ~isempty(instrfind)
     fclose(instrfind);
      delete(instrfind);
end

% Specify serial port
portName = 'COM5'   ;

deviceAddress = 0; % 0 = both, 1 = X-axis, 2 = Y-axis
axisNumber = 0; 

% Set up serial object
s = serial(portName);
set(s, 'BaudRate',115200, 'DataBits',8, 'FlowControl','none',...
    'Parity','none', 'StopBits',1, 'Terminator','CR/LF'); % Don't mess with these values

fopen(s);
disp('Stage connected.')
%% Connect to HandyScope
disp('Connecting to HandyScope...')
% Open LibTiePie and display library info if not yet opened:
import LibTiePie.Const.*
import LibTiePie.Enum.*

if ~exist('LibTiePie', 'var')
  % Open LibTiePie:
  LibTiePie = LibTiePie.Library
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
    MHz = 50;
    scp.SampleFrequency = MHz*1e6; %  MHz

    % Set record length:
    scp.RecordLength = 1e3; % n Samples: max = 33553920 ~ 3e7 (67107840?)

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
    
    % Set range on each channel (V)
    scp.Channels(1).Range = 5 ; 
    scp.Channels(2).Range = 5 ;
    
    % Print oscilloscope info:
    display(scp);
    
    else
    warning('No Scope Detected')
end

%% MANUALLY SET HOME WITH STAGE

% Tune the stage with Zaber GUI Interface (must quit Matlab and restart
% once homed)
% OR
% Use the function 'AF_moveToMili(s,pos_x,pos_y)' to set the position

% move to home
AF_moveToMili(s,25,25) % s is serial pos in mm

% SET CURRENT POSITION AS HOME
[home_pos] = AF_setHome(s) ;

%% DEFINE RASTER POINTS/AREA

raster_x_size = 20; % mm 
raster_y_size = 20; % mm 
step_size = 2.5; % mm
pause_time = 0.25; % seconds - Time for motion to stop before and after measurement - Oscilliscope will wait for itself

raster_x = (home_pos(1) - 0.5*(raster_x_size*20000)) : step_size*20000 : (home_pos(1) + 0.5*(raster_x_size*20000)) ;
raster_y = (home_pos(2) - 0.5*(raster_y_size*20000)) : step_size*20000 : (home_pos(2) + 0.5*(raster_y_size*20000)) ;

N_samples = length(raster_x)*length(raster_y);
Scan_time = 2*N_samples*(pause_time*2 + scp.RecordLength/scp.SampleFrequency);    %Very approximate

display(strcat('Rasters Defined, V.Approx Scan time =',num2str(Scan_time/60,3),'min'));

%% RASTER SCAN AND MEASURE WITH HYDROPHONE
disp('Scan running...')
tic;
targetFreq = 1e6;
record=1;
pause('on')

Xn = numel(raster_x);
Yn = numel(raster_y);
tn = scp.RecordLength;
Chn = 2;

scanData = zeros(length(raster_x),length(raster_y), tn, Chn);

for ii = 1 : numel(raster_x)
    for jj = 1: numel(raster_y)
        
        if mod(ii,2)==0 % iseven -> this makes the scan snake
            raster_yROW = flip(raster_y);
        else
            raster_yROW = raster_y;
        end
      
        AF_moveToPos(s, raster_x(ii), raster_yROW(jj))
        
        pause(pause_time) % can tweak these to spped up or slow down scan
        
        [scp, arData, darRangeMin, darRangeMax] = AF_takeMeasOscilloscope( scp );

        pause(pause_time) % Redundant?

        scanData(ii,jj,:,:) = arData;
        
    end
end


% Return to the centre
AF_moveToPos(s,home_pos(1),home_pos(2));

toc

disp('Scan Complete.');

%% Saving results

File_loc = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
File_name = 'TEST BARNEY';

Save_String=strcat(File_loc,File_name,'.mat');
save(Save_String,'scanData');
disp(strcat('File Saved: Data\',File_name,'.mat'));