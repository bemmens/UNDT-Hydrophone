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
portName = 'COM5'   ;    % CHECK

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
    MHz = 25;     % CHECK
    scp.SampleFrequency = MHz*1e6; %  MHz

    % Set record length:
    scp.RecordLength = 3e5; % n Samples: max = 33553920 ~ 3e7 (67107840?)    % CHECK
    record_time = scp.RecordLength/scp.SampleFrequency;

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
    scp.Channels(1).Range = 5 ;     % CHECK
    scp.Channels(2).Range = 5 ;     % CHECK
    
    % Print oscilloscope info:
    display(scp);
    
    else
    warning('No Scope Detected')
end

disp(strcat('Record time per measurement:',string(record_time),'s.'))
%% MANUALLY SET HOME WITH STAGE

% Tune the stage with Zaber GUI Interface (must quit Matlab and restart
% once homed)

% SET CURRENT POSITION AS HOME
%[home_pos] = AF_setHome(s) ;

% OR

% Use the function 'AF_moveToMili(s,pos_x,pos_y)' to set the position

% move to home
AF_moveToMili(s,25,25) % s is serial pos in mm     % CHECK

% SET CURRENT POSITION AS HOME
[home_pos] = AF_setHome(s) ;

%% DEFINE RASTER POINTS/AREA

raster_x_size = 50; % mm           % CHECK
raster_y_size = 50; % mm           % CHECK
step_size = 10; % mm               % CHECK
pause_time = 0.25; % seconds - Time for motion to stop before and after measurement - Oscilliscope will wait for itself     % CHECK

raster_x = (home_pos(1) - 0.5*(raster_x_size*20000)) : step_size*20000 : (home_pos(1) + 0.5*(raster_x_size*20000)) ;
raster_y = (home_pos(2) - 0.5*(raster_y_size*20000)) : step_size*20000 : (home_pos(2) + 0.5*(raster_y_size*20000)) ;

xs = (raster_x - home_pos(1))/20000;
ys = (raster_y- home_pos(2))/20000;

tn = scp.RecordLength;
Chn = 2;  % Number of chanels to be saved  
scanData = zeros(length(raster_x),length(raster_y), tn, Chn);
memReq = length(raster_x)*length(raster_y)*tn*Chn*8*1e-9;
disp(strcat('scanData requires approx:',string(memReq),'GB'))

% Save Important Variables
raster.xs = xs;
raster.ys = ys;
raster.pause = pause_time;
raster.stepsize = step_size;

N_samples = length(raster_x)*length(raster_y);
Scan_time = 2*N_samples*(pause_time*2 + scp.RecordLength/scp.SampleFrequency);    %Very approximate

display(strcat('Rasters Defined, V.Approx Scan time =',num2str(Scan_time/60,3),'min'));

%% RASTER SCAN AND MEASURE WITH HYDROPHONE
disp('Scan running...')
tic;
pause('on')

Xn = numel(raster_x);
Yn = numel(raster_y);


% Raster starts in bottom left corner, then travels +dy, reaches end travels
% +dx and then -dy until done. 
for ii = 1 : numel(raster_x)
    for jj = 1: numel(raster_y)
        
        if mod(ii,2)==0 % iseven ->  go backwards; this makes the scan snake
            raster_yROW = flip(raster_y);
            AF_moveToPos(s, raster_x(ii), raster_yROW(jj))
            pause(pause_time) % can tweak these to spped up or slow down scan
            [scp, arData, darRangeMin, darRangeMax] = AF_takeMeasOscilloscope( scp );
            pause(pause_time) % Redundant?
            scanData(ii,end+1-jj,:,:) = arData;
        else
            AF_moveToPos(s, raster_x(ii), raster_y(jj))
            pause(pause_time) % can tweak these to spped up or slow down scan
            [scp, arData, darRangeMin, darRangeMax] = AF_takeMeasOscilloscope( scp );
            pause(pause_time) % Redundant?
            scanData(ii,jj,:,:) = arData;
        end
    end
end


% Return to the centre
AF_moveToPos(s,home_pos(1),home_pos(2));

toc

disp('Scan Complete.');

%% Saving results
disp('Saving...');
File_loc = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\'; % CHECK
File_name = 'yuhan test'; % CHECK

Save_String=strcat(File_loc,File_name,'.mat');
save(Save_String,'scanData','raster','scp',"-v7.3");
disp(strcat('File Saved: Data\',File_name,'.mat'));