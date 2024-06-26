%% Details
% Which uses a trigger from channel 2 to initialise the recording at each point.
% Produces a single array with the raw hydrophone data at each point in the raster.
% Warning can create LARGE files

clear all;
close all;
clc;
fclose all;

%% CONNECT TO STAGE (may have to change serial port name)
disp('Connecting to stage...')
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
    MHz = 50;     % CHECK
    scp.SampleFrequency = MHz*1e6; %  MHz

    % Set record length:
    record_time = 0.5/1e3; % ms % CHECK
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
    scp.TriggerTimeOut = 5000 * 1e-3; % ms
    
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
    chTr.Levels(1) = 0.5; % Trigger Level (%: 0.5->50%)
    chTr.Hystereses(1) = 0.1; % Hystereses (%)
    
    % Release reference:
    clear chTr;
    
    % Set range on each channel (V)
    scp.Channels(1).Range = 2 ;     % CHECK
    scp.Channels(2).Range = 5 ;     % CHECK
    
    % Print oscilloscope info:
    display(scp);
    
    else
    warning('No Scope Detected')
end

disp(strcat('Record time per measurement:',string(record_time*1e3),'ms.'))
%% MANUALLY SET HOME WITH STAGE

% Tune the stage with Zaber GUI Interface (must quit Matlab and restart
% once homed)
% OR
% Use the function 'AF_moveToMili(s,pos_x,pos_y)' to set the position

% move to home
AF_moveToMili(s,28,23) % s is serial pos in mm     % CHECK

pause(5) % to allow stage to register new position
% SET CURRENT POSITION AS HOME
[home_pos] = AF_setHome(s) ;

%% DEFINE RASTER POINTS/AREA

raster_x_size = 20; % mm           % CHECK
raster_y_size = 20; % mm           % CHECK
step_size = 1; % mm               % CHECK
pause_time = 50/1000; % ms - Time for motion to stop before and after measurement - Oscilliscope will wait for itself     % CHECK

raster_x = (home_pos(1) - 0.5*(raster_x_size*20000)) : step_size*20000 : (home_pos(1) + 0.5*(raster_x_size*20000)) ;
raster_y = (home_pos(2) - 0.5*(raster_y_size*20000)) : step_size*20000 : (home_pos(2) + 0.5*(raster_y_size*20000)) ;

xs = (raster_x - home_pos(1))/20000;
ys = (raster_y- home_pos(2))/20000;

tn = scp.RecordLength; % n samples 
Chn = 2;  % Number of chanels to be saved   % CHECK
scanData = zeros(length(raster_x),length(raster_y), tn, Chn); % scan results saved here (inc. trigger chanel)
%scanData = zeros(length(raster_x),length(raster_y), tn); % NO TRIGGER CHANEL

memReq = length(raster_x)*length(raster_y)*tn*Chn*8*1e-9;
%memReq = length(raster_x)*length(raster_y)*tn*8*1e-9; % NO TRIGGER CHANEL
disp(strcat('scanData requires approx:',string(memReq),'GB'))

% Save Raster Data
raster.xs = xs;
raster.raster_x = raster_x;
raster.raster_y = raster_y;
raster.ys = ys;
raster.pause = pause_time;
raster.stepsize = step_size;

N_samples = length(raster_x)*length(raster_y);
Scan_time = 2*N_samples*(pause_time*2 + scp.RecordLength/scp.SampleFrequency);    %Very approximate
display(strcat('Rasters Defined, V.Approx Scan time =',num2str(Scan_time/60,3),'min'));

%% RASTER SCAN AND MEASURE WITH HYDROPHONE
disp('Scan Started')
tic;
pause('on')

Xn = numel(raster_x);
Yn = numel(raster_y);


% Raster starts in bottom left corner, then travels +dy, reaches end travels
% +dx and then -dy until done. 
prog = 0;
f = waitbar(0,'Scan Running...');
for ii = 1 : numel(raster_x)
    for jj = 1: numel(raster_y)
        if mod(ii,2)==0 % iseven ->  go backwards; this makes the scan snake
            %disp([raster_x(ii),raster_yROW(jj)])
            raster_yROW = flip(raster_y);
            AF_moveToPos(s, raster_x(ii), raster_yROW(jj))
            pause(pause_time) % can tweak these to spped up or slow down scan
            [scp, arData, darRangeMin, darRangeMax] = AF_takeMeasOscilloscope( scp );
            pause(pause_time) % Redundant?
            scanData(ii,end+1-jj,:,:) = arData;
            %scanData(ii,end+1-jj,:) = arData(:,1); % NO TRIGGER MEMORY
            %disp([ii,numel(raster_y)+1-jj])
        else
            %disp([raster_x(ii),raster_y(jj)])
            AF_moveToPos(s, raster_x(ii), raster_y(jj))
            pause(pause_time) % can tweak these to spped up or slow down scan
            [scp, arData, darRangeMin, darRangeMax] = AF_takeMeasOscilloscope( scp );
            pause(pause_time) % Redundant?
            scanData(ii,jj,:,:) = arData;
            %scanData(ii,jj,:) = arData(:,1); % NO TRIGGER MEMORY
            %disp([ii,jj])
        end
        prog = prog + 1;
        f = waitbar((prog/(Xn*Yn)),f,'Scan Running...');
    end
end


% Return to the centre
AF_moveToPos(s,home_pos(1),home_pos(2));

toc

close(f)
disp('Scan Complete.');

%% Saving results

scpSettings.RecordLength = scp.RecordLength;
scpSettings.SampleFrequency = scp.SampleFrequency;
scpSettings.timestamp = datetime;

disp('Saving...');
File_loc = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\'; % CHECK
File_name = 'MatchTest20'; % CHECK

Save_String=strcat(File_loc,File_name,'.mat');
save(Save_String,'scanData','raster','scpSettings',"-v7.3");
disp(strcat('File Saved: Data\',File_name,'.mat'));