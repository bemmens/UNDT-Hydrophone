clear all
%%

analysisVersion = 1;

MHz = 50;     % CHECK
record_time = 0.1/1e3; % seconds                % CHECK

% Set pre sample ratio:
scp.PreSampleRatio = 0; 

% Trigger settings
% Set trigger timeout: 
scp.TriggerTimeOut = 0 * 1e-3; % ms -> Long delay to physically indicate trigger not found

% Set range on each channel (V)
scp.Channels(1).Range = 2 ;     % CHECK
scp.Channels(2).Range = 5 ;     % CHECK

scpSettings.nRepeats = 5;           % CHECK

serialPort = 'COM5';    

inputVersion = 1;

%% Raster

ymin = 26;
ymax = 33;
xmin = 20;
xmax = 28;
zmin = 30;
zmax = 30;

xhome = mean([xmin,xmax]);
yhome = mean([ymin,ymax]);
zhome = mean([zmin,zmax]);

xsize = xmax - xmin;
ysize = ymax-ymin;
zsize = zmax-zmin;

raster.home = [xhome,yhome,zhome]; % home position [x,y,x] in mm     % CHECK
raster.size = [xsize ysize zsize]; % [X,Y,Z] in mm                      % CHECK

raster.home = [23.75,27,20]; % home position [x,y,x] in mm     % CHECK
raster.size = [5 5 5]; % [X,Y,Z] in mm                      % CHECK

raster.step = [1,1,1]; % [dx,dy,dx] mm - must be greater than zero          % CHECK

raster.pause_time = 50/1000; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK

%% Savefile anme and location
File_loc = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\'; % CHECK
File_name = 'repeats_test'; % CHECK

% Overwite protection
Save_String=strcat(File_loc,File_name,'.mat');
if isfile(Save_String)
    error('Use a unique savefile name.')
else
    disp(strcat('Output Loc: Data\',File_name,'.mat'));
end

%% Save Input
inputLoc = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataIn\ScanParametersTest';
save(inputLoc)
disp(strcat('Input Loc: DataIn\','ScanParametersTest','.mat'));