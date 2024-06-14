clear all
%% Initialise Zaber Satges

import zaber.motion.ascii.Connection;
import zaber.motion.Units;

connection = Connection.openSerialPort('COM5');
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

disp('Moving to raster.home')
raster.home = [10,10,10]; % home position [x,y,x] in mm
xAxis.moveAbsolute(raster.home(1), Units.LENGTH_MILLIMETRES)
yAxis.moveAbsolute(raster.home(2), Units.LENGTH_MILLIMETRES)
zAxis.moveAbsolute(raster.home(3), Units.LENGTH_MILLIMETRES)
disp('Homed')

raster.size = [3,3,3]; % [X,Y,Z] in mm           % CHECK
raster.step = 1; % mm                                 % CHECK
raster.pause_time = 50/1000; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK

raster.xs = (raster.home(1) - 0.5*(raster.size(1))) : raster.step : (raster.home(1) + 0.5*(raster.size(1))) ;
raster.ys = (raster.home(2) - 0.5*(raster.size(2))) : raster.step : (raster.home(2) + 0.5*(raster.size(2))) ;
raster.zs = (raster.home(3) - 0.5*(raster.size(3))) : raster.step : (raster.home(3) + 0.5*(raster.size(3))) ;

NPoints = length(raster.xs)*length(raster.ys)*length(raster.zs);
step_time = 0/1e3; % s

Scan_time = NPoints*(raster.pause_time +step_time);    %Very approximate
display(strcat('Rasters Defined, V.Approx Scan time =',num2str(Scan_time/60,3),'min'));

if min(raster.home - raster.size) < 0
    error('ERROR: raster.size too big')
end

%% Make scan snake
% Define the array to store the coordinates
snakeCoords = zeros(NPoints,3);

% Loop through the z-axis
for k = 1:length(raster.zs)
    % Check if the y-axis movement should be reversed
    if mod(k, 2) == 0
        ys = flip(raster.ys);
    else
        ys = raster.ys;
    end
    
    % Loop through the y-axis
    for j = 1:length(ys)
        % Check if the x-axis movement should be reversed
        if mod(j, 2) == 0
            xs = flip(raster.xs);
        else
            xs = raster.xs;
        end
        
        % Loop through the x-axis
        for i = 1:length(xs)
            % Calculate the index for the current coordinate
            index = (k-1)*length(ys)*length(xs) + (j-1)*length(xs) + i;
            
            % Add the current coordinate to the array
            snakeCoords(index,:) = [xs(i), ys(j), raster.zs(k)];
        end
    end
end


%% SCAN
disp('Scan Started')
tic;
pause('on')

prog = 0;
f = waitbar(0,'Scan Running...');

oldCoords = [0,0,0];

for n = 1: NPoints

    if snakeCoords(n,1) ~= oldCoords(1)
        xAxis.moveAbsolute(snakeCoords(n,1), Units.LENGTH_MILLIMETRES)
    elseif snakeCoords(n,2) ~= oldCoords(2)
        yAxis.moveAbsolute(snakeCoords(n,2), Units.LENGTH_MILLIMETRES)
    elseif snakeCoords(n,3) ~= oldCoords(3)
        zAxis.moveAbsolute(snakeCoords(n,3), Units.LENGTH_MILLIMETRES)
    end

    pause(raster.pause_time) % can tweak this to spped up or slow down scan

    % Take measurement

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

