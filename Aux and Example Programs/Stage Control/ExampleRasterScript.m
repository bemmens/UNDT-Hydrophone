clear all
close all
clc

% Specifying serial port examples
% Windows: 'COM5'
% Linux: '/dev/ttyUSB0'
% Mac: '/dev/tty.usbserial-DA00FT01'
portName = 'COM3';

deviceAddress = 0; % 0 = both, 1 = X-axis, 2 = Y-axis
axisNumber = 0; 

hm = 0; % set to 1 if you want the stages to go to hard home at start
rast = 1; % Raster scan?

% Set up serial object
s = serial(portName);
set(s, 'BaudRate',115200, 'DataBits',8, 'FlowControl','none',...
    'Parity','none', 'StopBits',1, 'Terminator','CR/LF'); % Don't mess with these values

fopen(s);
home_pos = sethome(s, axisNumber); % Sets its current position as it's starting point

% If hm = 1 the stage will go to a hard stop to reset the stages
if hm == 1
    sendCommand(s, deviceAddress, axisNumber, 'home');
    pollUntilIdle(s, deviceAddress, axisNumber);
    gohome(s, deviceAddress, axisNumber, home_pos)
else
end

% Perform a raster point scan
for ii = 1:1:3 % Y-axis
    if ii > 1 % Don't move the stage to next Y position until after the first set of X-axis readings
        deviceAddress = 2; % Set deviceAddress to Y-axis (2)
        [pos(ii,jj,:)] = raster_point(s, deviceAddress, axisNumber, 10.*20e3, 0.005); % Move to next Y-co-ordinate
    else
    end
    for jj = 1:1:3 % X-axis
        deviceAddress = 1; % Set deviceAddress to X-axis (1)
        [pos(ii,jj,:)] = raster_point(s, deviceAddress, axisNumber, 10.*20e3, 0.005);
    end
    gohome(s, deviceAddress, axisNumber, home_pos) % Go back to X-axis home
end
deviceAddress = 0;
gohome(s, deviceAddress, axisNumber, home_pos)
softstop(s)
