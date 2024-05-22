function [s, home_pos] = Connect2Stage(hm)
% CONNECT TO LINEAR STAGES %%
% Specifying serial port examples
% Windows: 'COM5'
% Linux: '/dev/ttyUSB0'goh
% Mac: '/dev/tty.usbserial-DA00FT01'
portName = 'COM10';

deviceAddress = 0; % 0 = both, 1 = X-axis, 2 = Y-axis
axisNumber = 0; % 0 = both, 1 = X-axis, 2 = Y-axis

% Set up serial object
s = serial(portName);
set(s, 'BaudRate', 115200, 'DataBits',8, 'FlowControl','none',...
    'Parity','none', 'StopBits',1, 'Terminator','CR/LF');

% Open command channels between computer and stages
fopen(s);
% set the current stage position as inpsection "home"
home_pos = sethome(s, deviceAddress, axisNumber);

%% RESET STAGES %%
if hm == 1
    % Performs a hard home calibration before returning to scan starting
    % point and begin raster scan
    sendCommand(s, deviceAddress, axisNumber, 'home'); % Go to hard-stop home
    pollUntilIdle(s, deviceAddress, axisNumber); % Wait till action complete
    
    gohome(s, deviceAddress, axisNumber, home_pos) % Go to scan home
else
end

end
