function [home_pos] = AF_setHome(s)
% SETHOME sets the current stage position as HOME for later reference so
% you can travel back there at the end of scans
% 
% INPUTS
%       S               -   is the serial port handle for the XY stage
%       axisNumber      -   Confusingly named, axisNumber should always = 0
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    % Find the stages current position and set it to be "home"
    % Get current X position
    reply = sendCommand(s, 1, 0, 'get pos');
    pos(1) = str2num(reply.data);
    % Get current Y position
    reply = sendCommand(s, 2, 0, 'get pos');
    pos(2) = str2num(reply.data);
    disp(['Current [X,Y] position is [' num2str(pos(1)) ', ' num2str(pos(2)) '].']);
    disp(['Current [x,y] position in mm is [' num2str(pos(1)./20e3, '%.2f') ', ' num2str(pos(2)./20e3, '%.2f') '] mm.']);
    
    home_pos = pos; % define the current position as inspection home
    
catch e
    % close port first if an exception occurs
    fclose(s);
    rethrow(e)
end
end