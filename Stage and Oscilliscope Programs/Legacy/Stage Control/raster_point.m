function [pos] = raster_point(s, deviceAddress, axisNumber, dx, dt)
% RASTER_POINT performs a single relative movement to the next raster point
% 
% INPUTS
%       S               -   is the serial port handle for the XY stage
%       deviceAddress   -   0-2 defining which axis or axes are to move.
%                           0 = both axes
%                           1 = x-axis
%                           2 = y-axis
%       axisNumber      -   Confusingly named, axisNumber should always = 0
%       dx              -   is the incremental distance to move
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    % Send stage to inspection point, ignoring the reply
    sendCommand(s, deviceAddress, axisNumber, ['move rel ', num2str(dx)]);
    pollUntilIdle(s, deviceAddress, axisNumber);
    
    % Get current position
    pos = sethome(s, axisNumber);
    pause(dt)
   
catch e
    % close port first if an exception occurs
    fclose(s);
    rethrow(e)
end
end