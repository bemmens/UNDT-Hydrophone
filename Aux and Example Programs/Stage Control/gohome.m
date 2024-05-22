function gohome(s, deviceAddress, axisNumber, home_pos)
% GOHOME instructs the XY stage to travel back to a predefined home
% position in either X, Y or both axis.
% 
% INPUTS
%       S               -   is the serial port handle for the XY stage
%       deviceAddress   -   0-2 defining which axis or axes are to go home.
%                           0 = both axes
%                           1 = x-axis
%                           2 = y-axis
%       axisNumber      -   Confusingly named, axisNumber should always = 0
%       home_pos        -   is the predefined [x,y] co-ordinate of HOME
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    if deviceAddress == 0
        % Send X-stage to inspection home command, ignoring the reply
        sendCommand(s, 1, axisNumber, ['move abs ', num2str(home_pos(1))]);
        % Send Y-stage to inspection home command, ignoring the reply
        sendCommand(s, 2, axisNumber, ['move abs ', num2str(home_pos(2))]);
    elseif deviceAddress == 1
        % Send X-stage to inspection home command, ignoring the reply
        sendCommand(s, 1, axisNumber, ['move abs ', num2str(home_pos(1))]);
    elseif deviceAddress == 2
        % Send X-stage to inspection home command, ignoring the reply
        sendCommand(s, 2, axisNumber, ['move abs ', num2str(home_pos(2))]);
    else
    end
    % Wait until axis finishes
    pollUntilIdle(s, deviceAddress, axisNumber);
    
    % Get current X position
    reply = sendCommand(s, 1, axisNumber, 'get pos');
    pos(1) = str2num(reply.data);
    % Get current Y position
    reply = sendCommand(s, 2, axisNumber, 'get pos');
    pos(2) = str2num(reply.data);
    
    disp(['Current [X,Y] position is [' num2str(pos(1)) ', ' num2str(pos(2)) '].']);
    disp(['Current [x,y] position in mm is [' num2str(pos(1).*5e-5, '%.2f') ', ' num2str(pos(2).*5e-5, '%.2f') '] mm.']);
    
catch e
    % close port first if an exception occurs
    fclose(s);
    rethrow(e)
end
end