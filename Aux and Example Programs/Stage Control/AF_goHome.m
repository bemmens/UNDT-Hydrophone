function AF_goHome(s, deviceAddress, home_pos)
% GOHOME instructs the XY stage to travel back to a predefined home
% position in either X, Y or both axis.
% 
% INPUTS
%       S -   is the serial port handle for the XY stage
%       deviceAddress   -   0-2 defining which axis or axes are to go home.
%                           0 = both axes
%                           1 = x-axis
%                           2 = y-axis
%       home_pos -   is the predefined [x,y] co-ordinate of HOME

try
    if deviceAddress == 0
        % Send X-stage to inspection home command, ignoring the reply
        sendCommand(s, 1, 0, ['move abs ', num2str(home_pos(1))]);
        % Send Y-stage to inspection home command, ignoring the reply
        sendCommand(s, 2, 0, ['move abs ', num2str(home_pos(2))]);
    elseif deviceAddress == 1
        % Send X-stage to inspection home command, ignoring the reply
        sendCommand(s, 1, 0, ['move abs ', num2str(home_pos(1))]);
    elseif deviceAddress == 2
        % Send X-stage to inspection home command, ignoring the reply
        sendCommand(s, 2, 0, ['move abs ', num2str(home_pos(2))]);
    else
    end
    % Wait until axis finishes
    pollUntilIdle(s, deviceAddress, 0);    
    
catch e
    % close port first if an exception occurs
    fclose(s);
    rethrow(e)
end

end