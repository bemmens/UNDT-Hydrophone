function [home_pos] = AF_setHome(s)
% SETHOME sets the current stage position as HOME for later reference so
% you can travel back there at the end of scans
% 
% INPUTS
%       S - serial port handle for the XY stage
try
    pos = AF_getPos(s) ;
    home_pos = pos; % define the current position as inspection home    
catch e
    % close port first if an exception occurs
    fclose(s);
    rethrow(e)
end

end