function AF_softstop(s)
% This function allows the serial port to be reopened easily for repeated
% scans. If an error occurs whilst the serial port (s) is open it will need
% to be closed before the same script can be run again. 
%
% SOFTSTOP performs a soft stop of the serial port communication link
% between the computer and the thing being controlled (serial handle s)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Close port and clean up serial object
fclose(s);
delete(s);

% Uncomment if the it says it cannot open COM10
% delete(instrfindall); 
% clear s
end