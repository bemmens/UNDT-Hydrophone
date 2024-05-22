function [pos] = AF_getPos(s)

% get x pos
reply = sendCommand(s, 1, 0, 'get pos');
pos(1) = str2num(reply.data);
% get y pos
reply = sendCommand(s, 2, 0, 'get pos');
pos(2) = str2num(reply.data);
% display pos
disp(['Current [X,Y] position is [' num2str(pos(1)) ', ' num2str(pos(2)) '].']);
disp(['Current [x,y] position in mm is [' num2str(pos(1).*5e-5, '%.2f') ', ' num2str(pos(2).*5e-5, '%.2f') '] mm.']);
    
end