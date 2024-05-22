function AF_moveToPos(s,posX,posY)

posX = posX*20000 ;
posY = posY*20000 ;
% move in x
sendCommand(s, 1, 0, ['move abs ', num2str(posX)])
% wait for stage to finish
pollUntilIdle(s, 1, 0);
% move in y
sendCommand(s, 2, 0, ['move abs ', num2str(posY)])
% wait for stage to finish
pollUntilIdle(s, 1, 0);

end