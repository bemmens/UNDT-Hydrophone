function [x,y,z] = getPos(xAxis,yAxis,zAxis) 

%import zaber.motion.ascii.Connection;
import zaber.motion.Units;

    xaxisSettings = xAxis.getSettings();
    x = xaxisSettings.get('pos', Units.LENGTH_MILLIMETRES);
    yaxisSettings = yAxis.getSettings();
    y = yaxisSettings.get('pos', Units.LENGTH_MILLIMETRES);
    zaxisSettings = zAxis.getSettings();
    z = zaxisSettings.get('pos', Units.LENGTH_MILLIMETRES);

end