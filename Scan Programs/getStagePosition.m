%% Initialise Zaber Satges

import zaber.motion.ascii.Connection;
import zaber.motion.Units;

connection = Connection.openSerialPort('COM4');                         %CHECK
try
    connection.enableAlerts();

    deviceList = connection.detectDevices();
    fprintf('Found %d Zaber devices.\n', deviceList.length);

    for i = 1:length(deviceList)
        device = deviceList(i);
        
        axis = device.getAxis(1);
        if ~axis.isHomed()
            axis.home();
            fprintf('Homing device with address %d.\n', device.getDeviceAddress());
        end
    end

    disp('All axes zeroed.')

    zAxis = deviceList(1).getAxis(1);
    xAxis = deviceList(2).getAxis(1);
    yAxis = deviceList(3).getAxis(1);


    xaxisSettings = xAxis.getSettings();
    x = xaxisSettings.get('pos', Units.LENGTH_MILLIMETRES);
    yaxisSettings = yAxis.getSettings();
    y = yaxisSettings.get('pos', Units.LENGTH_MILLIMETRES);
    zaxisSettings = zAxis.getSettings();
    z = zaxisSettings.get('pos', Units.LENGTH_MILLIMETRES);

    disp('[x,y,z]mm:')
    disp([x,y,z])

    connection.close();
catch exception
    connection.close();
    rethrow(exception);
end