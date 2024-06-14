clear all
%%
import zaber.motion.ascii.Connection;
import zaber.motion.Units;

connection = Connection.openSerialPort('COM5');
try
    connection.enableAlerts();

    deviceList = connection.detectDevices();
    fprintf('Found %d devices.\n', deviceList.length);

    for i = 1:length(deviceList)
    device = deviceList(i);
    fprintf('Homing all axes of device with address %d.\n', device.getDeviceAddress());
    device.getAllAxes().home();
    end

    deviceList(1).getAxis(1).moveAbsolute(10, Units.LENGTH_MILLIMETRES) % z
    deviceList(2).getAxis(1).moveAbsolute(20, Units.LENGTH_MILLIMETRES) % x
    deviceList(3).getAxis(1).moveAbsolute(30, Units.LENGTH_MILLIMETRES) % y

    %connection.close();
catch exception
    connection.close();
    rethrow(exception);
end
