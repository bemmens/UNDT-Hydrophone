
import zaber.motion.ascii.Connection;
import zaber.motion.Units;

connection = Connection.openSerialPort('COM5');                         %CHECK
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

%% DEFINE raster
repeat = 1;
while repeat == 1
raster.home = input('raster.home [x,y,z]:'); % home position [x,y,x] in mm     % CHECK
raster.size = input('raster.size [x,y,z]:'); % [X,Y,Z] in mm                      % CHECK
raster.step = input('raster.step :'); % mm                                       % CHECK
raster.pause_time = 50/1000; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK

if min(raster.home - raster.size/2) < 0
    error('ERROR: raster.size too big')
elseif min(raster.home - raster.size/2) == 0
    warning('RASTER LIMIT = AXIS LIMIT')
end

raster.xs = (raster.home(1) - 0.5*(raster.size(1))) : raster.step : (raster.home(1) + 0.5*(raster.size(1))) ;
raster.ys = (raster.home(2) - 0.5*(raster.size(2))) : raster.step : (raster.home(2) + 0.5*(raster.size(2))) ;
raster.zs = (raster.home(3) - 0.5*(raster.size(3))) : raster.step : (raster.home(3) + 0.5*(raster.size(3))) ;

raster.xlims = [min(raster.xs),max(raster.xs)];
raster.ylims = [min(raster.ys),max(raster.ys)];
raster.zlims = [min(raster.zs),max(raster.zs)];

disp('Moving to raster.home [x,y,zMin] ...')
xAxis.moveAbsolute(raster.home(1), Units.LENGTH_MILLIMETRES)
yAxis.moveAbsolute(raster.home(2), Units.LENGTH_MILLIMETRES)
zAxis.moveAbsolute(raster.zlims(2), Units.LENGTH_MILLIMETRES)
repeat = input('Change? [n=0/y=1]: ');
end

traceScanVolume(xAxis,yAxis,zAxis,raster)

connection.close();
catch exception
    connection.close();
    rethrow(exception);
end

disp(raster)
disp('DONE - nice job!')