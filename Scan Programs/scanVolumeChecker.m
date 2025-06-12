%% Use this to check quickly test your scan volume
% All inputs on command line
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


%% DEFINE raster

% make sure the full range of is movement is within the limits of the
% apparatus

ymin = 25;
ymax = 25;
xmin = 25;
xmax = 25;
zmin = 40-20;
zmax = 40;

xhome = mean([xmin,xmax]);
yhome = mean([ymin,ymax]);
zhome = mean([zmin,zmax]);

xsize = xmax - xmin;
ysize = ymax-ymin;
zsize = zmax-zmin;

raster.home = [xhome,yhome,zhome]; % home position [x,y,x] in mm     % CHECK
raster.size = [xsize ysize zsize]; % [X,Y,Z] in mm                      % CHECK

raster.home = [26.5,28.5,30]; % home position [x,y,x] in mm     % CHECK
raster.size = [30 30 5]; % [X,Y,Z] in mm - max [50,50,40]                  % CHECK

raster.step = [1,1,1]; % mm  [dx,dy,dz]                      % CHECK
raster.pause_time = 50/1000; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK

if min(raster.home - raster.size/2) < 0
    error('ERROR: raster.size too big')
elseif min(raster.home - raster.size/2) == 0
    warning('RASTER LIMIT = AXIS LIMIT')
end

raster.xs = (raster.home(1) - 0.5*(raster.size(1))) : raster.step(1) : (raster.home(1) + 0.5*(raster.size(1))) ;
raster.ys = (raster.home(2) - 0.5*(raster.size(2))) : raster.step(2) : (raster.home(2) + 0.5*(raster.size(2))) ;
raster.zs = (raster.home(3) - 0.5*(raster.size(3))) : raster.step(3) : (raster.home(3) + 0.5*(raster.size(3))) ;

raster.xlims = [min(raster.xs),max(raster.xs)];
raster.ylims = [min(raster.ys),max(raster.ys)];
raster.zlims = [min(raster.zs),max(raster.zs)];

% Move to x,y centre, z bottom
disp('Moving to raster.home [x,y,zMin] ...')
xAxis.moveAbsolute(raster.home(1), Units.LENGTH_MILLIMETRES)
yAxis.moveAbsolute(raster.home(2), Units.LENGTH_MILLIMETRES)
zAxis.moveAbsolute(raster.zlims(2), Units.LENGTH_MILLIMETRES)

cont = input('Continue? [Yes = enter, No = 0]:');
if cont == 0
    cont = 1; 
    error('Canceled.')
end

% UNCOMMENT TO TRACE SCAN VOLUME
traceScanVolume(xAxis,yAxis,zAxis,raster)   

% UNCOMMENT TO RETURN TO TRUE HOME
%disp('Returning to raster.home [x0,y0,z0]...')
%xAxis.moveAbsolute(raster.home(1), Units.LENGTH_MILLIMETRES)
%yAxis.moveAbsolute(raster.home(2), Units.LENGTH_MILLIMETRES)
%zAxis.moveAbsolute(raster.home(3), Units.LENGTH_MILLIMETRES)

connection.close();
catch exception
    connection.close();
    rethrow(exception);
end

disp(raster)
disp('DONE - nice job!')