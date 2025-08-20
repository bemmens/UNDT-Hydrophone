%% DEFINE raster
% Use scanVolumeChecker to quickly make sure that the raster parameters are
% correct without having to boot up HandyScope each time.

c_water = 1450; % speed of sound m/s
Hz = 2e6; % CHECK
wavelength = c_water*1e3/Hz; % in mm

raster.start = [6.6819   15.3985   25.0000-15/2]; % home position [x,y,z] in mm     % CHECK
raster.end = [0 0 0];
raster.resolution = 0.5e0*wavelength; % [dx,dy,dz] mm - must be greater than zero          % CHECK
raster.length = norm(raster.end-raster.start);
NPoints = round(raster.length/raster.resolution);
raster.xs = linspace(raster.start(1),raster.end(1),NPoints);
raster.ys = linspace(raster.start(2),raster.end(2),NPoints);
raster.zs = linspace(raster.start(3),raster.end(3),NPoints);

raster.pause_time = 20/1000; % s - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK

