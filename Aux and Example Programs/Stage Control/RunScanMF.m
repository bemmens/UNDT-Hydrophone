clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XYSTAGE_RASTER_AQUIRE performs a 2D raster point scan
%
% LINEAR STAGES: Zaber ?????
% AQUISITION: HF2LI Zurich Instruments Lock-in Amplifier (device990)
% 
% BEFORE RUNNING: Plug in, via USB, and power on the Zurich Instrument and 
% linear stage detailed above. Move the stage by manual control to where
% you wish the scan to start from.
% 
% Connects to the linear stages 'COM10' and sets the starting location to
% be the current position of the stage 'home_pos'.
%
% Users can choose to calibrate the stage by going to hard home, before
% returning to 'home_pos'. This is advised when the stage is first turned
% on.
% 
% The stage steps in x & y and aquires at each position, averaging over
% 'ave' number of measurements.
% 
% SCAN PARAMETERS:
% Ax = X-axis scan length
% Ay = Y-axis scan length
% res_mm = scan resolution (mm)
%
% Copyright 2016 Dr Robert Hughes, University of Bristol.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('\\ads.bris.ac.uk\filestore\MyFiles\Staff19\phzrrh\Documents\MATLAB\Scan File')
addpath('\\ads.bris.ac.uk\filestore\MyFiles\Staff19\phzrrh\Documents\MATLAB\XYstage')
addpath('\\ads.bris.ac.uk\filestore\MyFiles\Staff19\phzrrh\Documents\MATLAB\myfuncs')

eaddress = 'robert.hughes@bristol.ac.uk';

%% CODE OPTIONS %%
save_deats = 1; % 1 = save the data
hm = 0; % 1 = go to hard home at start
showoff = 0; % 1 = scan along the full x and y length before raster begins
rast = 1; % Raster Scan? (NOT ACTIVE)
lvfig = 0; % Live ifgure plotting? (NOT ACTIVE)
send_email =  1; % 1 = send an email upon completion of scan

%% MEASUREMENT PARAMETERS %%
freq = [5e6, 6.5e6];
harms = [1, 2, 3];
no_harms = length(harms);
no_frqs = length(freq).*no_harms;
fff = 1;
if no_frqs == 1
    frq = freq;
else
    for ff=1:2:no_frqs-1
        frq(ff:ff+1) = freq.*harms(fff); % Frequency in MHz (NOW ACTIVE CONTROL)
        fff = fff+1;
    end
end
amp = [225, 200, 175, 150, 125, 100].*1e-3; % Output amplitude in mV (NOW ACTIVE CONTROL)
%amp = 100.*1e-3; 
ave = 5; % Number of averages for each measurement
% frq = [14, 14.5, 15, 15.5, 16, 16.5].*1e6;

%% SCAN PARAMETERS %%
Ax = 20; % X-length of sample in mm
Ay = 20; % Y-width of sample in mm
res_mm = 0.2; % increment in mm
dt = 0.005; % wait delay on each position

%% SAVE SETTINGS %%
folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff19\phzrrh\Documents\MATLAB\Scans';

%% START SCAN %%
% START RASTER SCAN %
res = res_mm.*20e3; % increment size in motor steps
Nx = Ax/res_mm; % Number of X-axis points
Ny = Ay/res_mm; % Number of Y-axis points
totN = Nx.*Ny;
fprintf('Total number of points %. \n', totN)
t_est = 62.*totN./100;
fprintf('Estimated inspection time %.2f s\n', t_est)

% define length axis
x = [0:res_mm:Ax-res_mm];
y = [0:res_mm:Ay-res_mm];

% CONNECT TO LINEAR STAGES %%
% Specifying serial port examples
% Windows: 'COM5'
% Linux: '/dev/ttyUSB0'goh
% Mac: '/dev/tty.usbserial-DA00FT01'
portName = 'COM3';

deviceAddress = 0; % 0 = both, 1 = X-axis, 2 = Y-axis
axisNumber = 0; % 0 = both, 1 = X-axis, 2 = Y-axis

% Set up serial object
s = serial(portName);
set(s, 'BaudRate', 115200, 'DataBits',8, 'FlowControl','none',...
    'Parity','none', 'StopBits',1, 'Terminator','CR/LF');

% Open command channels between computer and stages
fopen(s);
% set the current stage position as inpsection "home"
home_pos = sethome(s, deviceAddress, axisNumber);

%% RESET STAGES %%
if hm == 1
    % Performs a hard home calibration before returning to scan starting
    % point and begin raster scan
    sendCommand(s, deviceAddress, axisNumber, 'home'); % Go to hard-stop home
    pollUntilIdle(s, deviceAddress, axisNumber); % Wait till action complete
    
    gohome(s, deviceAddress, axisNumber, home_pos) % Go to scan home
else
end
fclose(s);

% MF scan
clear ziDAQ
% Connect to the stage
fopen(s);
%     [s, home_pos] = Connect2Stage(hm);
% Connect to the Lock-in Amplifier

ziDAQ('connect');
device = autoDetect;

filename = (['BOT1_tb1_MF_', num2str(frq(1)./1e6), '-', num2str(frq(end)./1e6),'MHz_[', num2str(Ax), ',', num2str(Ay), ']_', num2str(res_mm), 'mm']);
%filename = (['TOP_tb5_MF_', num2str(frq./1e6), 'MHz_[', num2str(Ax), ',', num2str(Ay), ']_', num2str(res_mm), 'mm']);

tic
r = rasterScan2(s, device, home_pos, x, y, res, freq, amp, ave, totN, dt, no_frqs, no_harms);
g2 = toc;

% Save
if save_deats == 1
    cd(folder)
    saveAmp = amp;
    save([filename, '.mat'], 'r', 'x', 'y', 'home_pos', 'saveAmp')
    %         clearvars r filename saveAmp xy
else
end

clear ziDAQ
fclose(s);
pause(5)

softstop(s)

message = (['Dear Dr Hughes, All scans finished. Frequencies = ', num2str(frq./1e6), ' MHz']);
subject = ('All Scans Complete');
matlabmail(eaddress, message, subject)

return

r = r';
% y = y;
figure('name', 'Imaginary')
for ii = 1:1:6
    subplot(1,6,ii)
    surf(x, y, imag(r(:,:,ii))); view(2); shading interp; colormap hot; axis equal tight
    xlabel('X-axis, mm')
    ylabel('Y-axis, mm')
    title([num2str(frq(ii)./1e6), 'MHz'])
end
% axis tight equal
figure('name', 'Phase')
surf(x, y, angle(r)); view(2); shading interp; colormap hot
xlabel('X-axis, mm')
ylabel('Y-axis, mm')
% axis tight equal

cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff19\phzrrh\Documents\MATLAB\myfuncs')
[mag, pha, h1, h2] = removebackground(x, y, r, 1);




m = 10;
n = 10;
test = medfilt2(mag, [m n]);
figure; surf(x, y, test); view(2); shading interp; colormap hot

test2 = medfilt2(pha, [m n]);
figure; surf(x, y, test2); view(2); shading interp; colormap hot

figure; surf(x, y, mag-test); view(2); shading interp; colormap hot
figure; surf(x, y, pha-test2); view(2); shading interp; colormap hot

test = fft2(mag);
figure; surf(abs(test)); view(2); shading interp; colormap gray