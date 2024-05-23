clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XYSTAGE_RASTER_AQUIRE performs a 2D raster point scan
%
% LINEAR STAGES: Zaber ?????
% AQUISITION: HF2LI Zurich Instruments Lock-in Amplifier (device990)
% 
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

%% CODE OPTIONS %%
save_deats = 1; % 1 = save the data
hm = 0; % 1 = go to hard home at start
showoff = 0; % 1 = scan along the full x and y length before raster begins
rast = 1; % Raster Scan? (NOT ACTIVE)
lvfig = 0; % Live ifgure plotting? (NOT ACTIVE)
send_email =  1; % 1 = send an email upon completion of scan

%% MEASUREMENT PARAMETERS %%
frq = 20; % Frequency in MHz (NOT ACTIVE CONTROL)
amp = 100; % Output amplitude in mV (NOT ACTIVE CONTROL)
ave = 5; % Number of averages for each measurement

%% SCAN PARAMETERS %%
Ax = 20; % X-length of sample in mm
Ay = 20; % Y-width of sample in mm
res_mm = 0.2; % increment in mm
dt = 0.005; % wait delay on each position

%% SAVE SETTINGS %%
folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff19\phzrrh\Documents\MATLAB\Scans';
filename = (['2tb6_', num2str(frq), 'MHz_[', num2str(Ax), ',', num2str(Ay), ']_', num2str(res_mm), 'mm']);

%% CONNECT TO LINEAR STAGES %%
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

%% START SCAN %%
if showoff == 1
    % Moves the stage to the extremes of the scan area before raster scan
    
    % Move forward Ax mm along X-axis
    deviceAddress = 1;
    step = Ax; % step size in mm (1mm = 20,000 steps)
    step_n = step.*20e3; % step size in motor steps
    sendCommand(s, deviceAddress, axisNumber, ['move rel ', num2str(step_n)])
    pollUntilIdle(s, deviceAddress, axisNumber);
    gohome(s, deviceAddress, axisNumber, home_pos)
    
    % Move forward Ay mm along Y-axis
    deviceAddress = 2;
    step = Ay; % step size in mm (1mm = 20,000 steps)
    step_n = step.*20e3; % step size in motor steps
    sendCommand(s, deviceAddress, axisNumber, ['move rel ', num2str(step_n)])
    pollUntilIdle(s, deviceAddress, axisNumber);
    gohome(s, deviceAddress, axisNumber, home_pos)
else
end

% START RASTER SCAN %
if rast == 1;
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
    r = zeros(length(x), length(y));
    
    % counter starting values
    xx = 1;
    yy = 1;
    
    % Connect to the Lock-in Amplifier
    clear ziDAQ
    ziDAQ('connect');
    device = autoDetect;
    
    if lvfig == 1; % Open live figure
        scrsz = get(groot, 'ScreenSize');
        h1 = figure('Position', [1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]); imagesc(x,y,abs(r)); % Open Magnitude Plot
        h2 = figure('Position', [scrsz(4)/2 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]); imagesc(x,y,angle(r)); % Open Magnitude Plot
    else
    end
    
    progress = 0;
    h = waitbar(progress, 'Please wait...', 'Name', 'Scanning...',...
        'CreateCancelBtn',...
        'setappdata(gcbf, ''canceling'',1)');
    setappdata(h, 'canceling', 0)
    tic;
    for yy = 1:1:Ny
        % Check for Cancel button press
        if getappdata(h, 'canceling')
            break
            softftop(s)
        end
        for xx = 1:1:Nx
            % Check for Cancel button press
            if getappdata(h, 'canceling')
                break 
                softftop(s)
            end
            deviceAddress = 1;
            xy(xx,yy) = raster_point(s, deviceAddress, axisNumber, res, dt);
            r(xx,yy) = simple_cmplx(device, ave);
            prcnt = progress/totN;
            waitbar(prcnt, h, ['Please wait...', sprintf('%3.1f', 100*prcnt), '% Complete'])
            progress = progress+1;
        end
        gohome(s, deviceAddress, axisNumber, home_pos)
        
        deviceAddress = 2;
        xy(xx,yy) = raster_point(s, deviceAddress, axisNumber, res, dt);
        pause(2*dt)
    end
    toc;
    delete(h)
else
end

gohome(s, 0, axisNumber, home_pos)
softstop(s)

r = r';
% y = y;

figure('name', 'Magnitude')
surf(x, y, abs(r)); view(2); shading interp; colormap hot
xlabel('X-axis, mm')
ylabel('Y-axis, mm')
% axis tight equal
figure('name', 'Phase')
surf(x, y, angle(r)); view(2); shading interp; colormap hot
xlabel('X-axis, mm')
ylabel('Y-axis, mm')
% axis tight equal

cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff19\phzrrh\Documents\MATLAB\myfuncs')
[mag, pha, h1, h2] = removebackground(x, y, r, 1);

if save_deats == 1
    cd(folder)
    close gcf; close gcf;
    save([filename, '.mat'])
%     savefig(h1, [filename, '_NMagnitude.fig'])
%     savefig(h2, [filename, '_NPhase.fig'])
else
end

eaddress = 'robert.hughes@bristol.ac.uk';
message = (['Dear Dr Hughes, All scans finished. Frequencies = ', num2str(frq./1e6), ' MHz']);
subject = ('All Scans Complete');
matlabmail(eaddress, message, subject)

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