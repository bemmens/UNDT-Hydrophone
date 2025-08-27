function plotProgress(mVperMPa, scpSettings, scanData, xyz_switch, raster, i,j)

% scanData.ab [Na,Nb,Nsamples,Nchannels,nRepeats]
scanData_noTrigger.XY = squeeze(scanData.XY(:,:,:,1,:));
scanData_noTrigger.YZ = squeeze(scanData.YZ(:,:,:,1,:));
scanData_noTrigger.XZ = squeeze(scanData.XZ(:,:,:,1,:));

% scanData.ab [Na,Nb,Nsamples,nRepeats]
scanData_noBias.XY = scanData_noTrigger.XY - mean(scanData_noTrigger.XY,3);
scanData_noBias.YZ = scanData_noTrigger.YZ - mean(scanData_noTrigger.YZ,3);
scanData_noBias.XZ = scanData_noTrigger.XZ - mean(scanData_noTrigger.XZ,3);

% scanData.ab [Na,Nb,Nsamples,nRepeats]
Vrms.XY = mean(squeeze(rms(scanData_noBias.XY,3)),3); % rms voltage at [x,y,z0]
Vrms.YZ = mean(squeeze(rms(scanData_noBias.YZ,3)),3); % rms voltage at [x,y,z0]
Vrms.XZ = mean(squeeze(rms(scanData_noBias.XZ,3)),3); % rms voltage at [x,y,z0]

% Vrms.ab [Na,Nb]
MPa.XY = Vrms.XY*1e3/mVperMPa; 
MPa.YZ = Vrms.YZ*1e3/mVperMPa; 
MPa.XZ = Vrms.XZ*1e3/mVperMPa; 

%% Waveform Check
figure(1)
subplot(1,2,1)
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us
switch xyz_switch
    case 1
        plot(t,squeeze(scanData_noBias.XY(i,j,:)))
    case 2
        plot(t, squeeze(scanData_noBias.YZ(i,j,:)))
    case 3
        plot(t, squeeze(scanData_noBias.XZ(i,j,:)))
end
xlabel('Time [us]');
ylabel('Voltage [V]');
%% Plot 3D orthogonal views - CoPilot
subplot(1,2,2)
hold on

% XY plane
[X1, Y1] = meshgrid(raster.xs, raster.ys);
Z1 = ones(size(X1)) * raster.home(3);
surf(X1, Y1, Z1, MPa.XY', 'EdgeColor', 'none')
switch xyz_switch
    case 1
        view([0.00 90.00]) % XY
    case 2
        view([90.00 0.00]) % YZ
    case 3
        view([0.00 0.00]) % Set view for XZ plane
end

% YZ plane
[Y2, Z2] = meshgrid(raster.ys, raster.zs);
X2 = ones(size(Y2)) * raster.home(1);
surf(X2, Y2, Z2, MPa.YZ', 'EdgeColor', 'none')
switch xyz_switch
    case 1
        view([0.00 90.00]) % XY
    case 2
        view([90.00 0.00]) % YZ
    case 3
        view([0.00 0.00]) % Set view for XZ plane
end

% XZ plane
[X3, Z3] = meshgrid(raster.xs, raster.zs);
Y3 = ones(size(X3)) * raster.home(2);
surf(X3, Y3, Z3, MPa.XZ', 'EdgeColor', 'none')
switch xyz_switch
    case 1
        view([0.00 90.00]) % XY
    case 2
        view([90.00 0.00]) % YZ
    case 3
        view([0.00 0.00]) % Set view for XZ plane
end

cb = colorbar;
cb.Label.String = 'Pressure (MPa)';
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
%xlim([15,30])
%ylim([15,33])
title('Tri-Planar Scan of Acoustic Field')
subtitle('Stage Coordinates')
hold off


