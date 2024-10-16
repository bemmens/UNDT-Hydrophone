%% Define Scan Volume

raster.home = [xhome,yhome,zhome]; % home position [x,y,x] in mm     % CHECK
raster.size = [xsize ysize zsize]; % [X,Y,Z] in mm                      % CHECK

%raster.home = [23.75,27,20]; % home position [x,y,x] in mm     % CHECK
%raster.size = [10 0 40]; % [X,Y,Z] in mm                      % CHECK

raster.step = [1,1,1]; % [dx,dy,dx] mm - must be greater than zero          % CHECK

raster.pause_time = 50/1000; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK

raster.xs = (raster.home(1) - 0.5*(raster.size(1))) : raster.step(1) : (raster.home(1) + 0.5*(raster.size(1))) ;
raster.ys = (raster.home(2) - 0.5*(raster.size(2))) : raster.step(2) : (raster.home(2) + 0.5*(raster.size(2))) ;
raster.zs = (raster.home(3) - 0.5*(raster.size(3))) : raster.step(3) : (raster.home(3) + 0.5*(raster.size(3))) ;

raster.xlims = [min(raster.xs),max(raster.xs)];
raster.ylims = [min(raster.ys),max(raster.ys)];
raster.zlims = [min(raster.zs),max(raster.zs)];

%% Define Scan Sequence
NPoints = length(raster.xs)*length(raster.ys) + length(raster.ys)*length(raster.zs) + length(raster.xs)*length(raster.zs);

snakeCoords.XY = zeros(length(raster.xs)*length(raster.ys),3);
snakeCoords.YZ = zeros(length(raster.ys)*length(raster.zs),3);
snakeCoords.XZ = zeros(length(raster.xs)*length(raster.zs),3);

ys = raster.ys;
xs = raster.xs;
zs = flip(raster.zs); % invert z-axis: start at the bottom

% XY
index = 1;
for j = 1:length(ys)
    if xs(end) == max(xs)
        xs = flip(raster.xs);
    else
        xs = raster.xs;
    end
    for i = 1:length(xs)
        snakeCoords.XY(index,:) = [xs(i), ys(j), raster.home(3)];
        index = index +1;
    end
end

% YZ
index = 1;
for j = 1:length(zs)
    if ys(end) == max(xs)
        ys = flip(raster.ys);
    else
        ys = raster.ys;
    end
    for i = 1:length(ys)
        snakeCoords.YZ(index,:) = [raster.home(1), ys(i), zs(j)];
        index = index +1;
    end
end


% XZ
index = 1;
for j = 1:length(zs)
    if xs(end) == max(xs)
        xs = flip(raster.xs);
    else
        xs = raster.xs;
    end
    for i = 1:length(xs)
        snakeCoords.XZ(index,:) = [xs(i), raster.home(2), zs(j)];
        index = index +1;
    end
end

%% Create results struct

scp.RecordLength = 100;
scpSettings.nRepeats = 3;

scanData.XY = zeros(length(raster.xs),length(raster.ys),scp.RecordLength,2,scpSettings.nRepeats); % [x,y,wvfm,chanel,nth repeat]
scanData.YZ = zeros(length(raster.ys),length(raster.zs),scp.RecordLength,2,scpSettings.nRepeats); % [y,z,wvfm,chanel,nth repeat]
scanData.XZ = zeros(length(raster.xs),length(raster.zs),scp.RecordLength,2,scpSettings.nRepeats); % [x,z,wvfm,chanel,nth repeat]


%% SCAN
disp('Scan Started')
tStart = tic;
pause('on')

prog = 0;
f = waitbar(0,'Scan Starting...');

oldCoords = raster.home;

NPointsXY = length(snakeCoords.XY(:,1));
NPointsYZ = length(snakeCoords.YZ(:,1));
NPointsXZ = length(snakeCoords.XZ(:,1));

% XY Scan
for n = 1: NPointsXY
    tStartStep = tic;
    
    % Move Sensor
    % Only engage axis if position has changed
    if snakeCoords.XY(n,1) ~= oldCoords(1)
        xAxis.moveAbsolute(snakeCoords.XY(n,1), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords.XY(n,2) ~= oldCoords(2)
        yAxis.moveAbsolute(snakeCoords.XY(n,2), Units.LENGTH_MILLIMETRES)
    end

    pause(raster.pause_time) % can tweak this to speed up or slow down scan: risk of shaky sensor

    % Calculate the indices for the current coordinate
    i = find(raster.xs == snakeCoords.XY(n,1));
    j = find(raster.ys == snakeCoords.XY(n,2));

    % Take measurement
    for r = 1:scpSettings.nRepeats
        
        [scp, measurement] = takeMeasOscilloscope( scp );
      
        % Store the measurement in the data array
        scanData.XY(i,j,:,:,r) = measurement;
    end

    % Admin
    oldCoords = snakeCoords.XY(n,:);

    % Progress tracking
    prog = prog + 1;
    dtStep = toc(tStartStep);
    progFrac = prog/NPoints; 
    NPointsRemaining = NPoints - prog;
    estTimeRemaining = round(NPointsRemaining*dtStep/60); % minutes
    f = waitbar((progFrac),f,strcat("Scan Running... Estimated Time Remaining: ", string(estTimeRemaining),'mins'));

end

% YZ Scan
for n = 1: NPointsYZ
    tStartStep = tic;
    
    % Move Sensor
    % Only engage axis if position has changed
    if snakeCoords.YZ(n,2) ~= oldCoords(2)
        yAxis.moveAbsolute(snakeCoords.YZ(n,2), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords.YZ(n,3) ~= oldCoords(3)
        zAxis.moveAbsolute(snakeCoords.YZ(n,3), Units.LENGTH_MILLIMETRES)
    end

    pause(raster.pause_time) % can tweak this to speed up or slow down scan: risk of shaky sensor

    % Calculate the indices for the current coordinate
    i = find(raster.ys == snakeCoords(n,2));
    j = find(raster.zs == snakeCoords(n,3));

    % Take measurement
    for r = 1:scpSettings.nRepeats
        
        [scp, measurement] = takeMeasOscilloscope( scp );
      
        % Store the measurement in the data array
        scanData.YZ(i,j,:,:,r) = measurement;
    end

    % Admin
    oldCoords = snakeCoords.YZ(n,:);

    % Progress tracking
    prog = prog + 1;
    dtStep = toc(tStartStep);
    progFrac = prog/NPoints; 
    NPointsRemaining = NPoints - prog;
    estTimeRemaining = round(NPointsRemaining*dtStep/60); % minutes
    f = waitbar((progFrac),f,strcat("Scan Running... Estimated Time Remaining: ", string(estTimeRemaining),'mins'));

end

% XZ Scan
for n = 1: NPointsXZ
    tStartStep = tic;
    
    % Move Sensor
    % Only engage axis if position has changed
    if snakeCoords.XZ(n,1) ~= oldCoords(1)
        xAxis.moveAbsolute(snakeCoords.XY(n,1), Units.LENGTH_MILLIMETRES)
    end
    if snakeCoords.XZ(n,3) ~= oldCoords(3)
        zAxis.moveAbsolute(snakeCoords.XZ(n,3), Units.LENGTH_MILLIMETRES)
    end

    pause(raster.pause_time) % can tweak this to speed up or slow down scan: risk of shaky sensor

    % Calculate the indices for the current coordinate
    i = find(raster.xs == snakeCoords.XZ(n,1));
    j = find(raster.zs == snakeCoords.XZ(n,3));

    % Take measurement
    for r = 1:scpSettings.nRepeats
        
        [scp, measurement] = takeMeasOscilloscope( scp );
      
        % Store the measurement in the data array
        scanData.XZ(i,j,:,:,r) = measurement;
    end

    % Admin
    oldCoords = snakeCoords.XZ(n,:);

    % Progress tracking
    prog = prog + 1;
    dtStep = toc(tStartStep);
    progFrac = prog/NPoints; 
    NPointsRemaining = NPoints - prog;
    estTimeRemaining = round(NPointsRemaining*dtStep/60); % minutes
    f = waitbar((progFrac),f,strcat("Scan Running... Estimated Time Remaining: ", string(estTimeRemaining),'mins'));

end


