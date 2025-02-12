%% Mode List
%{
mode 1: Limit Mode Cuboidal
mode 2: Centre-Size Mode Cuboidal
mode 3: Random XYZ Line
mode 4: Tri-Planar
mode 5: Random Cloud
%}

mode = 5;

switch mode
    case 1

        %% Limit Mode Cuboidal
        ymin = 0;
        ymax = 50;
        xmin = 0;
        xmax = 50;
        zmin = 0;
        zmax = 40;
        
        xhome = mean([xmin,xmax]);
        yhome = mean([ymin,ymax]);
        zhome = mean([zmin,zmax]);
        
        xsize = xmax - xmin;
        ysize = ymax-ymin;
        zsize = zmax-zmin;
        
        raster.home = [xhome,yhome,zhome]; % home position [x,y,x] in mm     % CHECK
        raster.size = [xsize ysize zsize]; % [X,Y,Z] in mm                      % CHECK
        
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
        
        raster.xlims = [xmin,xmax];
        raster.ylims = [ymin,ymax];
        raster.zlims = [zmin,zmax];

        NPoints = length(raster.xs)*length(raster.ys)*length(raster.zs);

        % Make scan snake
        % Define the array to store the coordinates
        snakeCoords = zeros(NPoints,3);
        ys = raster.ys;
        xs = raster.xs;
        zs = flip(raster.zs); % invert z-axis: -ve is with gravity
        
        % Loop through the z-axis
        for k = 1:length(raster.zs)
            % Check if the y-axis movement should be reversed
            if ys(end) == max(ys)
                ys = flip(raster.ys);
            else
                ys = raster.ys;
            end
            
            % Loop through the y-axis
            for j = 1:length(ys)
                % Check if the x-axis movement should be reversed
                if xs(end) == max(xs)
                    xs = flip(raster.xs);
                else
                    xs = raster.xs;
                end
        
                % Loop through the x-axis
                for i = 1:length(xs)
                    % Calculate the index for the current coordinate
                    index = (k-1)*length(ys)*length(xs) + (j-1)*length(xs) + i;
                    
                    % Add the current coordinate to the array
                    snakeCoords(index,:) = [xs(i), ys(j), zs(k)];
                end
            end
        end
        raster.scanList = snakeCoords;
        raster.NPoints = NPoints;

    case 2
        wavelength = 0.657; % in mm
        
        raster.home = [25,25,30]; % home position [x,y,x] in mm     % CHECK
        raster.size = [0,20,20]; % [X,Y,Z] in mm                      % CHECK
        raster.step = [1,1,1/8]*wavelength; % [dx,dy,dx] mm - must be greater than zero          % CHECK
        raster.pause_time = 50/1000; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK
        
        raster.xs = (raster.home(1) - 0.5*(raster.size(1))) : raster.step(1) : (raster.home(1) + 0.5*(raster.size(1))) ;
        raster.ys = (raster.home(2) - 0.5*(raster.size(2))) : raster.step(2) : (raster.home(2) + 0.5*(raster.size(2))) ;
        raster.zs = (raster.home(3) - 0.5*(raster.size(3))) : raster.step(3) : (raster.home(3) + 0.5*(raster.size(3))) ;
        
        raster.xlims = [min(raster.xs),max(raster.xs)];
        raster.ylims = [min(raster.ys),max(raster.ys)];
        raster.zlims = [min(raster.zs),max(raster.zs)];

        NPoints = length(raster.xs)*length(raster.ys)*length(raster.zs);

         % Make scan snake
        % Define the array to store the coordinates
        snakeCoords = zeros(NPoints,3);
        ys = raster.ys;
        xs = raster.xs;
        zs = flip(raster.zs); % invert z-axis: -ve is with gravity
        
        % Loop through the z-axis
        for k = 1:length(raster.zs)
            % Check if the y-axis movement should be reversed
            if ys(end) == max(ys)
                ys = flip(raster.ys);
            else
                ys = raster.ys;
            end
            
            % Loop through the y-axis
            for j = 1:length(ys)
                % Check if the x-axis movement should be reversed
                if xs(end) == max(xs)
                    xs = flip(raster.xs);
                else
                    xs = raster.xs;
                end
        
                % Loop through the x-axis
                for i = 1:length(xs)
                    % Calculate the index for the current coordinate
                    index = (k-1)*length(ys)*length(xs) + (j-1)*length(xs) + i;
                    
                    % Add the current coordinate to the array
                    snakeCoords(index,:) = [xs(i), ys(j), zs(k)];
                end
            end
        end
        raster.scanList = snakeCoords;
        raster.NPoints = NPoints;


    case 3
        wavelength = 1; % in mm
        
        raster.home = [25,25,20]; % home position [x,y,x] in mm     % CHECK
        raster.size = [50,50,40]; % [X,Y,Z] in mm                      % CHECK
        raster.step = [1,1,1]*wavelength; % [dx,dy,dx] mm - must be greater than zero          % CHECK
        raster.pause_time = 50/1000; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK
        
        raster.xlims = [raster.home(1)-raster.size(1)/2,raster.home(1)+raster.size(1)/2];
        raster.ylims = [raster.home(2)-raster.size(2)/2,raster.home(2)+raster.size(2)/2];
        raster.zlims = [raster.home(3)-raster.size(3)/2,raster.home(3)+raster.size(3)/2];

        raster.Nx = round(raster.size/raster.step);
        raster.Ny = round(raster.size/raster.step);
        raster.Nz = round(raster.size/raster.step);

        rng(0,'twister'); % Make repeatably random

        raster.xs = sort((raster.xlims(2)-raster.xlims(1)).*rand(raster.Nx,1) + raster.xlims(1));
        raster.ys = sort((raster.ylims(2)-raster.ylims(1)).*rand(raster.Ny,1) + raster.ylims(1));
        raster.zs = sort((raster.zlims(2)-raster.zlims(1)).*rand(raster.Nz,1) + raster.zlims(1));

        NPoints = length(raster.xs)*length(raster.ys)*length(raster.zs);
        raster.NPoints = NPoints;
        % Make scan snake
        % Define the array to store the coordinates
        snakeCoords = zeros(NPoints,3);
        ys = raster.ys;
        xs = raster.xs;
        zs = flip(raster.zs); % invert z-axis: -ve is with gravity
        
        % Loop through the z-axis
        for k = 1:length(raster.zs)
            % Check if the y-axis movement should be reversed
            if ys(end) == max(ys)
                ys = flip(raster.ys);
            else
                ys = raster.ys;
            end
            
            % Loop through the y-axis
            for j = 1:length(ys)
                % Check if the x-axis movement should be reversed
                if xs(end) == max(xs)
                    xs = flip(raster.xs);
                else
                    xs = raster.xs;
                end
        
                % Loop through the x-axis
                for i = 1:length(xs)
                    % Calculate the index for the current coordinate
                    index = (k-1)*length(ys)*length(xs) + (j-1)*length(xs) + i;
                    
                    % Add the current coordinate to the array
                    snakeCoords(index,:) = [xs(i), ys(j), zs(k)];
                end
            end
        end
        raster.scanList = snakeCoords;

    case 4
        wavelength = 0.657; % in mm
        
        ymin = 0;
        ymax = 50;
        xmin = 0;
        xmax = 50;
        zmin = 0;
        zmax = 10;
        
        xhome = mean([xmin,xmax]);
        yhome = mean([ymin,ymax]);
        zhome = mean([zmin,zmax]);
        
        xsize = xmax - xmin;
        ysize = ymax-ymin;
        zsize = zmax-zmin;
        
        raster.home = [xhome,yhome,zhome]; % home position [x,y,x] in mm     % CHECK
        raster.size = [xsize ysize zsize]; % [X,Y,Z] in mm                      % CHECK
        
        raster.home = [25,25,30]; % home position [x,y,x] in mm     % CHECK
        raster.size = [30 30 20]; % [X,Y,Z] in mm                      % CHECK
        raster.step = [2,2,1/2]*wavelength; % [dx,dy,dx] mm - must be greater than zero          % CHECK
        
        raster.pause_time = 50/1000; % ms - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK
        
        raster.xs = (raster.home(1) - 0.5*(raster.size(1))) : raster.step(1) : (raster.home(1) + 0.5*(raster.size(1))) ;
        raster.ys = (raster.home(2) - 0.5*(raster.size(2))) : raster.step(2) : (raster.home(2) + 0.5*(raster.size(2))) ;
        raster.zs = (raster.home(3) - 0.5*(raster.size(3))) : raster.step(3) : (raster.home(3) + 0.5*(raster.size(3))) ;
        
        raster.xlims = [min(raster.xs),max(raster.xs)];
        raster.ylims = [min(raster.ys),max(raster.ys)];
        raster.zlims = [min(raster.zs),max(raster.zs)];
        
        NPoints = length(raster.xs)*length(raster.ys) + length(raster.ys)*length(raster.zs) + length(raster.xs)*length(raster.zs);

        if min(raster.home - raster.size/2) < 0
            error('ERROR: raster.size too big')
        elseif min(raster.home - raster.size/2) == 0
            warning('RASTER LIMIT = AXIS LIMIT')
        end

        % Define Scan Sequence
        snakeCoords = struct;
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
            if ys(end) == max(ys)
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

        raster.scanList = snakeCoords.XY;
        raster.scanList = [raster.scanList;snakeCoords.YZ];
        raster.scanList = [raster.scanList;snakeCoords.XZ];
        raster.NPoints = NPoints;

    case 5
        wavelength = 1; % in mm
        
        raster.home = [25,25,20]; % home position [x,y,x] in mm     % CHECK
        raster.size = [50,50,40]; % [X,Y,Z] in mm                      % CHECK
    
        raster.pause_time = 50/1000; % s - Time for motion to stop before  measurement - Oscilliscope will wait for itself     % CHECK
        
        raster.xlims = [raster.home(1)-raster.size(1)/2,raster.home(1)+raster.size(1)/2];
        raster.ylims = [raster.home(2)-raster.size(2)/2,raster.home(2)+raster.size(2)/2];
        raster.zlims = [raster.home(3)-raster.size(3)/2,raster.home(3)+raster.size(3)/2];

        raster.NPoints = 50000;

        %rng(0,'twister'); % Make repeatably random
        raster.scanList_rand = rand(raster.NPoints+1,3);
        
        % scale to volume limits
        raster.scanList_rand(:,1) = (raster.xlims(2)-raster.xlims(1)).*raster.scanList_rand(:,1) + raster.xlims(1);
        raster.scanList_rand(:,2) = (raster.ylims(2)-raster.ylims(1)).*raster.scanList_rand(:,2) + raster.ylims(1);
        raster.scanList_rand(:,3) = (raster.zlims(2)-raster.zlims(1)).*raster.scanList_rand(:,3) + raster.zlims(1);
        
        %{ 
        lowest z to highest
        raster.scanList = flip(sortrows(raster.scanList_rand,3)); % zsnaked
        %}

        %{ 
        outwards from central point
        homestack = repmat(raster.home,[raster.NPoints,1]);
        raster.scanList_rand(:,4) = sqrt(sum((raster.scanList_rand - homestack).^2,2));
        raster.scanList = sortrows(raster.scanList_rand,4);
        %}
 
        % Move to the nearest point
        raster.scanList = zeros(raster.NPoints,3);

        prog = 0;
        f = waitbar(0,'Finding Shortest Route...');

        for n = 1:raster.NPoints
 
            tStartStep = tic;

            currentLoc = raster.scanList_rand(1,:);
            
            raster.scanList_rand(:,4) = sqrt(sum((raster.scanList_rand - currentLoc).^2,2)); % Distance to all other points
            raster.scanList_rand = sortrows(raster.scanList_rand,4);

            nextLoc = raster.scanList_rand(2,:);
            
            raster.scanList(n,1:4) = nextLoc;

            raster.scanList_rand(1,:) = []; %remove scaned point from list

            prog = prog + 1;
            dtStep = toc(tStartStep);
            progFrac = prog/raster.NPoints; 
            NPointsRemaining = raster.NPoints - prog;
            estTimeRemaining = round(NPointsRemaining*dtStep/60); % minutes
            f = waitbar((progFrac),f,strcat("Finding Shortest Route... Estimated Time Remaining: ", string(estTimeRemaining),'mins'));

        end
        close(f)
        raster.NPoints = size(raster.scanList,1);
end
disp(raster.NPoints )

%% Time Estimate
scpSettings.record_time = 100*1e-6;
scpSettings.nRepeats = 5;

measureTime = raster.pause_time + scpSettings.record_time*scpSettings.nRepeats;

speed = 26;% mm/s (max speed)
meanDX = mean(abs(gradient(raster.scanList(:,3)))); %mm
moveTime = meanDX/speed;

timePerPoint = moveTime + measureTime;
scanTime = raster.NPoints*timePerPoint/60;
%% Plot Scan

cmap = colormap(jet(raster.NPoints));
scatter3(raster.scanList(:,1),raster.scanList(:,2),raster.scanList(:,3),ones(1,raster.NPoints)*3,cmap,"filled")
%plot3(raster.scanList(:,1),raster.scanList(:,2),raster.scanList(:,3))

xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title('Raster Pattern')
cb = colorbar;
cb.Label.String = 'Scan Progress';