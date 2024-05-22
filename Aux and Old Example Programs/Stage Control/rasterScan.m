function r = rasterScan(s, device, home_pos, x, y, res, freq, Amp, ave, totN, dt, no_frqs, no_harms)
    lvfig = 0;

    deviceAddress = 2;
    axisNumber = 0; % 0 = both, 1 = X-axis, 2 = Y-axis
    r = ones(length(x), length(y)).*NaN;
    
    % Set the Lock-in Amp parameters
    example_connect_config_rob(device, freq, Amp, no_harms)
    
    % counter starting values
    xx = 1;
    yy = 1;
    
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
    
    for yy = 1:1:length(y)
        % Check for Cancel button press
        if getappdata(h, 'canceling')
            break
            delete(h)
            gohome(s, 0, axisNumber, home_pos)
            softstop(s)
        end
        for xx = 1:1:length(x)
            % Check for Cancel button press
            if getappdata(h, 'canceling')
                break
                delete(h)
                gohome(s, 0, axisNumber, home_pos)
                softstop(s)
            end
            deviceAddress = 1;
            xy(xx,yy) = raster_point(s, deviceAddress, axisNumber, res, dt);
            r(xx,yy) = multiFreqCmplx(device, ave, no_frqs);
            
            prcnt = progress/totN;
            waitbar(prcnt, h, ['Please wait...', sprintf('%3.1f', 100*prcnt), '% Complete'])
            progress = progress+1;
        end
        gohome(s, deviceAddress, axisNumber, home_pos)
        
        deviceAddress = 2;
        xy(xx,yy) = raster_point(s, deviceAddress, axisNumber, res, dt);
        pause(2*dt)
    end

    if getappdata(h, 'canceling')
        return
        message = ('An error occured whilst scanning');
        reply = datestr(datetime);
        subject = (['Error Occured: ', reply]);
        matlabmail(eaddress, message, subject);
        
    end
    
    delete(h)
    gohome(s, 0, axisNumber, home_pos)

end