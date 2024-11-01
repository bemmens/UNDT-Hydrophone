function [scp, arData] = takeMeasOscilloscope( scp )

    % Start measurement:
    scp.start();

    % Wait for measurement to complete:
    while ~scp.IsDataReady
        pause(10e-3) % 10 ms delay, to save CPU time.
    end

    % Get data:
    arData = scp.getData();

end

