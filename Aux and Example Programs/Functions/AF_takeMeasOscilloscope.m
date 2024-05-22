function [ scp, arData, darRangeMin, darRangeMax ] = AF_takeMeasOscilloscope( scp )

    % Start measurement:
    scp.start();

    % Wait for measurement to complete:
    while ~scp.IsDataReady
        pause(10e-3) % 10 ms delay, to save CPU time.
    end

    % Get data:
    arData = scp.getData();

    % Get all channel data value ranges (which are compensated for probe gain/offset):
    clear darRangeMin;
    clear darRangeMax;
    for i = 1 : length(scp.Channels)
        [darRangeMin(i), darRangeMax(i)] = scp.Channels(i).getDataValueRange();
    end

end

