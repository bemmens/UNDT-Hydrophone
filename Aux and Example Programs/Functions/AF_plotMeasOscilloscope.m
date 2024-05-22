function AF_plotMeasOscilloscope(scp, arData, darRangeMin, darRangeMax)

    % Plot results:
    figure(500);
    plot((1:scp.RecordLength) / scp.SampleFrequency, arData);
    axis([0 (scp.RecordLength / scp.SampleFrequency) min(darRangeMin) max(darRangeMax)]);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');

    % Close oscilloscope:
    clear scp;
    
end