% OscilloscopeBlock.m - for LibTiePie 0.5+
%
% This example performs a block mode measurement and plots the data.
%
% Find more information on http://www.tiepie.com/LibTiePie .

% Open LibTiePie and display library info if not yet opened:
%LibTiePieNeeded

% Search for devices:
LibTiePie.DeviceList.update();

% Try to open an oscilloscope with block measurement support:
clear scp;
for k = 0 : LibTiePie.DeviceList.Count - 1
    item = LibTiePie.DeviceList.getItemByIndex(k);
    if item.canOpen(DEVICETYPE.OSCILLOSCOPE)
        scp = item.openOscilloscope();
        if ismember(MM.BLOCK, scp.MeasureModes)
            break;
        else
            %clear scp;
        end
    end
end
clear item

if exist('scp', 'var')
    
      
  % Set measure mode:
    scp.MeasureMode = MM.BLOCK;

    Ch1 =scp.Channels(1);
    Ch1.Enabled = true;

    Ch2 =scp.Channels(2);
    Ch2.Enabled = false;

    % Set range:
    Ch1.Range = 0.4; % 8 V

    % Set coupling:
    Ch1.Coupling = CK.DCV; % DC Volt

    % Release reference:
    clear Ch1;
    clear Ch2;

    
    % Set sample frequency:
    scp.SampleFrequency =500e6; %

    % Set record length:
    scp.RecordLength = 65e6; 

    % Set pre sample ratio:
    scp.PreSampleRatio = 0; % 0 %
    
    % Set resolution :
    scp.Resolution=8;
    
    
    
    % Disable all channel trigger sources:
    

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

    % Plot results:
    figure(500);
    plot((1:scp.RecordLength) / scp.SampleFrequency, arData);
    axis([0 (scp.RecordLength / scp.SampleFrequency) min(darRangeMin) max(darRangeMax)]);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');

    % Close oscilloscope:
   % clear scp;
else
    error('No oscilloscope available with block measurement support!');
end