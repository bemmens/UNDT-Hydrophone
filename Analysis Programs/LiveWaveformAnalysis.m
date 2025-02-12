File_loc = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\'; % CHECK
File_name = 'TankConnectorMk5_14'; % CHECK
Save_String=strcat(File_loc,File_name,'.mat');
load(Save_String)

wvfms_raw = saveData.data'; 
wvfms = wvfms_raw(:, any(wvfms_raw, 1)); % remove empmty recordings
timestamps = saveData.timestamps(:,any(saveData.timestamps,1));

%%
centred = noBias(wvfms);
mean_wvfm = mean(centred,2);
mag = abs(centred);
mean_mag = mean(mag);

t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

%%


figure(1)
subplot(2,1,1)
plot(t,wvfms(:,1))
xlabel('Time [us]')
ylabel('Volts')
title('Waveform @ t = 0mins')

subplot(2,1,2)
plot(t,mag(:,1))
xlabel('Time [us]')
ylabel('abs(Volts)')
title('Waveform @ t = 0mins')

figure(2)
plot(timestamps/60,mean_mag/mean_mag(1))
xlabel('Time [Minutes]')
ylabel('Fraction of V(t=0)')
title('Mean(abs(waveform voltage) over ~2hrs ')

%%
figure(5)
mVperMPa = 153.23; % CHECK
plot(t,mean_wvfm*1e3/mVperMPa)
xlabel('Time [us]')
%ylabel('Displacement [Volts]')
ylabel('Pressure [MPa]')
%ylim([-0.8,0.8])
title('Time averaged waveform 30 cycles')

%%
src.Value = 1;
figure(3)
slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', size(wvfms, 2), 'sliderstep', [1/size(wvfms, 2),1/10] ,'Value', 1, 'Position', [0.4*figure(3).Position(3), 0.05*figure(3).Position(4), 0.2*figure(3).Position(3), 0.04*figure(3).Position(4)]);
slider_label = uicontrol('Style', 'text', 'Position', [0.4*figure(3).Position(3), 0.1*figure(3).Position(4), 0.2*figure(3).Position(3), 0.04*figure(3).Position(4)], 'String', 'Time Slider');
slider.Callback = @(src, event) updatePlot(src, t, wvfms,saveData);
updatePlot(src, t, wvfms,saveData)



%% Fourier Analysis
figure(4)
slider_freq = uicontrol('Style', 'slider', 'Min', 1, 'Max', size(wvfms, 2), 'sliderstep', [1/size(wvfms, 2),1/10] ,'Value', 1, 'Position', [0.4*figure(4).Position(3), 0.05*figure(4).Position(4), 0.2*figure(4).Position(3), 0.04*figure(4).Position(4)]);
slider_label_freq = uicontrol('Style', 'text', 'Position', [0.4*figure(4).Position(3), 0.1*figure(4).Position(4), 0.2*figure(4).Position(3), 0.04*figure(4).Position(4)], 'String', 'Frequency Slider');
slider_freq.Callback = @(src, event) updateFreqPlot(src, t, wvfms, scpSettings.SampleFrequency);
updateFreqPlot(slider_freq, t, wvfms, scpSettings.SampleFrequency)

%% 
% Bandpass filter design
Fs = scpSettings.SampleFrequency; % Sampling Frequency
F0 = 0.420*1e6; % Centre
width = 0.1*1e6;
Fpass1 = F0-width; % First Passband Frequency
Fpass2 = F0+width; % Second Passband Frequency

wvfms_biases = mean(wvfms,1);

% Apply the bandpass filter
wvfms_filtered = bandpass(wvfms, [Fpass1 Fpass2], Fs);

% Plot the filtered and unfiltered waveform
figure(6)
plot(t, wvfms(:, 1)-wvfms_biases(:, 1), t, wvfms_filtered(:, 1))
xlabel('Time [us]')
ylabel('Volts')
title('Waveform @ t = 0mins (Unfiltered and 0.4MHz Bandpass Filtered)')
legend('Unfiltered', 'Filtered')

% Add a slider to control the time index for the filtered waveform plot
src_filtered.Value = 1;
slider_filtered = uicontrol('Style', 'slider', 'Min', 1, 'Max', size(wvfms_filtered, 2), 'sliderstep', [1/size(wvfms_filtered, 2),1/10] ,'Value', 1, 'Position', [0.4*figure(6).Position(3), 0.05*figure(6).Position(4), 0.2*figure(6).Position(3), 0.04*figure(6).Position(4)]);
slider_label_filtered = uicontrol('Style', 'text', 'Position', [0.4*figure(6).Position(3), 0.1*figure(6).Position(4), 0.2*figure(6).Position(3), 0.04*figure(6).Position(4)], 'String', 'Filtered Time Slider');
slider_filtered.Callback = @(src, event) updateFilteredPlot(src, t, wvfms, wvfms_filtered, saveData,wvfms_biases);
updateFilteredPlot(src_filtered, t, wvfms, wvfms_filtered, saveData,wvfms_biases)

function updateFilteredPlot(src, t, wvfms, wvfms_filtered, saveData,wvfms_biases)
    plot(t, wvfms(:, round(src.Value))-wvfms_biases(:, round(src.Value)), 'Color', [0 0 1 0.25])
    hold on
    plot(t, wvfms_filtered(:, round(src.Value)))
    hold off
    xlabel('Time [us]')
    ylabel('Volts')
    ylim([-0.05,0.05])
    %xlim([0,25])
    time = saveData.timestamps(round(src.Value));
    title(strcat('Waveform @ t =',string(time),'s (Unfiltered and 0.419MHz Bandpass Filtered)'))
    legend('Unfiltered', 'Filtered')
end
%% Funcitons

function updatePlot(src, t, wvfms,saveData)
    plot(t, wvfms(:, round(src.Value)))
    xlabel('Time [us]')
    ylabel('Volts')
    ylim([-0.15,0.15])
    %xlim([0,25])
    time = saveData.timestamps(round(src.Value));
    title(strcat('Waveform @ t =',string(time/60),'mins'))
end

function updateFreqPlot(src, t, wvfms, SampleFrequency)
    idx = round(src.Value);
    L = length(t);
    Y = fft(wvfms(:, idx));
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = SampleFrequency*(0:(L/2))/L;
    
    plot(f*1e-6, P1)
    xlim([0,1])
    xlabel('Frequency (MHz)')
    ylabel('|P1(f)|')
    title('Single-Sided Amplitude Spectrum of X(t)')
end