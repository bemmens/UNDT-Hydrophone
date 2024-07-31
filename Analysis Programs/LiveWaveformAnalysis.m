File_loc = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\'; % CHECK
File_name = 'RingArray2Hrs'; % CHECK
Save_String=strcat(File_loc,File_name,'.mat');
load(Save_String)

wvfms_raw = saveData.data'; 
wvfms = wvfms_raw(:, any(wvfms_raw, 1)); % remove empmty recordings
timestamps = saveData.timestamps(:,any(saveData.timestamps,1));

%%
centred = noBias(wvfms);
mag = abs(centred);
mean_mag = mean(mag);

%%
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

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
src.Value = 1;
figure(3)
slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', size(wvfms, 2), 'sliderstep', [1/size(wvfms, 2),1/10] ,'Value', 1, 'Position', [0.4*figure(3).Position(3), 0.05*figure(3).Position(4), 0.2*figure(3).Position(3), 0.04*figure(3).Position(4)]);
slider_label = uicontrol('Style', 'text', 'Position', [0.4*figure(3).Position(3), 0.1*figure(3).Position(4), 0.2*figure(3).Position(3), 0.04*figure(3).Position(4)], 'String', 'Time Slider');
slider.Callback = @(src, event) updatePlot(src, t, wvfms,saveData);
updatePlot(src, t, wvfms,saveData)

function updatePlot(src, t, wvfms,saveData)
    plot(t, wvfms(:, round(src.Value)))
    xlabel('Time [us]')
    ylabel('Volts')
    ylim([-1,1.5])
    xlim([0,25])
    time = saveData.timestamps(round(src.Value));
    title(strcat('Waveform @ t =',string(time/60),'mins'))
end