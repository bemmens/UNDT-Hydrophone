clearvars

folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';
file_name = strcat('PinCyclePin',string(3));
path = strcat(folder_path,file_name,'.mat');
load(path)

bigness = size(scanData);

size(scanData)
x_index = 13;
y_index = 13;
data = squeeze(scanData(x_index,y_index,:,1));

pinNums = [1,2,3,4,6,7,8,9];
waveforms = zeros(bigness(3),length(pinNums));

folder_path = 'C:\Users\gv19838\OneDrive - University of Bristol\PhD\Hydrophone\UNDT-Hydrophone\DataOut\';

for i = 1:length(pinNums)
    pinNum = pinNums(i);
    file_name = strcat('PinCyclePin',string(pinNum));
    path = strcat(folder_path,file_name,'.mat');
    load(path)
    waveforms(:,i) = squeeze(scanData(x_index,y_index,:,:));;
end
%%
%waveformsflat = reshape(waveforms, [], bigness(3), length(pinNums));

%%
t = (1:scpSettings.RecordLength)*1e6/scpSettings.SampleFrequency; % us

figure(1)
for i = 1:length(waveforms(1,:))
    plot(t,waveforms(:,i))
    xlim([25,60])
    hold on
end
xlabel('Time [us]')
ylabel('Voltage [V]')
title('Waveform for Single Pin Activation')
legend(string(pinNums))
hold off

%%
figure(2)
variation = std(waveforms,0,2);
plot(t,variation)
xlim([25,60])
xlabel('Time [us]')
ylabel('Voltage [V]')
title('Standard Deviation of Waveforms Across All Pins')

%% Mean
figure(3)
avrg = mean(waveforms,2);
plot(t,avrg)
xlabel('Time [us]')
ylabel('Voltage [V]')
title('Mean of Waveforms Across All Pins')
xlim([25,60])

%% Band pass filter
bpass = bandpass(avrg,[0.9e6,2e6],scpSettings.SampleFrequency); % 1 +/- 0.1MHz
figure(4)
plot(t,bpass)
xlim([25,60])
xlabel('Time [us]')
ylabel('Voltage [V]')
title('Mean Waveform with Bandpass Filter: 1+/-0.1MHz')

%% Spectrogram
figure(5)
spectrogram(avrg,7e2,[],1e3,scpSettings.SampleFrequency,'yaxis');
ylim([0,5])
colorbar('off')
title('Bandpassed Waveform Spectrogram')
%xlabel('Frequency (Hz)')
%ylabel('Time (us)')

%% Listen
soundsc(bpass,10e3) % waveform,playback sample freqency